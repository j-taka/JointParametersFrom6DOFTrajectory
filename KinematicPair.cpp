/*! \brief  Estimate parameters of various joints from observation
 *           (If you know details, please see Section 5 in doctor dissertation of Jun Takamatsu) 
 *  \file   KinematicPair.cpp
 *  \author Jun Takamatsu <j-taka@cvl.iis.u-tokyo.ac.jp>
 *  \date   21 May 2004
 *  \date   17 Feb 2015 (re-implement)
 */

#include <cmath>
#include <cassert>
#include <iostream>
#include "KinematicPair.h"
#include "eigen_matrix_utility.h"
#include <gsl/gsl_multimin.h>
#include <PolynomialSolver.h>
#include <MonomialParser.h>

//
void KinematicPair::SaveBase(std::ofstream &ofs) const
{
	ofs << param.size() << std::endl;
	for (size_t i(0); i < param.size(); ++i) {
		ofs << param[i][0] << " " << param[i][1] << " " << param[i][2] << std::endl;
	}
}

void KinematicPair::LoadBase(std::ifstream &ifs)
{
	size_t num;
	ifs >> num;
	param.resize(num);
	for (size_t i(0); i < num; ++i) {
		ifs >> param[i][0] >> param[i][1] >> param[i][2];
	}
}



/*! \biref wmin in Singular Value Decomposition */
static const double SVDEPSILON = 1.0e-6;

/*! \brief  Calculate differnce between two orientations
 *          The difference can be represented \theta-radian rotation.
 *          This function returns \theta
 *  
 *  \return the difference [radian]
 *  \param  o1 one orientation
 *  \param  o2 the other orientation 
 */
static double DiffOri(const Eigen::Matrix3d &o1, const Eigen::Matrix3d &o2)
{
    const Eigen::Matrix3d temp = o2 * o1.transpose(); // diffenece orientation	
    return acos(ProperTriFunc((temp.trace() - 1.0) / 2.0)); /* 1+2 cos \theta */
}

/*! \brief calculate the correctness of the estimation in revolute joint
 * \return value of the evaluation function 
 * \param xvec_ptr [0]: azimus angle of Axis A
 *                 [1]: elevation angle of Axis A
 *                 [2]: azimus angle of Axis B
 *                 [3]: elevation angle of Axis B 
 * \param params pointer to motion data
 */
double RotKPEstimate(const ::gsl_vector *xvec_ptr, void *params)
{
	const std::vector<MotionMatrixd> *param_ptr = (std::vector<MotionMatrixd>*) params;
	 
	const double sa = sin(::gsl_vector_get(xvec_ptr, 0)); 
	const double ca = cos(::gsl_vector_get(xvec_ptr, 0));
    const double sb = sin(::gsl_vector_get(xvec_ptr, 1)); 
	const double cb = cos(::gsl_vector_get(xvec_ptr, 1));
    const double sc = sin(::gsl_vector_get(xvec_ptr, 2));
	const double cc = cos(::gsl_vector_get(xvec_ptr, 2));
    const double sd = sin(::gsl_vector_get(xvec_ptr, 3)); 
	const double cd = cos(::gsl_vector_get(xvec_ptr, 3));

	double res(0.0);
	std::vector<MotionMatrixd>::const_iterator it;
	for (it = param_ptr->begin(); it != param_ptr->end(); it++){
        const Eigen::Matrix3d r = it->R();
        /* evaluation function is \sum_{i} (1 - \cos \theta) */
        res += 1.0 - r(0, 0) * sa * cb * sc * cd - r(0, 1) * sa * sb * sc * cd - r(0, 2) * ca * sc * cd
					- r(1, 0) * sa * cb * sc * sd - r(1, 1) * sa * sb * sc * sd - r(1, 2) * ca * sc * sd
					- r(2, 0) * sa * cb * cc - r(2, 1) * sa * sb * cc - r(2, 2) * ca * cc;
    }
    return res;
}

/*! \brief calculate the gradient of the evaluation function of the revolute joint
 * \param xvec_ptr [0]: azimus angle of Axis A
 *                 [1]: elevation angle of Axis A
 *                 [2]: azimus angle of Axis B
 *                 [3]: elevation angle of Axis B 
 * \param params pointer to motion data
 * \retval df_ptr [0]: difference along azimus angle of Axis A
 *                [1]: difference along elevation angle of Axis A
 *                [2]: difference along azimus angle of Axis B
 *                [3]: difference elevation angle of Axis B 
 */
void DRotKPEstimate(const ::gsl_vector *xvec_ptr, void *params, ::gsl_vector *df_ptr)
{
	std::vector<MotionMatrixd> *param_ptr = (std::vector<MotionMatrixd> *) params;

	const double sa = sin(::gsl_vector_get(xvec_ptr, 0)); 
	const double ca = cos(::gsl_vector_get(xvec_ptr, 0));
    const double sb = sin(::gsl_vector_get(xvec_ptr, 1)); 
	const double cb = cos(::gsl_vector_get(xvec_ptr, 1));
    const double sc = sin(::gsl_vector_get(xvec_ptr, 2));
	const double cc = cos(::gsl_vector_get(xvec_ptr, 2));
    const double sd = sin(::gsl_vector_get(xvec_ptr, 3)); 
	const double cd = cos(::gsl_vector_get(xvec_ptr, 3));

	double dest[4] = { 0.0, 0.0, 0.0, 0.0 };

	std::vector<MotionMatrixd>::const_iterator it;
    for (it = param_ptr->begin(); it != param_ptr->end(); it++){
        const Eigen::Matrix3d r = it->R();
        /* evaluation function is \sum_{i} (1 - \cos \theta) */
        /* difference wrt \alpha*/
        dest[0] += -r(0, 0) * ca * cb * sc * cd - r(0, 1) * ca * sb * sc * cd + r(0, 2) * sa * sc * cd
					- r(1, 0) * ca * cb * sc * sd - r(1, 1) * ca * sb * sc * sd + r(1, 2) * sa * sc * sd
					- r(2, 0) * ca * cb * cc - r(2, 1) * ca * sb * cc + r(2, 2) * sa * cc;
        /* difference wrt \beta*/
        dest[1] += r(0, 0) * sa * sb * sc * cd - r(0, 1) * sa * cb * sc * cd
					+ r(1, 0) * sa * sb * sc * sd - r(1, 1) * sa * cb * sc * sd
					+ r(2, 0) * sa * sb * cc - r(2, 1) * sa * cb * cc;
        /* difference wrt \phi*/
        dest[2] += -r(0, 0) * sa * cb * cc * cd - r(0, 1) * sa * sb * cc * cd - r(0, 2) * ca * cc * cd
					- r(1, 0) * sa * cb * cc * sd - r(1, 1) * sa * sb * cc * sd - r(1, 2) * ca * cc * sd
					+ r(2, 0) * sa * cb * sc + r(2, 1) * sa * sb * sc + r(2, 2) * ca * sc;
        /* difference wrt \gamma */
        dest[3] += r(0, 0) * sa * cb * sc * sd + r(0, 1) * sa * sb * sc * sd + r(0, 2) * ca * sc * sd
					- r(1, 0) * sa * cb * sc * cd - r(1, 1) * sa * sb * sc * cd - r(1, 2) * ca * sc * cd;
    }
	::gsl_vector_set(df_ptr, 0, dest[0]);
	::gsl_vector_set(df_ptr, 1, dest[1]);
	::gsl_vector_set(df_ptr, 2, dest[2]);
	::gsl_vector_set(df_ptr, 3, dest[3]);

    return;
}


/*! \brief calculate the evaluation function and its gradient of the revolute joint
 * \param xvec_ptr [0]: azimus angle of Axis A
 *                 [1]: elevation angle of Axis A
 *                 [2]: azimus angle of Axis B
 *                 [3]: elevation angle of Axis B 
 * \param params pointer to motion data
 * \retval f_prt value of the evaluation function
 * \retval df_ptr [0]: difference along azimus angle of Axis A
 *                [1]: difference along elevation angle of Axis A
 *                [2]: difference along azimus angle of Axis B
 *                [3]: difference elevation angle of Axis B 
 */
void FDRotKPEstimate(const ::gsl_vector *xvec_ptr, void *params_ptr, double *f_ptr, ::gsl_vector *df_ptr)
{
	*f_ptr = RotKPEstimate(xvec_ptr, params_ptr);
	DRotKPEstimate(xvec_ptr, params_ptr, df_ptr);
}

struct StKPParam
{
	Eigen::Matrix3d baseR;
	const std::vector<MotionMatrixd> *motion;
};

/*! \brief calculate the correctness of the estimation in fixed and prismatic joint
 * \return value of the evaluation function 
 * \param xvec_ptr orientation (RPY)
 * \param params pointer to motion data
 */
double StKPEstimate(const ::gsl_vector *xvec_ptr, void *params)
{
	const Eigen::Matrix3d baseR = ((StKPParam*)params)->baseR;
	const std::vector<MotionMatrixd> *param_ptr = ((StKPParam*)params)->motion;
    const Eigen::Matrix3d ro = RPY2Mat(Eigen::Vector3d(::gsl_vector_get(xvec_ptr, 0), ::gsl_vector_get(xvec_ptr, 1), ::gsl_vector_get(xvec_ptr, 2))) * baseR;
 
	double res(0.0);
	std::vector<MotionMatrixd>::const_iterator it;
    for (it = param_ptr->begin(); it != param_ptr->end(); it++){
        const Eigen::Matrix3d r = it->R();
        /* evaluation function is \sum_{i} (1 - \cos \theta) */
		res += 0.5 * (3 - (r * ro.transpose()).trace());
		/*
		res += 0.5 * (3 - r(0, 0) * ro(0, 0) - r(0, 1) * ro(0, 1) - r(0, 2) * ro(0, 2)
						- r(1, 0) * ro(1, 0) - r(1, 1) * ro(1, 1) - r(1, 2) * ro(1, 2)
						- r(2, 0) * ro(2, 0) - r(2, 1) * ro(2, 1) - r(2, 2) * ro(2, 2));
						*/
	}
    return res;
}

/* 変数が与えられたときの評価関数のグラディエントを求める
 * param src[0] : ロール
 * param src[1] : ピッチ
 * param src[2] : ヨー
 * param dest[0] : ロールの偏微分
 * param dest[1] : ピッチの偏微分
 * param dest[2] : ヨーの偏微分 */
void DStKPEstimate(const ::gsl_vector *xvec_ptr, void *params, ::gsl_vector *df_ptr)
{
	const Eigen::Matrix3d baseR = ((StKPParam*)params)->baseR;
	const std::vector<MotionMatrixd> *param_ptr = ((StKPParam*)params)->motion;
	
	const Eigen::Vector3d temp(::gsl_vector_get(xvec_ptr, 0), ::gsl_vector_get(xvec_ptr, 1), ::gsl_vector_get(xvec_ptr, 2));
    const Eigen::Matrix3d dr = RPYdR2Mat(temp);
    const Eigen::Matrix3d dp = RPYdP2Mat(temp);
    const Eigen::Matrix3d dy = RPYdY2Mat(temp);

	double dest[3] = { 0.0, 0.0, 0.0 };

	std::vector<MotionMatrixd>::const_iterator it;
    for (it = param_ptr->begin();it != param_ptr->end(); it++){
        /* Evaluatinon function is \sum_{i} (1 - \cos \theta) */
        const Eigen::Matrix3d r = it->R();
		dest[0] += 0.5 * (3 - (r * (dr * baseR).transpose()).trace());
		dest[1] += 0.5 * (3 - (r * (dp * baseR).transpose()).trace());
		dest[2] += 0.5 * (3 - (r * (dy * baseR).transpose()).trace());
		/*
		dest[0] += 0.5 * (3 - r(0, 0) * dr(0, 0) - r(0, 1) * dr(0, 1) - r(0, 2) * dr(0, 2)
							- r(1, 0) * dr(1, 0) - r(1, 1) * dr(1, 1) - r(1, 2) * dr(1, 2)
							- r(2, 0) * dr(2, 0) - r(2, 1) * dr(2, 1) - r(2, 2) * dr(2, 2));
		dest[1] += 0.5 * (3 - r(0, 0) * dp(0, 0) - r(0, 1) * dp(0, 1) - r(0, 2) * dp(0, 2)
						- r(1, 0) * dp(1, 0) - r(1, 1) * dp(1, 1) - r(1, 2) * dp(1, 2)
							- r(2, 0) * dp(2, 0) - r(2, 1) * dp(2, 1) - r(2, 2) * dp(2, 2));
		dest[2] += 0.5 * (3 - r(0, 0) * dy(0, 0) - r(0, 1) * dy(0, 1) - r(0, 2) * dy(0, 2)
							- r(1, 0) * dy(1, 0) - r(1, 1) * dy(1, 1) - r(1, 2) * dy(1, 2)
							- r(2, 0) * dy(2, 0) - r(2, 1) * dy(2, 1) - r(2, 2) * dy(2, 2));
		*/
	}
	::gsl_vector_set(df_ptr, 0, dest[0]);
	::gsl_vector_set(df_ptr, 1, dest[1]);
	::gsl_vector_set(df_ptr, 2, dest[2]);

    return;
}

void FDStKPEstimate(const ::gsl_vector *xvec_ptr, void *params_ptr, double *f_ptr, ::gsl_vector *df_ptr)
{
	*f_ptr = StKPEstimate(xvec_ptr, params_ptr);
	DStKPEstimate(xvec_ptr, params_ptr, df_ptr);
}

Eigen::Matrix3d KinematicPair::RotXY(double src1, double src2)
{
	return Eigen::AngleAxisd(src1, Eigen::Vector3d::UnitX()).matrix() * Eigen::AngleAxisd(src2, Eigen::Vector3d::UnitY()).matrix();
}

Eigen::Matrix3d KinematicPair::RotXZ(double src1, double src2)
{
	return Eigen::AngleAxisd(src1, Eigen::Vector3d::UnitX()).matrix() * Eigen::AngleAxisd(src2, Eigen::Vector3d::UnitZ()).matrix();
}

/* 変数が与えられたときの評価関数のグラディエントを求める
 * 物体Aとリンク座標間での  Rx Rz
 * 物体Bとリンク座標間での  Rx Ry */
double SRotKPEstimate(const ::gsl_vector *xvec_ptr, void *params)
{
	const std::vector<MotionMatrixd> *param_ptr = (std::vector<MotionMatrixd>*) params;
	const Eigen::Vector3d y(0, 1, 0), z(0, 0, 1);
	const Eigen::Matrix3d al = KinematicPair::RotXZ(::gsl_vector_get(xvec_ptr, 0), ::gsl_vector_get(xvec_ptr, 1));
	const Eigen::Matrix3d bl = KinematicPair::RotXY(::gsl_vector_get(xvec_ptr, 2), ::gsl_vector_get(xvec_ptr, 3));
	double res(0.0);
	std::vector<MotionMatrixd>::const_iterator it;
    for (it=param_ptr->begin(); it!=param_ptr->end(); it++){
        const Eigen::Matrix3d r = it->R();
        /* evaluation function is \sum_{i} (\sin \theta)^2 */
        res += (bl * z).dot(r * al * y) * (bl * z).dot(r * al * y);
    }
    return res;
}

/* 変数が与えられたときの評価関数のグラディエントを求める
 * 物体Aとリンク座標間での Rx Rz
 * 物体Bとリンク座標間での Rx Ry */
void DSRotKPEstimate(const ::gsl_vector *xvec_ptr, void *params, ::gsl_vector *df_ptr)
{
 	const std::vector<MotionMatrixd> *param_ptr = (std::vector<MotionMatrixd>*) params;
	const Eigen::Vector3d y(0, 1, 0), z(0, 0, 1);
	const Eigen::Matrix3d al   = KinematicPair::RotXZ(::gsl_vector_get(xvec_ptr, 0), ::gsl_vector_get(xvec_ptr, 1));
	const Eigen::Matrix3d d1al = dRotX2Mat(::gsl_vector_get(xvec_ptr, 0)) * Eigen::AngleAxisd(::gsl_vector_get(xvec_ptr, 1), Eigen::Vector3d::UnitZ()).matrix();
	const Eigen::Matrix3d d2al = Eigen::AngleAxisd(::gsl_vector_get(xvec_ptr, 0), Eigen::Vector3d::UnitX()).matrix() * dRotZ2Mat(::gsl_vector_get(xvec_ptr, 1));

	const Eigen::Matrix3d bl   = KinematicPair::RotXY(::gsl_vector_get(xvec_ptr, 2), ::gsl_vector_get(xvec_ptr, 3));
	const Eigen::Matrix3d d3bl = dRotX2Mat(::gsl_vector_get(xvec_ptr, 2)) * Eigen::AngleAxisd(::gsl_vector_get(xvec_ptr, 3), Eigen::Vector3d::UnitY()).matrix();
	const Eigen::Matrix3d d4bl = Eigen::AngleAxisd(::gsl_vector_get(xvec_ptr, 2), Eigen::Vector3d::UnitX()).matrix() * dRotY2Mat(::gsl_vector_get(xvec_ptr, 3));

	double dest[4] = { 0.0, 0.0, 0.0, 0.0 };
    
	std::vector<MotionMatrixd >::const_iterator it;
	for (it = param_ptr->begin(); it != param_ptr->end(); it++){
		const Eigen::Matrix3d r = it->R();
        /* evaluation function is \sum_{i} (\sin \theta)^{2} */
        dest[0] += 2 * (bl * z).dot(r * d1al * y) * (bl * z).dot(r * al * y);
        dest[1] += 2 * (bl * z).dot(r * d2al * y) * (bl * z).dot(r * al * y);
        dest[2] += 2 * (d3bl * z).dot(r * al * y) * (bl * z).dot(r * al * y);
        dest[3] += 2 * (d4bl * z).dot(r * al * y) * (bl * z).dot(r * al * y);
    }
	::gsl_vector_set(df_ptr, 0, dest[0]);
	::gsl_vector_set(df_ptr, 1, dest[1]);
	::gsl_vector_set(df_ptr, 2, dest[2]);
	::gsl_vector_set(df_ptr, 3, dest[3]);

    return;
}

void FDSRotKPEstimate(const ::gsl_vector *xvec_ptr, void *params_ptr, double *f_ptr, ::gsl_vector *df_ptr)
{
	*f_ptr = SRotKPEstimate(xvec_ptr, params_ptr);
	DSRotKPEstimate(xvec_ptr, params_ptr, df_ptr);
}

std::ostream& operator<<(std::ostream &o, const EstError &src)
{
    o << "Translation error: " << src.transDiv << std::endl
      << "Rotation error: " << src.rotDiv;
    return o;
}

void KinematicPair::EstOri0Param(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd> &src, double &dest2)
{
    assert(!src.empty());
	/* initialize for gsl */
	::gsl_multimin_function_fdf my_func;
	my_func.n = 3;
	my_func.f = &StKPEstimate;
	my_func.df = &DStKPEstimate;
	my_func.fdf = &FDStKPEstimate;
	StKPParam stkp_param;
	stkp_param.baseR = src.front().R();
	stkp_param.motion = &src;
	my_func.params = (void*) &stkp_param;

	/* set initial guess */
	::gsl_vector *vec_p = ::gsl_vector_alloc(3);
	::gsl_vector_set(vec_p, 0, 0);
	::gsl_vector_set(vec_p, 1, 0);
	::gsl_vector_set(vec_p, 2, 0);
	
	const ::gsl_multimin_fdfminimizer_type *type_ptr = ::gsl_multimin_fdfminimizer_conjugate_fr;
	::gsl_multimin_fdfminimizer *minimizer_ptr = ::gsl_multimin_fdfminimizer_alloc(type_ptr, my_func.n);

	// set the tolerance and starting step size
	double step_size = 0.01;
	double tolerance = 1.0e-4;
	::gsl_multimin_fdfminimizer_set(minimizer_ptr, &my_func, vec_p, step_size, tolerance);

	size_t iteration(0);
	const size_t max_iteration(100);

	int status(0);
	do{
		iteration++;
		status = ::gsl_multimin_fdfminimizer_iterate(minimizer_ptr);
		
		if (status){ 
			break;
		}

		status = ::gsl_multimin_test_gradient(minimizer_ptr->gradient, tolerance);
	} while (status == ::GSL_CONTINUE && iteration < max_iteration);
	// set answer
	const Eigen::Vector3d vec(::gsl_vector_get(minimizer_ptr->x, 0), ::gsl_vector_get(minimizer_ptr->x, 1), ::gsl_vector_get(minimizer_ptr->x, 2));
	const Eigen::Matrix3d ro = RPY2Mat(vec) * stkp_param.baseR;
	Eigen::Vector3d vec2 = Mat2RPY(ro);

    param.push_back(vec2);
    RefineOri0(dest, src, dest2);
	// free 
	::gsl_multimin_fdfminimizer_free(minimizer_ptr);
	::gsl_vector_free(vec_p);
}

/* 1軸回転の軸向きパラメータ推定
 * param dest  : 物体aからみた物体bの姿勢修正後の位置姿勢の配列
 * param src   : 物体aからみた物体bのオリジナルの位置姿勢の配列
 * param dest2 : 誤差の平均 */
void KinematicPair::EstOri1Param(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd> &src, double &dest2)
{
	assert(!src.empty());
	/* initialize for gsl */
	::gsl_multimin_function_fdf my_func;
	my_func.n = 4;
	my_func.f = &RotKPEstimate;
	my_func.df = &DRotKPEstimate;
	my_func.fdf = &FDRotKPEstimate;
	my_func.params = (void*) &src;

	/* set initial guess */
	const size_t _NUM_TRIAL = 10;
	int mint1(0), mint2((int) src.size() - 1);
    double minang = fabs((src[mint1].R().transpose() * src[mint2].R()).trace() - 1.0);
    for (size_t i(0); i < _NUM_TRIAL; ++i){
        // select two configurations
        int t1 = rand() % src.size();
        int t2 = rand() % src.size();
        // calculate angle
        const double ang = fabs((src[t1].R().transpose() * src[t2].R()).trace() - 1.0);
        if (ang < minang){
            minang = ang; 
			mint1 = t1; 
			mint2 = t2;
        }
    }
    Eigen::Vector3d la;
	double ang;
	Mat2AngleAxis(ang, la, Eigen::Matrix3d(src[mint1].R().transpose() * src[mint2].R()));
	::gsl_vector *vec_p = ::gsl_vector_alloc(4);
    ::gsl_vector_set(vec_p, 0, acos(ProperTriFunc(la[2])));
    ::gsl_vector_set(vec_p, 1, atan2(la[1], la[0]));
    const Eigen::Vector3d lb = src[mint1].R() * la;
    ::gsl_vector_set(vec_p, 2, acos(ProperTriFunc(lb[2])));
    ::gsl_vector_set(vec_p, 3, atan2(lb[1], lb[0]));

	const ::gsl_multimin_fdfminimizer_type *type_ptr = ::gsl_multimin_fdfminimizer_conjugate_fr;
	::gsl_multimin_fdfminimizer *minimizer_ptr = ::gsl_multimin_fdfminimizer_alloc(type_ptr, my_func.n);

	// set the tolerance and starting step size
	double step_size = 0.01;
	double tolerance = 1.0e-4;
	::gsl_multimin_fdfminimizer_set(minimizer_ptr, &my_func, vec_p, step_size, tolerance);

	size_t iteration(0);
	const size_t max_iteration(100);

	int status(0);
	do{
		iteration++;
		status = ::gsl_multimin_fdfminimizer_iterate(minimizer_ptr);
		
		if (status){ 
			break;
		}

		status = ::gsl_multimin_test_gradient(minimizer_ptr->gradient, tolerance);
	} while (status == ::GSL_CONTINUE && iteration < max_iteration);
	// set answer
	Eigen::Vector3d vec;

	// set answer
    vec[0] = sin(::gsl_vector_get(minimizer_ptr->x, 0)) * cos(::gsl_vector_get(minimizer_ptr->x, 1));
    vec[1] = sin(::gsl_vector_get(minimizer_ptr->x, 0)) * sin(::gsl_vector_get(minimizer_ptr->x, 1));
    vec[2] = cos(::gsl_vector_get(minimizer_ptr->x, 0));
    vec /= vec.norm(); // normalize
    param.push_back(vec);
    vec[0] = sin(::gsl_vector_get(minimizer_ptr->x, 2)) * cos(::gsl_vector_get(minimizer_ptr->x, 3));
    vec[1] = sin(::gsl_vector_get(minimizer_ptr->x, 2)) * sin(::gsl_vector_get(minimizer_ptr->x, 3));
    vec[2] = cos(::gsl_vector_get(minimizer_ptr->x, 2));
    vec /= vec.norm(); // normalize
    param.push_back(vec);
    RefineOri1(dest, src, dest2); // refine data

	// free 
	::gsl_multimin_fdfminimizer_free(minimizer_ptr);
	::gsl_vector_free(vec_p);
}

bool KinematicPair::isValidSolutionInOri2(const std::vector<std::complex<double> > &src) const
{
	for (size_t i(0); i <src.size(); ++i){
		if (fabs(src[i].imag()) > _NEARLY_ZERO){
			return false;
		}
	}
	return true;
}

bool KinematicPair::CalculateAnotherAxis(Eigen::Vector3d &dest, const std::vector<Eigen::Vector3d> &src) const
{
	Eigen::Matrix3d mat = Eigen::Matrix3d::Zero();
	for (size_t i(0); i < src.size(); ++i){
		mat += src[i] * src[i].transpose();
	}
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(mat);
	double max_eig = eig.eigenvalues()[0];
	for (int i(1); i < 3; ++i){
		if (max_eig < eig.eigenvalues()[i]){
			max_eig = eig.eigenvalues()[i];
		}
	}
	double min_eig = max_eig * _NEARLY_ZERO;
	bool set(false);
	for (int i(0); i < 3; ++i){
		if (eig.eigenvalues()[i] < min_eig){
			if (set){ return false; }
			dest = eig.eigenvectors().col(i).transpose();
			set = true;
		}
	}
	return set;
}

void KinematicPair::SetEquationInInitOri2(std::vector<Polynomial<double> > &I, const std::vector<size_t> &samp, const std::vector<MotionMatrixd> &src) const
{
	std::vector<std::string> names(3);
	names[0] = std::string("a_x");
	names[1] = std::string("a_y");
	names[2] = std::string("a_z");
	MonomialParser mp(names);
	I.resize(3);
	I[0] = 1.0 * mp("a_x^2") + 1.0 * mp("a_y^2") + 1.0 * mp("a_z^2") - 1.0 * mp("");
	Polynomial<double> vec[3] = { 1.0 * mp("a_x"), 1.0 * mp("a_y"), 1.0 * mp("a_z") };
	for (size_t i(1); i < 3; i++){
		const Eigen::Matrix3d R1 = src[samp[i]].R() * src[samp[0]].R().transpose();
		const Eigen::Matrix3d R2 = src[samp[i + 1]].R() * src[samp[0]].R().transpose();
		Polynomial<double> vec1[3] = { R1(0, 0) * mp("a_x") + R1(0, 1) * mp("a_y") + R1(0, 2) * mp("a_z"),  
									   R1(1, 0) * mp("a_x") + R1(1, 1) * mp("a_y") + R1(1, 2) * mp("a_z"), 
									   R1(2, 0) * mp("a_x") + R1(2, 1) * mp("a_y") + R1(2, 2) * mp("a_z") }; 
		Polynomial<double> vec2[3] = { R2(0, 0) * mp("a_x") + R2(0, 1) * mp("a_y") + R2(0, 2) * mp("a_z"),  
									   R2(1, 0) * mp("a_x") + R2(1, 1) * mp("a_y") + R2(1, 2) * mp("a_z"), 
									   R2(2, 0) * mp("a_x") + R2(2, 1) * mp("a_y") + R2(2, 2) * mp("a_z") }; 
		I[i] = ::RemoveSmallValue(vec2[0] * (vec[1] * vec1[2] - vec[2] * vec1[1]) 
								+ vec2[1] * (vec[2] * vec1[0] - vec[0] * vec1[2])
								+ vec2[2] * (vec[0] * vec1[1] - vec[1] * vec1[0]), _NEARLY_ZERO);
	}
}

void KinematicPair::SetSolutionInOri2(std::vector<double> &dest, const Eigen::Vector3d &la, const Eigen::Vector3d &lb) const
{
	dest.resize(4);
	if (fabs(la[1]) < _NEARLY_ZERO && fabs(la[2]) < _NEARLY_ZERO){
		dest[0] = 0.0;
		dest[1] = atan2(-la[0], 0);
	}
	else{
		dest[0] = atan2(la[2], la[1]);
		const double cz = sqrt(pow(la[1], 2) + pow(la[2], 2)); 
		double sign;
		if (fabs(la[1]) > fabs(la[2])){
			sign = (la[1] * cos(dest[0]) > 0 ? 1 : -1);
		}
		else{
			sign = (la[2] * sin(dest[0]) > 0 ? 1 : -1);
		}
		dest[1] = atan2(-la[0], sign * cz);
	}
	if (fabs(lb[1]) < _NEARLY_ZERO && fabs(lb[2]) < _NEARLY_ZERO){
		dest[2] = 0.0;
		dest[3] = atan2(lb[0], 0);
	}
	else{
		dest[2] = atan2(-lb[1], lb[2]);
		const double cy = sqrt(pow(lb[1], 2) + pow(lb[2], 2));
		double sign;
		if (fabs(lb[1]) > fabs(lb[2])){
			sign = (lb[1] * -sin(dest[2]) > 0 ? 1 : -1);
		}
		else{
			sign = (lb[2] * cos(dest[2]) > 0 ? 1 : -1);
		}
		dest[3] = atan2(lb[0], sign * cy);
	}
}

void KinematicPair::ChooseRandom(std::vector<size_t> &dest, size_t max_val)
{
	assert(dest.size() <= max_val);
	for (size_t i(0); i < dest.size(); i++){
		size_t r = (size_t)((rand() / ((double) RAND_MAX + 1.0)) * max_val);
		size_t j;
		for (j = 0; j < i; ++j){
			if (r >= dest[j]){
				r++;
			}
			else{
				for (size_t k(0); k < i - j; ++k){
					dest[i - k] = dest[i - k - 1];
				}
				break;
			}
		}
		dest[j] = r;
		max_val--;
	}
}

void KinematicPair::InitOri2Param(std::vector<double> &dest, const std::vector<MotionMatrixd> &src) const
{
	std::vector<std::string> names(3);
	names[0] = std::string("a_x");
	names[1] = std::string("a_y");
	names[2] = std::string("a_z");

	const int MAX_LOOP = 10;
	bool isfirst(true);
	double eval(0.0);
	for (int l(0); l < MAX_LOOP; ++l){
		// choose 4 randomly
		std::vector<size_t> samp(4);
		ChooseRandom(samp, src.size());
		std::vector<Polynomial<double> > I;
		SetEquationInInitOri2(I, samp, src);
		for (size_t i(0); i < I.size(); ++i){
			I[i].SetOrder(Polynomial<double>::_GREVLEX);
		}
		std::vector<Polynomial<double> > G, Gr;
		PolynomialSolver<double> ps;
		ps.BuchBergerAlgorithm(G, I);
		ps.Reduced(Gr, G);
		if (ps.isConstantIncluded(Gr)){
			continue; // ...
		}
		ps.SetEquation(Gr);
		std::vector<PolynomialSolver<double>::Answer> ans;
		ps.Solve(ans);
		for (size_t i(0); i < ans.size(); ++i){
	#if 0
			for (size_t j(0); j < names.size(); ++j){
				std::cout << names[j] << " = " << ans[i][j] << std::endl;
			}
			for (size_t j(0); j < I.size(); ++j){
				std::cout << I[j].toString(names) << " = " << I[j].Apply(ans[i]) << std::endl;
			}
			std::cout << std::endl;
	#endif
			if (isValidSolutionInOri2(ans[i])){
				std::vector<Eigen::Vector3d> orths(4);
				orths[0] = Eigen::Vector3d(ans[i][0].real(), ans[i][1].real(), ans[i][2].real());
				for (size_t j(0); j < 3; ++j){
					const Eigen::Matrix3d R = src[samp[j]].R() * src[samp[0]].R().transpose();
					orths[j + 1] = R * orths[0];
				}
				Eigen::Vector3d lb;
				if (CalculateAnotherAxis(lb, orths)){
					// convert
					Eigen::Vector3d la = src[samp[0]].R() * orths[0];
					std::vector<double> tmp_ans;
					SetSolutionInOri2(tmp_ans, la, lb);
					// choose the best
					::gsl_vector *vec_p = ::gsl_vector_alloc(4);
					::gsl_vector_set(vec_p, 0, tmp_ans[0]);
					::gsl_vector_set(vec_p, 1, tmp_ans[1]);
					::gsl_vector_set(vec_p, 2, tmp_ans[2]);
					::gsl_vector_set(vec_p, 3, tmp_ans[3]);
					const double tmp_eval = ::SRotKPEstimate(vec_p, (void*) &src);
					if (isfirst || eval > tmp_eval){
						eval = tmp_eval;
						dest = tmp_ans;
						isfirst = false;
					}
				}
			}
		}
	}
	if (isfirst){
		std::cerr << "Fail to initialize" << std::endl;
		exit(-1);
	}
}

/* 2軸回転の軸向きパラメータ推定
 * param dest  : 物体aからみた物体bの姿勢修正後の位置姿勢の配列
 * param src   : 物体aからみた物体bのオリジナルの位置姿勢の配列
 * param dest2 : 誤差の平均 */
void KinematicPair::EstOri2Param(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd> &src, double &dest2)
{
 	assert(!src.empty());
	/* initialize for gsl */
	::gsl_multimin_function_fdf my_func;
	my_func.n = 4;
	my_func.f = &SRotKPEstimate;
	my_func.df = &DSRotKPEstimate;
	my_func.fdf = &FDSRotKPEstimate;
	my_func.params = (void*) &src;

	// set initial guess
	std::vector<double> init_sol;
	InitOri2Param(init_sol, src);
	::gsl_vector *vec_p = ::gsl_vector_alloc(4);
	::gsl_vector_set(vec_p, 0, init_sol[0]);
	::gsl_vector_set(vec_p, 1, init_sol[1]);
	::gsl_vector_set(vec_p, 2, init_sol[2]);
	::gsl_vector_set(vec_p, 3, init_sol[3]);

	const ::gsl_multimin_fdfminimizer_type *type_ptr = ::gsl_multimin_fdfminimizer_conjugate_fr;
	::gsl_multimin_fdfminimizer *minimizer_ptr = ::gsl_multimin_fdfminimizer_alloc(type_ptr, my_func.n);

	// set the tolerance and starting step size
	double step_size = 1.0e-6;
	double tolerance = 1.0e-8;
	::gsl_multimin_fdfminimizer_set(minimizer_ptr, &my_func, vec_p, step_size, tolerance);

	size_t iteration(0);
	const size_t max_iteration(100);

	int status(0);
	do{
		iteration++;
		status = ::gsl_multimin_fdfminimizer_iterate(minimizer_ptr);
		
		if (status){ 
			break;
		}

		status = ::gsl_multimin_test_gradient(minimizer_ptr->gradient, tolerance);
	} while (status == ::GSL_CONTINUE && iteration < max_iteration);
	// set answer
	const Eigen::Vector3d vec1(::gsl_vector_get(minimizer_ptr->x, 0), ::gsl_vector_get(minimizer_ptr->x, 1), 0.0); 
	param.push_back(vec1);
	const Eigen::Vector3d vec2(::gsl_vector_get(minimizer_ptr->x, 2), ::gsl_vector_get(minimizer_ptr->x, 3), 0.0); 
	param.push_back(vec2);
    RefineOri2(dest, src, dest2); // refine data

	// free 
	::gsl_multimin_fdfminimizer_free(minimizer_ptr);
	::gsl_vector_free(vec_p);
}

void KinematicPair::EstOri3Param(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd> &src, double &dest2)
{
    // nothing
    RefineOri3(dest, src, dest2);
}

/* 主成分分析により直進方向を求める
 * param dest  : 物体aからみた物体bの誤差修正後の位置姿勢の配列
 * param src   : 物体aからみた物体bの姿勢修正後の位置姿勢の配列
 * param dest2 : 誤差の平均 */
void KinematicPair::EstOri0Loc1Param(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd> &src, double &dest2)
{
	// calculate center
	Eigen::Vector3d cen = Eigen::Vector3d::Zero();
	std::vector<MotionMatrixd>::const_iterator it;
	for (it = src.begin(); it!=src.end(); it++){ 
		cen += it->T(); 
	}
    cen /= static_cast<double>(src.size());
    /* calculate covariance matrix */
	Eigen::Matrix3d a = Eigen::Matrix3d::Zero();
    for (it = src.begin(); it != src.end(); it++){
        const Eigen::Vector3d temp = it->T() - cen;
        a(0, 0) += temp[0] * temp[0];
        a(0, 1) += temp[0] * temp[1];
        a(0, 2) += temp[0] * temp[2];
        a(1, 1) += temp[1] * temp[1];
        a(1, 2) += temp[1] * temp[2];
        a(2, 2) += temp[2] * temp[2];
    }
    a(1, 0) = a(0, 1); a(2, 0) = a(0, 2); a(2, 1) = a(1, 2);
    // calculate eigen vector
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(a);
	double dmax(eig.eigenvalues()[0]);
	int check(0);
	for (int i(1);i < 3; ++i){
		if (dmax < eig.eigenvalues()[i]){
			dmax = eig.eigenvalues()[i]; 
			check = i;
		}
    }
	Eigen::Vector3d temp = eig.eigenvectors().col(check).transpose();
	temp /= temp.norm();
    Eigen::Vector3d temp2 = src.front().R().transpose() * temp;
    temp2 /= temp2.norm();
    param.push_back(temp2);
    param.push_back(temp);
    RefineOri0Loc1(dest, src, dest2); /* 推定値を用いて修正する */
}

template <typename t_matrix>
static t_matrix PseudoInverse(const t_matrix& m, const double &tolerance=1.e-6)
{
	using namespace Eigen;
	typedef JacobiSVD<t_matrix> TSVD;
	unsigned int svd_opt(ComputeThinU | ComputeThinV);
	if(m.RowsAtCompileTime!=Dynamic || m.ColsAtCompileTime!=Dynamic)
	svd_opt= ComputeFullU | ComputeFullV;
	TSVD svd(m, svd_opt);
	const typename TSVD::SingularValuesType &sigma(svd.singularValues());
	typename TSVD::SingularValuesType sigma_inv(sigma.size());
	for(long i=0; i<sigma.size(); ++i)
	{
		if(sigma(i) > tolerance)
			sigma_inv(i)= 1.0/sigma(i);
		else
			sigma_inv(i)= 0.0;
	}
	return svd.matrixV()*sigma_inv.asDiagonal()*svd.matrixU().transpose();
}
/* 1軸or2軸回転機械リンクモデルの軸位置パラメータ推定
 * param src    : 物体aからみた物体bの位置姿勢の配列
 *                ただし, 姿勢に関しては推定値であること
 * param dest2  : 誤差の平均 */
void KinematicPair::EstLoc0Param(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd> &src, double &dest2)
{
	Eigen::MatrixXd a = Eigen::MatrixXd::Zero(6, 6); 
	Eigen::VectorXd b = Eigen::VectorXd::Zero(6);
	// solve by least square minimization	
    std::vector<MotionMatrixd >::const_iterator it;
    for (it = src.begin();it != src.end(); it++){
		Eigen::MatrixXd mat = Eigen::MatrixXd(3, 6);
		mat.block(0, 0, 3, 3) = -it->R();
		mat.block(0, 3, 3, 3) = Eigen::MatrixXd::Identity(3, 3);
		a += mat.transpose() * mat;
		b += mat.transpose() * it->T();
    }
	// solve equation using svd
	Eigen::VectorXd ans = ::PseudoInverse(a) * b;
	Eigen::Vector3d temp1(ans[0], ans[1], ans[2]);
    param.push_back(temp1);
	Eigen::Vector3d temp2(ans[3], ans[4], ans[5]);
	param.push_back(temp2);
    RefineLoc0(dest, src, dest2); // refinement
}

void KinematicPair::EstOri1Hel1Param(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd> &src, double &dest2)
{
	std::vector<double> thetas;
	// estimate amount of angles
	EstOriAmount(thetas, src);

	const Eigen::Vector3d axis = RPY2Mat(param.back()) * Eigen::Vector3d::UnitZ();
	Eigen::MatrixXd a = Eigen::MatrixXd::Zero(7, 7);
	Eigen::VectorXd b = Eigen::VectorXd::Zero(7);
	// solve least square minimization
    for (size_t i(0); i < src.size(); ++i){
        Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(3, 7);
		mat.block(0, 0, 3, 3) = -src[i].R();
		mat.block(3, 0, 3, 3) = Eigen::MatrixXd::Identity(3, 3);
		for (int j(0); j < 3; ++j){ mat(j, 6) = axis[j] * thetas[i]; }
		a += mat.transpose() * mat;
		b += mat.transpose() * src[i].T();
    }
	// solve equation using svd
	Eigen::VectorXd ans = ::PseudoInverse(a) * b;
	Eigen::Vector3d temp1(ans[0], ans[1], ans[2]);
    param.push_back(temp1);
	Eigen::Vector3d temp2(ans[3], ans[4], ans[5]);
    param.push_back(temp2);
	Eigen::Vector3d temp3(ans[6], 0.0, 0.0);
    param.push_back(temp3);
	// refinement
	RefineOri1Hel1(dest, src, thetas, dest2);
}

void KinematicPair::EstLoc3Param(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd> &src, double &dest2)
{
    // nothing 
    RefineLoc3(dest, src, dest2);
}

/* \pi 以上回転してたら\pm 2 \pi することにより補正する。*/
void KinematicPair::EstOriAmount(std::vector<double> &dest, const std::vector<MotionMatrixd> &src) const
{
    const Eigen::Matrix3d r_al = RPY2Mat(param[param.size() - 2]);
	const Eigen::Matrix3d r_lb = RPY2Mat(param.back()).transpose(); 

	std::vector<MotionMatrixd >::const_iterator it;
    for (it = src.begin(); it != src.end(); it++){
        const Eigen::Matrix3d r_ba = it->R();
        const Eigen::Matrix3d r = r_lb * r_ba * r_al;
        Eigen::Vector3d temp;
		double theta;
		Mat2AngleAxis(theta, temp, r);
		if (temp.dot(Eigen::Vector3d::UnitZ()) < 0){ 
			theta =- theta;  // change sign
		}
        dest.push_back(theta);
    }
    for (size_t i(1); i < dest.size(); i++){
        const double dt = fabs(dest[i] - dest[i - 1]);
        if (dt > M_PI){
            const int k = (int)((dt + M_PI) / (2 * M_PI));
			if (dest[i] - dest[i - 1] < 0){ 
				dest[i] += 2 * k * M_PI; 
			}
			else{ 
				dest[i] -= 2 * k * M_PI; 
			}
        }
    }
}

// eigenを使ったバージョンに変更
void KinematicPair::ZAxis2Coordinate(Eigen::Matrix3d &dest, const Eigen::Vector3d &src)
{
	// covariance matrix
	Eigen::Matrix3d u = src * src.transpose();
	// eigen
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(u);
	double wmax(0.0);
	for (int i(0); i < 3; ++i){ 
		if (eig.eigenvalues()[i] > wmax){ 
			wmax = eig.eigenvalues()[i]; 
		}
	}
    const double wmin = wmax * SVDEPSILON;
    // row element of the matrix v, corresponding to eigenvalue = 0
    Eigen::Vector3d xyz[3];
	int count(0);
	for (int i(0); i < 3; ++i){
		if (eig.eigenvalues()[i] < wmin){ // eigenvalue == 0 
            assert(count<2);
			xyz[count] = eig.eigenvectors().col(i);
            xyz[count] /= xyz[count].norm();
            count++;
        }
		else{
			xyz[2] = eig.eigenvectors().col(i);
			xyz[2] /= xyz[2].norm();
		}
	}
    assert(count == 2); // double check
    if (xyz[2].dot(xyz[0].cross(xyz[1])) > 0){
        dest(0, 0) = xyz[0][0]; dest(1, 0) = xyz[0][1]; dest(2, 0) = xyz[0][2];
        dest(0, 1) = xyz[1][0]; dest(1, 1) = xyz[1][1]; dest(2, 1) = xyz[1][2];
        dest(0, 2) = xyz[2][0]; dest(1, 2) = xyz[2][1]; dest(2, 2) = xyz[2][2];
    }
    else{
        dest(0, 0) = xyz[0][0]; dest(1, 0) = xyz[0][1]; dest(2, 0) = xyz[0][2];
        dest(0, 1) = -xyz[1][0]; dest(1, 1) = -xyz[1][1]; dest(2, 1) = -xyz[1][2];
        dest(0, 2) = xyz[2][0]; dest(1, 2) = xyz[2][1]; dest(2, 2) = xyz[2][2];
    }
}

/* z軸だけが指定された方向を向いた座標系を作り,
 * その変位をロールピッチヨー表現にして返す */
Eigen::Vector3d KinematicPair::SetCoordinate(const Eigen::Vector3d &src)
{
    Eigen::Matrix3d r;
    ZAxis2Coordinate(r, src);
    return Mat2RPY(r);
}
    
/* 推定された姿勢に直す
 * param dest  : 姿勢誤差修正されたデータ
 * param src   : オリジナルデータ
 * param dest2 : 誤差の平均 */
void KinematicPair::RefineOri0(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd> &src, double &dest2) const
{
    MotionMatrixd pos;
    pos.R() = RPY2Mat(param.back());
    dest2 = 0.0;
    std::vector<MotionMatrixd>::const_iterator it;
    for (it = src.begin();it != src.end(); it++){
        dest2 += DiffOri(pos.R(), it->R());
        pos.T() = it->T();
        dest.push_back(pos);
    }
    dest2 /= dest.size();
}

/* 推定された軸向きから姿勢誤差を取り除く
 * (ちなみにLFO+Van/error_correct.cppの中のMoveMultRootと同じことをやっている)
 * param dest  : 姿勢誤差修正されたデータ
 * param src   : オリジナルデータ
 * param dest2 : 誤差の平均 */
void KinematicPair::RefineOri1(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd> &src, double &dest2) const
{
    dest2 = 0.0;
    std::vector<MotionMatrixd>::const_iterator it;
    for (it = src.begin(); it != src.end(); it++){
        const Eigen::Matrix3d R = it->R();
        const Eigen::Vector3d a_b = R * param[0];
        Eigen::Vector3d axis = a_b.cross(param[1]);
        double angle = asin(ProperTriFunc(axis.norm()));
        dest2 += fabs(angle);
		if (param[1].dot(a_b) < 0){ 
			angle = -angle; 
		}
        axis /= axis.norm();
		const Eigen::Matrix3d R1 = Eigen::AngleAxisd(angle, axis).matrix(); 
        MotionMatrixd pos;
		pos.R() = R1 * R;
		pos.T() = it->T();
        dest.push_back(pos);
    }
    dest2 /= dest.size();
}

/* 推定された2軸向きから姿勢誤差を取り除く
 * (ちなみにLFO+Van/error_correct.cppの中のMoveMultRootと同じことをやっている)
 * param dest  : 姿勢誤差修正されたデータ
 * param src   : オリジナルデータ
 * param dest2 : 誤差の平均 */
void KinematicPair::RefineOri2(std::vector<MotionMatrixd > &dest, const std::vector<MotionMatrixd > &src, double &dest2) const
{
	const Eigen::Matrix3d al = KinematicPair::RotXZ(param[0][0], param[0][1]);
	const Eigen::Matrix3d bl = KinematicPair::RotXY(param[1][0], param[1][1]);
    dest2 = 0.0;
    std::vector<MotionMatrixd >::const_iterator it;
    for (it = src.begin(); it != src.end(); it++){
		const Eigen::Matrix3d R = it->R();
        const Eigen::Vector3d a1 = R * al * Eigen::Vector3d::UnitY();
        const Eigen::Vector3d a2 = bl * Eigen::Vector3d::UnitZ();
        Eigen::Vector3d axis = a1.cross(a2);
        double angle = asin(ProperTriFunc(axis.norm())) - M_PI / 2.0;
        dest2 += fabs(angle);
		if (a1.dot(a2) < 0){ 
			angle = -angle; 
		}
        axis /= axis.norm();
		const Eigen::Matrix3d R1 = Eigen::AngleAxisd(angle, axis).matrix();
		MotionMatrixd pos;
		pos.R() = R1 * R;
		pos.T() = it->T();
        dest.push_back(pos);
    }
    dest2 /= dest.size();
}

void KinematicPair::RefineOri3(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd> &src, double &dest2) const
{
    dest2 = 0.0; // no error
    dest = src; // just copy
}

/* ロケーションに関するノイズを取り除く */
void KinematicPair::RefineOri0Loc1(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd> &src, double &dest2) const
{
    Eigen::Vector3d cen = Eigen::Vector3d::Zero();
    std::vector<MotionMatrixd >::const_iterator it;
	for (it = src.begin(); it != src.end(); it++){ 
		cen += it->T(); 
	}
    cen /= static_cast<double>(src.size());
    dest2 = 0.0;
    for (it = src.begin(); it != src.end(); it++){
		MotionMatrixd pos;
		pos.R() = it->R();
        const Eigen::Vector3d temp = it->T() - cen;
        const double dot = param.back().dot(temp);
		pos.T() = cen + dot * param.back();
		dest2 += (pos.T() - it->T()).norm();
        dest.push_back(pos);
    }
    dest2 /= dest.size();
}

/* 推定された軸位置から位置誤差を取り除く
 * param dest  : 位置姿勢誤差修正されたデータ
 * param src   : 姿勢誤差修正されたデータ
 * param dest2 : 誤差の平均 */
void KinematicPair::RefineLoc0(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd> &src, double &dest2) const
{
    dest2 = 0.0;
	dest.resize(src.size());
    for (size_t i(0); i < src.size(); ++i){
		dest[i].R() = src[i].R();
		dest[i].T() = param.back() - src[i].R() * param[param.size() - 2];
		dest2 += (dest[i].T() - src[i].T()).norm();
    }
    dest2 /= dest.size();
}

void KinematicPair::RefineOri1Hel1(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd> &src, const std::vector<double> &src2, double &dest2) const
{
    assert(src.size() == src2.size());
	const Eigen::Vector3d axis = RPY2Mat(param[1]) * Eigen::Vector3d::UnitZ();
    dest2 = 0.0;
    for (size_t i(0); i < src.size(); ++i){
		MotionMatrixd pos;
		pos.R() = src[i].R();
		pos.T() = param[3] + param[4][0] * src2[i] * axis - pos.R() * param[2];
		dest2 += (pos.T() - src[i].T()).norm();
        dest.push_back(pos);
    }
    dest2 /= dest.size();
}

void KinematicPair::RefineLoc3(std::vector<MotionMatrixd > &dest, const std::vector<MotionMatrixd > &src, double &dest2) const
{
    dest2 = 0.0; // エラーなんて存在しないよ
    dest = src; // 単なるコピー
}


#if 0
ScrewKP::ScrewKP() : KinematicPair()
{
}

void ScrewKP::Estimation(std::vector<MotionMatrixd > &dest,
						 const std::vector<MotionMatrixd > &src, EstError &dest2)
{
    std::vector<MotionMatrixd > interTrj;
    /* データの削除 */
    dest.clear();
    param.clear();
    /* 軸向きパラメータ推定 */
    EstOri1Param(interTrj, src, dest2. rotDiv);
    /* 適当な座標系を設定する */
    param[0] = SetCoordinate(param[0]);
    param[1] = SetCoordinate(param[1]);    
    /* */
    EstOri1Hel1Param(dest, interTrj, dest2.transDiv);
}

/* ネジの向き
 * @param src : 1=基準物体 0=対象物体 */
Eigen::Matrix3d ScrewKP::Dic(int src) const
{
    assert(0<=src && src<2);
    Eigen::Matrix3d R = RPY2Mat(param[src]);
	Eigen::Matrix3d z;
	z[0] = 0.0; z[1] = 0.0; z[2] = 1.0;
    return R * z;
}

/* ネジ座標中心
 * @param src : 1=基準物体 0=対象物体 */
Eigen::Matrix3d ScrewKP::Cen(int src) const
{
    assert(0 <= src && src < 2);
    return param[src + 2];
}

/* ネジ座標姿勢(rpy) 
 * @param src : 0=基準物体 1=対象物体 */
Eigen::Matrix3d ScrewKP::Ori(int src) const
{
    assert(0 <= src && src < 2);
    return param[src];
}

/* 並進/回転 */
double ScrewKP::Ratio() const
{
    return param[4][0];
}

void ScrewKP::Print(std::ostream &o) const
{
    Eigen::Matrix3d z;
	z[0] = 0.0; z[1] = 0.0; z[2] = 1.0;
    o << "Target Object: " << param[2] << " : " << param[0] << " : "
      << RPY2Mat(param[0]) * z << std::endl;
    o << "Base Object: " << param[3] << " : " << param[1] << " : "
      << RPY2Mat(param[1]) * z << std::endl;
    o << "Trans/rot: " << param[4][0] << std::endl;
}

/* パラメータを表示する */
std::ostream& operator<<(std::ostream &o, const ScrewKP &src)
{
    src.Print(o);
    return o;
}

/* 固定点中心回転対偶用 */
Ori1DLoc3DKP::Ori1DLoc3DKP() : KinematicPair()
{
}

void Ori1DLoc3DKP::Estimation(std::vector<MotionMatrixd > &dest,
								const std::vector<MotionMatrixd > &src, EstError &dest2)
{
    /* データの削除 */
    dest.clear();
    param.clear();
    /* 軸向きパラメータ推定 */
    EstOri1Param(dest, src, dest2.rotDiv);
	dest2.transDiv = 0.0; // エラーなしということ
}

/* 回転軸向き
 * @param src : 0=基準物体 1=対象物体 */
Eigen::Matrix3d Ori1DLoc3DKP::Dic(int src) const
{
    assert(0 <= src && src < 2);
    return param[src];
}

/* ダミー(使えないし使わないで)
 * @param src : 0=基準物体 1=対象物体 */
Eigen::Matrix3d Ori1DLoc3DKP::Cen(int src) const
{
    assert(0 == 1);
    return Eigen::Matrix3d(0);
}

/* ダミー(使えないし使わないで)
 * @param src : 0=基準物体 1=対象物体 */
Eigen::Matrix3d Ori1DLoc3DKP::Ori(int src) const
{
    assert(0 == 1);
    return Eigen::Matrix3d(0);
}

/* ダミー(使えないし使わないで) */
double Ori1DLoc3DKP::Ratio() const
{
    assert(0 == 1);
    return 0.0;
}

void Ori1DLoc3DKP::Print(std::ostream &o) const
{
    o << "Target Object: " << param[0] << std::endl;
    o << "Base Object: " << param[1] << std::endl;
}

/* パラメータを表示する */
std::ostream& operator<<(std::ostream &o, const Ori1DLoc3DKP &src)
{
    src.Print(o);
    return o;
}

/* 固定姿勢対偶用 */
Ori0DLoc3DKP::Ori0DLoc3DKP() : KinematicPair()
{
}

void Ori0DLoc3DKP::Estimation(std::vector<MotionMatrixd > &dest,
								const std::vector<MotionMatrixd > &src, EstError &dest2)
{
    /* データの削除 */
    dest.clear();
    param.clear();
    /* 軸向きパラメータ推定 */
    EstOri0Param(dest, src, dest2.rotDiv);
	dest2.transDiv = 0.0; // エラーなしということ
}

/* ダミー(使えないし使わないで)
 * @param src : 0=基準物体 1=対象物体 */
Eigen::Matrix3d Ori0DLoc3DKP::Dic(int src) const
{
    assert(0 == 1);
    return Eigen::Matrix3d(0);
}

/* ダミー(使えないし使わないで)
 * @param src : 0=基準物体 1=対象物体 */
Eigen::Matrix3d Ori0DLoc3DKP::Cen(int src) const
{
    assert(0 == 1);
    return Eigen::Matrix3d(0);
}

/* ダミー(使えないし使わないで)
 * @param src : 0=基準物体 1=対象物体 */
Eigen::Matrix3d Ori0DLoc3DKP::Ori(int src) const
{
    assert(0 == 1);
    return Eigen::Matrix3d(0);
}

/* ダミー(使えないし使わないで)
 * @param src : 0=基準物体 1=対象物体 */
double Ori0DLoc3DKP::Ratio() const
{
    assert(0 == 1);
    return 0.0;
}

void Ori0DLoc3DKP::Print(std::ostream &o) const
{
    o << "Target Object: " << param[0] << std::endl;
    o << "Base Object: " << param[1] << std::endl;
}

/* パラメータを表示する */
std::ostream& operator<<(std::ostream &o, const Ori0DLoc3DKP &src)
{
    src.Print(o);
    return o;
}
#endif