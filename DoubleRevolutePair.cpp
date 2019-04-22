// DoubleRevolutePair.cpp
#include "DoubleRevolutePair.h"
#include "eigen_matrix_utility.h"
#include <cassert>
#include <iostream>
#include <gsl/gsl_multimin.h>
#include <cstdlib>

void DoubleRevolutePair::Estimation(std::vector<MotionMatrixd > &dest, const std::vector<MotionMatrixd > &src, EstError &dest2)
{
    std::vector<MotionMatrixd> interTrj;
    /* clear */
    dest.clear();
    param.clear();
    /* estimate orientation */
    EstOri2Param(interTrj, src, dest2.rotDiv);
    /* estimate location */
    EstLoc0Param(dest, interTrj, dest2.transDiv);
}

void DoubleRevolutePair::EstimationConsideringFreeDOF(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd> & src, const std::vector<MotionMatrixd> &orig, double &trans_err, size_t MAX_LOOP)
{
	dest.resize(src.size());
	for (size_t i(0); i < src.size(); ++i){
		dest[i].R() = src[i].R();
	}
	for (size_t i(0); i < MAX_LOOP; ++i){
		std::vector<MotionMatrixd> tmp;
		/* deep estimation */
		EstOriFreeParam(tmp, dest, orig);
		param.pop_back();
		param.pop_back();
		EstLoc0Param(dest, tmp, trans_err);
	}
}

void DoubleRevolutePair::RotZY2Angles(Eigen::Vector2d &dest, const Eigen::Matrix3d &src) const
{
	assert(fabs(src(2, 1)) < _NEARLY_ZERO); 
	dest[0] = atan2(-src(0, 1), src(1, 1));
	dest[1] = atan2(-src(2, 0), src(2, 2));
}

void DoubleRevolutePair::GetCurrentRotation(Eigen::Vector2d &dest, const MotionMatrixd &src) const
{
	Eigen::Matrix3d al = this->Orientation(KinematicPair::TARGET);
	Eigen::Matrix3d lb = this->Orientation(KinematicPair::BASE).transpose();
	Eigen::Matrix3d tmp = lb * src.R() * al;
	RotZY2Angles(dest, tmp);
}

struct NecessaryInfoForEval
{
	MotionMatrixd pos;
	Eigen::Matrix3d orig_R;
	Eigen::Matrix3d al;
	Eigen::Matrix3d bl;
	Eigen::Vector3d ca;
	Eigen::Vector3d cb;
	double w;
};

static double EstOriFreeParamEval(const ::gsl_vector *xvec_ptr, void *params)
{
	const NecessaryInfoForEval *param_ptr = (NecessaryInfoForEval*) params;
	
	const double tz = ::gsl_vector_get(xvec_ptr, 0);
	const double ty = ::gsl_vector_get(xvec_ptr, 1);
	const Eigen::Matrix3d Rzy = Eigen::AngleAxisd(tz, Eigen::Vector3d::UnitZ()).matrix() * Eigen::AngleAxisd(ty, Eigen::Vector3d::UnitY()).matrix();
	const Eigen::Matrix3d R = param_ptr->bl * Rzy * param_ptr->al.transpose();
	const Eigen::Vector3d dt = param_ptr->cb - R * param_ptr->ca - param_ptr->pos.T();
	const double ang = (R * param_ptr->orig_R.transpose()).trace(); // 1 + 2 cos theta
	return 0.5 * dt.squaredNorm() - param_ptr->w * ang; 
}

static void DEstOriFreeParamEval(const ::gsl_vector *xvec_ptr, void *params, ::gsl_vector *df_ptr)
{
	const NecessaryInfoForEval *param_ptr = (NecessaryInfoForEval*) params;

	const double tz = ::gsl_vector_get(xvec_ptr, 0);
	const double ty = ::gsl_vector_get(xvec_ptr, 1);
	const Eigen::Matrix3d Rzy = Eigen::AngleAxisd(tz, Eigen::Vector3d::UnitZ()).matrix() * Eigen::AngleAxisd(ty, Eigen::Vector3d::UnitY()).matrix();
	const Eigen::Matrix3d Rdzy = ::dRotZ2Mat(tz) * Eigen::AngleAxisd(ty, Eigen::Vector3d::UnitY()).matrix();
	const Eigen::Matrix3d Rzdy = Eigen::AngleAxisd(tz, Eigen::Vector3d::UnitZ()).matrix() * ::dRotY2Mat(ty);
	const Eigen::Matrix3d R = param_ptr->bl * Rzy * param_ptr->al.transpose();
	const Eigen::Matrix3d R1 = param_ptr->bl * Rdzy * param_ptr->al.transpose();
	const Eigen::Matrix3d R2 = param_ptr->bl * Rzdy * param_ptr->al.transpose();

	const Eigen::Vector3d dt = param_ptr->cb - R * param_ptr->ca - param_ptr->pos.T();
	const Eigen::Vector3d dt1 = - R1 * param_ptr->ca;
	const Eigen::Vector3d dt2 = - R2 * param_ptr->ca;
	
	const double ang1 = (R1 * param_ptr->orig_R.transpose()).trace();
	const double ang2 = (R2 * param_ptr->orig_R.transpose()).trace();

	const double dest1 = dt.dot(dt1) - param_ptr->w * ang1;
	const double dest2 = dt.dot(dt2) - param_ptr->w * ang2;
	::gsl_vector_set(df_ptr, 0, dest1);
	::gsl_vector_set(df_ptr, 1, dest2);

	return;
}

static void FDEstOriFreeParamEval(const ::gsl_vector *xvec_ptr, void *params_ptr, double *f_ptr, ::gsl_vector *df_ptr)
{
	*f_ptr = EstOriFreeParamEval(xvec_ptr, params_ptr);
	DEstOriFreeParamEval(xvec_ptr, params_ptr, df_ptr);
}

void DoubleRevolutePair::EstOriFreeParam(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd> &src, const std::vector<MotionMatrixd> &orig) const
{
	// reset translation
	dest.resize(src.size());
	for (size_t i(0); i < dest.size(); ++i){		
		dest[i].R() = src[i].R();
		dest[i].T() = orig[i].T();
	}

	NecessaryInfoForEval param;
	param.al = Orientation(KinematicPair::TARGET);
	param.bl = Orientation(KinematicPair::BASE);
	param.ca = CenterOfRotation(KinematicPair::TARGET);
	param.cb = CenterOfRotation(KinematicPair::BASE);
	param.w = 1;	

	::gsl_multimin_function_fdf my_func;
	my_func.n = 2;
	my_func.f = &EstOriFreeParamEval;
	my_func.df = &DEstOriFreeParamEval;
	my_func.fdf = &FDEstOriFreeParamEval;
	my_func.params = (void*) &param;
	dest.resize(src.size());
	for (size_t i(0); i < dest.size(); ++i){
		param.pos = dest[i];
		param.orig_R = orig[i].R();
		/* set initial guess */
		Eigen::Vector2d ang;
		GetCurrentRotation(ang, dest[i]);
		::gsl_vector *vec_p = ::gsl_vector_alloc(2);
		::gsl_vector_set(vec_p, 0, ang[0]);
		::gsl_vector_set(vec_p, 1, ang[1]);
		const ::gsl_multimin_fdfminimizer_type *type_ptr = ::gsl_multimin_fdfminimizer_conjugate_fr;
		::gsl_multimin_fdfminimizer *minimizer_ptr = ::gsl_multimin_fdfminimizer_alloc(type_ptr, my_func.n);

		// set the tolerance and starting step size
		double step_size = 1.0e-6;
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
		const double tz = ::gsl_vector_get(minimizer_ptr->x, 0);
		const double ty = ::gsl_vector_get(minimizer_ptr->x, 1);
		const Eigen::Matrix3d Rzy = Eigen::AngleAxisd(tz, Eigen::Vector3d::UnitZ()).matrix() * Eigen::AngleAxisd(ty, Eigen::Vector3d::UnitY()).matrix();
		dest[i].R() = param.bl * Rzy * param.al.transpose();
		// free 
		::gsl_multimin_fdfminimizer_free(minimizer_ptr);
		::gsl_vector_free(vec_p);
	}
}

Eigen::Vector3d DoubleRevolutePair::AxisDirection(int src) const
{
	assert(0 <= src && src < 2);
	if (src == 0){
		return Orientation(0) * Eigen::Vector3d::UnitY();
	}
	else{
		return Orientation(1) * Eigen::Vector3d::UnitZ();
	}
}

/* 回転中心
 * @param src : 1=基準物体 0=対象物体 */
const Eigen::Vector3d& DoubleRevolutePair::CenterOfRotation(int src) const
{
    assert(0 <= src && src < 2);
    return param[src + 2];
}

Eigen::Matrix3d DoubleRevolutePair::Orientation(int src) const
{
	assert(0 <= src && src < 2);
	if (src == 0){
		return KinematicPair::RotXZ(param[0][0], param[0][1]);
	}
	else{
		return KinematicPair::RotXY(param[1][0], param[1][1]);;
	}
}

void DoubleRevolutePair::Print(std::ostream &o) const
{
    Eigen::Vector3d y, z;
	y[0] = 0.0; y[1] = 1.0; y[2] = 0.0;
	z[0] = 0.0; z[1] = 0.0; z[2] = 1.0;
    o << "Target Object: y-Axis: " << RPY2Mat(param[0]) * y << " : " << param[2] << std::endl;
    o << "Base Object: z-Axis: " << RPY2Mat(param[1]) * z << " : " << param[3] << std::endl;
}

/* パラメータを表示する */
std::ostream& operator<<(std::ostream &o, const DoubleRevolutePair &src)
{
    src.Print(o);
    return o;
}

