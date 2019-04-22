// SphericalPair.cpp

#include "SphericalPair.h"
#include "eigen_matrix_utility.h"
#include <gsl/gsl_multimin.h>
#include <iostream>

void SphericalPair::Estimation(std::vector<MotionMatrixd > &dest, const std::vector<MotionMatrixd > &src, EstError &dest2)
{
    std::vector<MotionMatrixd> interTrj;
    /* clear */
    dest.clear();
    param.clear();
    /* estimate orientation part */
    EstOri3Param(interTrj, src, dest2.rotDiv);
    /* estimate center of rotation */
    EstLoc0Param(dest, interTrj, dest2.transDiv);
}

void SphericalPair::EstimationConsideringFreeDOF(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd> & src, const std::vector<MotionMatrixd> &orig, double &trans_err, size_t MAX_LOOP)
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

struct NecessaryInfoForEval
{
	MotionMatrixd pos;
	Eigen::Matrix3d orig_R;
	Eigen::Vector3d ca;
	Eigen::Vector3d cb;
	double w;
};

static double EstOriFreeParamEval(const ::gsl_vector *xvec_ptr, void *params)
{
	const NecessaryInfoForEval *param_ptr = (NecessaryInfoForEval*) params;
	
	const Eigen::Vector3d rpy(::gsl_vector_get(xvec_ptr, 0), ::gsl_vector_get(xvec_ptr, 1), ::gsl_vector_get(xvec_ptr, 2));
	const Eigen::Matrix3d dR = ::RPY2Mat(rpy);
	const Eigen::Matrix3d R = dR * param_ptr->orig_R;
	const Eigen::Vector3d dt = param_ptr->cb - R * param_ptr->ca - param_ptr->pos.T();
	const double ang = dR.trace(); // 1 + 2 cos theta
	return 0.5 * dt.squaredNorm() - param_ptr->w * ang; 
}

static void DEstOriFreeParamEval(const ::gsl_vector *xvec_ptr, void *params, ::gsl_vector *df_ptr)
{
	const NecessaryInfoForEval *param_ptr = (NecessaryInfoForEval*) params;

	const Eigen::Vector3d rpy(::gsl_vector_get(xvec_ptr, 0), ::gsl_vector_get(xvec_ptr, 1), ::gsl_vector_get(xvec_ptr, 2));
	const Eigen::Matrix3d dR  = ::RPY2Mat(rpy);
	const Eigen::Matrix3d dR1 = ::RPYdR2Mat(rpy);
	const Eigen::Matrix3d dR2 = ::RPYdP2Mat(rpy);
	const Eigen::Matrix3d dR3 = ::RPYdY2Mat(rpy);
	const Eigen::Matrix3d R  = dR * param_ptr->orig_R;
	const Eigen::Matrix3d R1 = dR1 * param_ptr->orig_R;
	const Eigen::Matrix3d R2 = dR2 * param_ptr->orig_R;
	const Eigen::Matrix3d R3 = dR3 * param_ptr->orig_R;

	const Eigen::Vector3d dt  = param_ptr->cb - R * param_ptr->ca - param_ptr->pos.T();
	const Eigen::Vector3d dt1 = -R1 * param_ptr->ca;
	const Eigen::Vector3d dt2 = -R2 * param_ptr->ca;
	const Eigen::Vector3d dt3 = -R3 * param_ptr->ca;

	
	const double ang1 = dR1.trace();
	const double ang2 = dR2.trace();
	const double ang3 = dR3.trace();

	const double dest1 = dt.dot(dt1) - param_ptr->w * ang1;
	const double dest2 = dt.dot(dt2) - param_ptr->w * ang2;
	const double dest3 = dt.dot(dt3) - param_ptr->w * ang3;
	::gsl_vector_set(df_ptr, 0, dest1);
	::gsl_vector_set(df_ptr, 1, dest2);
	::gsl_vector_set(df_ptr, 2, dest3);

	return;
}

static void FDEstOriFreeParamEval(const ::gsl_vector *xvec_ptr, void *params_ptr, double *f_ptr, ::gsl_vector *df_ptr)
{
	*f_ptr = EstOriFreeParamEval(xvec_ptr, params_ptr);
	DEstOriFreeParamEval(xvec_ptr, params_ptr, df_ptr);
}

void SphericalPair::EstOriFreeParam(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd> &src, const std::vector<MotionMatrixd> &orig) const
{
	// reset translation
	dest.resize(src.size());
	for (size_t i(0); i < dest.size(); ++i){		
		dest[i].R() = src[i].R();
		dest[i].T() = orig[i].T();
	}

	NecessaryInfoForEval param;
	param.ca = CenterOfRotation(KinematicPair::TARGET);
	param.cb = CenterOfRotation(KinematicPair::BASE);
	param.w = 1;	

	::gsl_multimin_function_fdf my_func;
	my_func.n = 3;
	my_func.f = &EstOriFreeParamEval;
	my_func.df = &DEstOriFreeParamEval;
	my_func.fdf = &FDEstOriFreeParamEval;
	my_func.params = (void*) &param;
	dest.resize(src.size());
	for (size_t i(0); i < dest.size(); ++i){
		param.pos = dest[i];
		param.orig_R = orig[i].R();
		/* set initial guess */
		::gsl_vector *vec_p = ::gsl_vector_alloc(3);
		::gsl_vector_set(vec_p, 0, 0.0);
		::gsl_vector_set(vec_p, 1, 0.0);
		::gsl_vector_set(vec_p, 2, 0.0);
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
		const Eigen::Vector3d rpy(::gsl_vector_get(minimizer_ptr->x, 0), ::gsl_vector_get(minimizer_ptr->x, 1), ::gsl_vector_get(minimizer_ptr->x, 2));
		const Eigen::Matrix3d dR = ::RPY2Mat(rpy);
		dest[i].R() = dR * orig[i].R();
		// free 
		::gsl_multimin_fdfminimizer_free(minimizer_ptr);
		::gsl_vector_free(vec_p);
	}
}

/* 回転中心
 * @param src : 0=基準物体 1=対象物体 */
Eigen::Vector3d SphericalPair::CenterOfRotation(int src) const
{
    assert(0 <= src && src < 2);
    return param[src];
}

void SphericalPair::Print(std::ostream &o) const
{
    o << "Target Object: " << param[0] << std::endl;
    o << "Base Object: " << param[1] << std::endl;
}

/* パラメータを表示する */
std::ostream& operator<<(std::ostream &o, const SphericalPair &src)
{
    src.Print(o);
    return o;
}
