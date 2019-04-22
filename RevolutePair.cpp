// RevolutePair.cpp
#include "RevolutePair.h"
#include "eigen_matrix_utility.h"
#include <iostream>
#include <gsl/gsl_multimin.h>

void RevolutePair::Estimation(std::vector<MotionMatrixd > &dest, const std::vector<MotionMatrixd > &src, EstError &dest2)
{
    /* clear */
    dest.clear();
    param.clear();
    /* estimate rotation axis */
    std::vector<MotionMatrixd> interTrj;
    EstOri1Param(interTrj, src, dest2.rotDiv);
    /* estimate center of rotation */
    EstLoc0Param(dest, interTrj, dest2.transDiv);
}

void RevolutePair::EstimationConsideringFreeDOF(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd> & src, const std::vector<MotionMatrixd> &orig, double &trans_err, size_t MAX_LOOP)
{
	dest.resize(src.size());
	for (size_t i(0); i < src.size(); ++i){
		dest[i].R() = src[i].R();
	}
	for (size_t i(0); i < MAX_LOOP; ++i){
		std::vector<MotionMatrixd> tmp;
		/* deep estimation */
		for (size_t j(0); j < dest.size(); ++j){		
			dest[j].T() = orig[j].T();
		}
		EstOriFreeParam2(tmp, dest);
		param.pop_back();
		param.pop_back();
		EstLoc0Param(dest, tmp, trans_err);
	}
}

struct NecessaryInfoForEval
{
	MotionMatrixd pos;
	Eigen::Vector3d lb;
	Eigen::Vector3d ca;
	Eigen::Vector3d cb;
	double w;
};

static double EstOriFreeParamEval(const ::gsl_vector *xvec_ptr, void *params)
{
	const NecessaryInfoForEval *param_ptr = (NecessaryInfoForEval*) params;

	const double theta = ::gsl_vector_get(xvec_ptr, 0);
	const Eigen::Matrix3d dr = Eigen::AngleAxisd(theta, param_ptr->lb).matrix(); 

	Eigen::Vector3d dt = param_ptr->cb - dr * param_ptr->pos.R() * param_ptr->ca - param_ptr->pos.T();
	return 0.5 * (dt.squaredNorm() + param_ptr->w * theta * theta); 
}

static void DEstOriFreeParamEval(const ::gsl_vector *xvec_ptr, void *params, ::gsl_vector *df_ptr)
{
	const NecessaryInfoForEval *param_ptr = (NecessaryInfoForEval*) params;

	const double theta = ::gsl_vector_get(xvec_ptr, 0);
	const Eigen::Matrix3d dr = Eigen::AngleAxisd(theta, param_ptr->lb).matrix(); 
	const Eigen::Matrix3d ddr = dAngleAxis2Mat(param_ptr->lb, theta);

	Eigen::Vector3d dt = param_ptr->cb - dr * param_ptr->pos.R() * param_ptr->ca - param_ptr->pos.T();
	Eigen::Vector3d ddt = -ddr * param_ptr->pos.R() * param_ptr->ca;
	double dest = dt.dot(ddt) + param_ptr->w * theta;
	::gsl_vector_set(df_ptr, 0, dest);

	return;
}

static void FDEstOriFreeParamEval(const ::gsl_vector *xvec_ptr, void *params_ptr, double *f_ptr, ::gsl_vector *df_ptr)
{
	*f_ptr = EstOriFreeParamEval(xvec_ptr, params_ptr);
	DEstOriFreeParamEval(xvec_ptr, params_ptr, df_ptr);
}

void RevolutePair::EstOriFreeParam2(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd> &src) const
{
	NecessaryInfoForEval param;
	param.ca = CenterOfRotation(KinematicPair::TARGET);
	param.cb = CenterOfRotation(KinematicPair::BASE);
	param.lb = AxisDirection(KinematicPair::BASE);
	param.w = 1.0;	

	::gsl_multimin_function_fdf my_func;
	my_func.n = 1;
	my_func.f = &EstOriFreeParamEval;
	my_func.df = &DEstOriFreeParamEval;
	my_func.fdf = &FDEstOriFreeParamEval;
	my_func.params = (void*) &param;
	dest.resize(src.size());
	for (size_t i(0); i < dest.size(); ++i){
		param.pos = src[i];
		/* set initial guess */
		::gsl_vector *vec_p = ::gsl_vector_alloc(1);
		::gsl_vector_set(vec_p, 0, 0.0);
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
		Eigen::Matrix3d rot = Eigen::AngleAxisd(::gsl_vector_get(minimizer_ptr->x, 0), param.lb).matrix(); 
		dest[i].R() = rot * src[i].R();
		dest[i].T() = src[i].T();
		// free 
		::gsl_multimin_fdfminimizer_free(minimizer_ptr);
		::gsl_vector_free(vec_p);
	}
}

void RevolutePair::EstOriFreeParam(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd> &src) const 
{
	const double _NEARLYZERO = 1.0e-6;
	const double degree_thresh = M_PI / 36.0; // 5 degrees
	Eigen::Matrix3d trans;
	const Eigen::Vector3d lb = AxisDirection(KinematicPair::BASE); 
	ZAxis2Coordinate(trans, lb);
	trans.transposeInPlace();
	const Eigen::Vector3d ca = CenterOfRotation(KinematicPair::TARGET);
	const Eigen::Vector3d cb = CenterOfRotation(KinematicPair::BASE);
	dest.resize(src.size());
	for (size_t i(0); i < dest.size(); ++i){
		Eigen::Vector3d cat = trans * src[i].R() * ca;
		Eigen::Vector3d cbt = trans * (cb - src[i].T());
		cat[2] = cbt[2] = 0.0;
		if (cat.norm() < _NEARLYZERO || cbt.norm() < _NEARLYZERO){
			dest[i] = src[i];
			continue;
		}
		double diff = atan2(cat[1], cat[0]) - atan2(cbt[1], cbt[0]);
		if (fabs(diff) >= degree_thresh){
			dest[i] = src[i];
			continue;
		}
		Eigen::Matrix3d rot = Eigen::AngleAxisd(diff, lb).matrix(); 
		dest[i].R() = rot * src[i].R();
		dest[i].T() = src[i].T();
#if 0
		// debug
		const Eigen::Vector3d cat2 = trans * dest[i].R() * ca;
		const Eigen::Vector3d cbt2 = trans * (cb - dest[i].T());
		double diff2 = atan2(cat2[1], cat2[0]) - atan2(cbt2[1], cbt2[0]);
		std::cout << diff << " " << diff2 << " " << (cat - cbt).norm() << " " << (cat2 - cbt2).norm() << std::endl;
		std::cout << (src[i].R() * la).transpose() << " " << lb.transpose() << std::endl;
		std::cout << (dest[i].R() * la).transpose() << " " << lb.transpose() << std::endl;
		getchar();
#endif
	}
}

/* ‰ñ“]Ž²‚ÌŒü‚«
 * @param src : 1=Šî€•¨‘Ì 0=‘ÎÛ•¨‘Ì */
const Eigen::Vector3d& RevolutePair::AxisDirection(int src) const
{
    assert(0 <= src && src < 2);
    return param[src];
}

/* ‰ñ“]’†S
 * @param src : 1=Šî€•¨‘Ì 0=‘ÎÛ•¨‘Ì */
const Eigen::Vector3d& RevolutePair::CenterOfRotation(int src) const
{
    assert(0 <= src && src < 2);
    return param[src + 2];
}

/* print parameter */
void RevolutePair::Print(std::ostream &o) const
{
    o << "Target Object: " << param[0] << " : " << param[2] << std::endl;
    o << "Base Object: " << param[1] << " : " << param[3] << std::endl;
}

std::ostream& operator<<(std::ostream &o,const RevolutePair &src)
{
    src.Print(o);
    return o;
}
