// SphericalPair.h

#pragma once

#include "KinematicPair.h"

/* Spherical pair */
class SphericalPair : public KinematicPair
{
public:
	SphericalPair() : KinematicPair(){}
    void Estimation(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd> &src, EstError &dest2);
	void EstimationConsideringFreeDOF(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd> & src, const std::vector<MotionMatrixd> &orig, double &trans_err, size_t MAX_LOOP = 10);
    void Print(std::ostream &dest) const;
    Eigen::Vector3d CenterOfRotation(int src) const; /* center of roration */
	KinematicPair::TypeOfJoint GetType() const{ return _SPHERICAL; }  	
	friend std::ostream& operator<<(std::ostream &dest,const SphericalPair &src);
private:
	void EstOriFreeParam(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd> &src, const std::vector<MotionMatrixd> &orig) const;
};
