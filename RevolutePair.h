// RevolutePair.h
#pragma once

#include "KinematicPair.h"

/* revolute pair */
class RevolutePair : public KinematicPair
{
public:
	RevolutePair() : KinematicPair(){}
    void Estimation(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd> &src, EstError &dest2);
	void EstimationConsideringFreeDOF(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd> & src, const std::vector<MotionMatrixd> &orig, double &trans_err, size_t MAX_LOOP = 10);
	void Print(std::ostream &dest) const;
    const Eigen::Vector3d& AxisDirection(int src) const; /* ‰ñ“]²‚ÌŒü‚« */
    const Eigen::Vector3d& CenterOfRotation(int src) const; /* ‰ñ“]’†S */
	KinematicPair::TypeOfJoint GetType() const{ return _REVOLUTE; }  
    friend std::ostream& operator<<(std::ostream &dest, const RevolutePair &src);
private:
	void EstOriFreeParam(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd> &src) const;
	void EstOriFreeParam2(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd> &src) const;
};


