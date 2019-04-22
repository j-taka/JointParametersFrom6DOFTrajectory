// FixedPair.h
#pragma once

#include "KinematicPair.h"

/*! \brief 3D-3D calbibration */
class FixedPair : public KinematicPair
{
public:
	FixedPair() : KinematicPair(){}
    void Estimation(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd> &src, EstError &dest2);
    void Print(std::ostream &dest) const;
    Eigen::Vector3d Translation(int src) const; /* Translation */
    Eigen::Vector3d Rotation(int src) const; /* Rotation */
	KinematicPair::TypeOfJoint GetType() const{ return _FIXED; }  
	friend std::ostream& operator<<(std::ostream &dest, const FixedPair &src);
};
