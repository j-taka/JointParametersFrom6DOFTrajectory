// PrismaticPair.h
#pragma once

#include "KinematicPair.h"

class PrismaticPair : public KinematicPair
{
public:
	PrismaticPair() : KinematicPair() {}
    void Estimation(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd> &src, EstError &dest2);
    void Print(std::ostream &dest) const;
    Eigen::Vector3d TranslationDirection(int src) const; /* ï¿êiï˚å¸ */
	KinematicPair::TypeOfJoint GetType() const{ return _PRISMATIC; }	
	friend std::ostream& operator<<(std::ostream &dest, const PrismaticPair &src);
};
