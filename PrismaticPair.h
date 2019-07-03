// PrismaticPair.h
#pragma once

#include "KinematicPair.h"

/* prismatic pair */
class PrismaticPair : public KinematicPair
{
public:
	PrismaticPair() : KinematicPair() {}
	int Save(const std::string &filename) const;
	int Load(const std::string &filename);
	void Estimation(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd> &src, EstError &dest2);
	void Print(std::ostream &dest) const;
    Eigen::Vector3d TranslationDirection(int src) const; /* translation direction */
	KinematicPair::TypeOfJoint GetType() const{ return _PRISMATIC; }	
	friend std::ostream& operator<<(std::ostream &dest, const PrismaticPair &src);
};
