// DoubleRevolutePair.h

#pragma once

#include "KinematicPair.h"

/* pair with two orthogonal revolute joints */
class DoubleRevolutePair : public KinematicPair
{
public:
	DoubleRevolutePair() : KinematicPair(){}
    void Estimation(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd> &src, EstError &dest2);
	void EstimationConsideringFreeDOF(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd> & src, const std::vector<MotionMatrixd> &orig, double &trans_err, size_t MAX_LOOP = 10);
    void Print(std::ostream &dest) const;
	Eigen::Vector3d AxisDirection(int src) const;
	const Eigen::Vector3d& CenterOfRotation(int src) const;
	Eigen::Matrix3d Orientation(int src) const;
	KinematicPair::TypeOfJoint GetType() const { return _DOUBLE_REVOLUTE; }
	friend std::ostream& operator<<(std::ostream &dest, const DoubleRevolutePair &src);
private:
	void EstOriFreeParam(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd> &src, const std::vector<MotionMatrixd> &orig) const;
	void GetCurrentRotation(Eigen::Vector2d &dest, const MotionMatrixd &src) const;
	void RotZY2Angles(Eigen::Vector2d &dest, const Eigen::Matrix3d &src) const;
};
