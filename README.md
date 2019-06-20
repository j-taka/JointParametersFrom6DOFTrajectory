This program is for estimating joint parameters from 6DOF relative trajectory between two links.

# Introduction

The joint to be hundled are
 - Revolute joint
 - Double revolute joint (two orthogonal joints are connected)
 - Spherical joint

In near future, I will open the function for prismatic joint.
I hope to implement the function for screw joint in future...

# Requirement
 - CMake
 - [PolynomialSolver](https://github.com/j-taka/PolynomialSolver)
 - OpenCV (used in PolynomialSolver)
 - [UtilForEigen] (https://github.com/j-taka/UtilsForEigen)
 - [GNU Scientific Library] (https://www.bruot.org/hp/libraries/)

For compiling samples
 - Boost
 - ArUco (openCV contrib)
 - Point Cloud Library (optional) (downloaded from the web, to use visualization)

# Compile
1. Use cmake to make project file
2. just compile it

# Test samples
Just run the execution with no command parameters from sample1 to sample3
> EstimateRevoluteJoint.exe

> EstimateDoubleRevoluteJoint.exe

> EstimateSphericalJoint.exe

For executing sample4
> EstimateRevoluteJointUsingMarker.exe -d=10 -l=0.141 -c="calib.xml" -v="door.mp4"
