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
 - PolynomialSolver (get from git server)
 - OpenCV (used in PolynomialSolver)
 - UtilForEigen (get from git server)
 - GNU Scientific Library (downloaded from the web)

For compiling samples
 - Point Cloud Library (downloaded from the web, to use visualization)
 - Boost

# Compile
1. Use cmake to make project file
2. just compile it

# Test samples
Just run the execution with no command parameters
