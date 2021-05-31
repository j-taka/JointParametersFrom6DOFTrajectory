This program is for estimating joint parameters from 6DOF relative trajectory between two links.

See the details from 
1. Yoshihiro Sato, Jun Takamatsu, Hiroshi Kimura, Katsushi Ikeuchi, “Recognition of a Mechanical Linkage Based on Occlusion-Robust Object Tracking,” IEEE International Conference on Multisensor Fusion and Integration for Intelligent Systems (MFI2003) , 2003.
2. Jun Takamatsu, “Abstraction of Manipulation Tasks to Automatically Generate Robot Motion from Observation,” Ph.D. Thesis , the University of Tokyo, Feb. 2004

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
 - [UtilForEigen](https://github.com/j-taka/UtilsForEigen)
 - [GNU Scientific Library](https://www.bruot.org/hp/libraries/)

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
> EstimateRevoluteJointUsingMarker.exe -d=10 -l=0.141 -c="calib.xml" -v="door.mp4" -p="est_param.txt"

![output_of_sample4](https://user-images.githubusercontent.com/11922392/60444303-e6d02800-9c57-11e9-84a0-abaa4d8602bc.gif)

After the estimation using sample4
> DrawRevoluteJointFromMarker.exe -d=10 -l=0.141 -c="calib.xml" -v="door2.mp4" -p="est_param.txt"

![output_of_sample5](https://user-images.githubusercontent.com/11922392/60444304-e6d02800-9c57-11e9-86ca-5ec04adf5eda.gif)

