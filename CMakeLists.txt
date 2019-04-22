PROJECT(EstimateJointParameters)

CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

#PCL
find_package(PCL 1.6 REQUIRED)
include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})
add_definitions(/EHsc)

if (COMMAND cmake_policy)
	cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

#gsl
set(GSL_DIR C:/gsl-1.16)
include_directories(${GSL_DIR})
link_directories(${GSL_DIR}/Win32/Release)
set(GSL_LIBS gsl.lib cblas.lib)

#my motion analysis program
set(MY_UTILS ../UtilsForEigen)
set(MY_POLY ../PolynomialSolver)
include_directories(${MY_UTILS} ${MY_POLY})

#opencv
find_package(OpenCV REQUIRED)

# for getting learning data
ADD_EXECUTABLE(
	EstimateRevoluteJoint
	sample1.cpp
	KinematicPair.cpp
	KinematicPair.h
	RevolutePair.cpp
	RevolutePair.h
	${MY_POLY}/Monomial.cpp
	${MY_POLY}/MonomialParser.cpp
	${MY_UTILS}/MotionMatrix.h
	${MY_UTILS}/eigen_matrix_utility.h
)

ADD_EXECUTABLE(
	EstimateDoubleRevoluteJoint
	sample2.cpp
	KinematicPair.cpp
	KinematicPair.h
	DoubleRevolutePair.cpp
	DoubleRevolutePair.h
	${MY_POLY}/Monomial.cpp
	${MY_POLY}/MonomialParser.cpp
	${MY_UTILS}/MotionMatrix.h
	${MY_UTILS}/eigen_matrix_utility.h
)

ADD_EXECUTABLE(
	EstimateSphericalJoint
	sample3.cpp
	KinematicPair.cpp
	KinematicPair.h
	SphericalPair.cpp
	SphericalPair.h
	${MY_POLY}/Monomial.cpp
	${MY_POLY}/MonomialParser.cpp
	${MY_UTILS}/MotionMatrix.h
	${MY_UTILS}/eigen_matrix_utility.h
)

# compile option for msvc
IF(MSVC)
	SET(CMAKE_EXE_LINKER_FLAGS "/NODEFAULTLIB:LIBCMT.lib")
ENDIF()

TARGET_LINK_LIBRARIES(EstimateRevoluteJoint ${GSL_LIBS} ${OpenCV_LIBS} ${PCL_LIBRARIES})
TARGET_LINK_LIBRARIES(EstimateDoubleRevoluteJoint ${GSL_LIBS} ${OpenCV_LIBS} ${PCL_LIBRARIES})
TARGET_LINK_LIBRARIES(EstimateSphericalJoint ${GSL_LIBS} ${OpenCV_LIBS} ${PCL_LIBRARIES})

