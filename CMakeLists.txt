PROJECT(EstimateJointParameters)

CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

set(USE_PCL 0)  # 0 is not used, 1 is used	

if (USE_PCL)
	# PCL
	find_package(PCL 1.6 REQUIRED)
	include_directories(${PCL_INCLUDE_DIRS})
	link_directories(${PCL_LIBRARY_DIRS})
	add_definitions(${PCL_DEFINITIONS})
	add_definitions(/EHsc)
	add_definitions(/D_USE_PCL)
else(USE_PCL)
	#boost
	set(Boost_INCLUDE_DIR $ENV{Boost_INCLUDE_DIR})
	set(Boost_LIBRARY_DIR $ENV{Boost_LIBRARY_DIR})
	find_package(Boost REQUIRED)
	include_directories(${Boost_INCLUDE_DIR})
	link_directories(${Boost_LIBRARY_DIR})
	# Eigen
	set(CMAKE_MODULE_PATH $ENV{Eigen3_DIR}/cmake)
	set(EIGEN3_INCLUDE_DIR $ENV{Eigen3_DIR})
	FIND_PACKAGE(Eigen3 REQUIRED)
	include_directories(${EIGEN3_INCLUDE_DIR})
endif(USE_PCL)

if (COMMAND cmake_policy)#PCL
	cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

#gsl
set(GSL_DIR C:/local/gsl-2.2/msvc2017_64)
include_directories(${GSL_DIR}/include)
link_directories(${GSL_DIR}/lib/gsl)
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

ADD_EXECUTABLE(
	EstimateRevoluteJointUsingMarker
	sample4.cpp
	KinematicPair.cpp
	KinematicPair.h
	RevolutePair.cpp
	RevolutePair.h
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
TARGET_LINK_LIBRARIES(EstimateRevoluteJointUsingMarker ${GSL_LIBS} ${OpenCV_LIBS})

