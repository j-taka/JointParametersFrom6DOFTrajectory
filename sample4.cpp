// sample4.cpp

#include <iostream>
#include "RevolutePair.h"
#include <opencv2/opencv.hpp>
#include <opencv2/aruco.hpp>
#include <opencv2/core/eigen.hpp>

using namespace std;
using namespace cv;

namespace {
	const char* about = "Estimate a revolute joint parameter using a marker";
	const char* keys =
		"{d        |       | dictionary: DICT_4X4_50=0, DICT_4X4_100=1, DICT_4X4_250=2,"
		"DICT_4X4_1000=3, DICT_5X5_50=4, DICT_5X5_100=5, DICT_5X5_250=6, DICT_5X5_1000=7, "
		"DICT_6X6_50=8, DICT_6X6_100=9, DICT_6X6_250=10, DICT_6X6_1000=11, DICT_7X7_50=12,"
		"DICT_7X7_100=13, DICT_7X7_250=14, DICT_7X7_1000=15, DICT_ARUCO_ORIGINAL = 16,"
		"DICT_APRILTAG_16h5=17, DICT_APRILTAG_25h9=18, DICT_APRILTAG_36h10=19, DICT_APRILTAG_36h11=20}"
		"{v        |       | Input from video file, if ommited, input comes from camera }"
		"{c        |       | Camera intrinsic parameters. Needed for camera pose }"
		"{l        | 0.1   | Marker side lenght (in meters). Needed for correct scale in camera pose }"
		"{dp       |       | File of marker detector parameters }"
		"{r        |       | show rejected candidates too }"
		"{refine   |       | Corner refinement: CORNER_REFINE_NONE=0, CORNER_REFINE_SUBPIX=1,"
		"CORNER_REFINE_CONTOUR=2, CORNER_REFINE_APRILTAG=3}";
}

/**
 */
static bool readCameraParameters(const std::string &filename, cv::Mat &camMatrix, cv::Mat &distCoeffs) 
{
	cv::FileStorage fs(filename, cv::FileStorage::READ);
	if (!fs.isOpened()) {
		return false;
	}
	fs["camera_matrix"] >> camMatrix;
	fs["distortion_coefficients"] >> distCoeffs;
	return true;
}

/**
 */
static bool readDetectorParameters(const std::string &filename, cv::Ptr<cv::aruco::DetectorParameters> &params) 
{
	cv::FileStorage fs(filename, cv::FileStorage::READ);
	if (!fs.isOpened()) {
		return false;
	}
	fs["adaptiveThreshWinSizeMin"] >> params->adaptiveThreshWinSizeMin;
	fs["adaptiveThreshWinSizeMax"] >> params->adaptiveThreshWinSizeMax;
	fs["adaptiveThreshWinSizeStep"] >> params->adaptiveThreshWinSizeStep;
	fs["adaptiveThreshConstant"] >> params->adaptiveThreshConstant;
	fs["minMarkerPerimeterRate"] >> params->minMarkerPerimeterRate;
	fs["maxMarkerPerimeterRate"] >> params->maxMarkerPerimeterRate;
	fs["polygonalApproxAccuracyRate"] >> params->polygonalApproxAccuracyRate;
	fs["minCornerDistanceRate"] >> params->minCornerDistanceRate;
	fs["minDistanceToBorder"] >> params->minDistanceToBorder;
	fs["minMarkerDistanceRate"] >> params->minMarkerDistanceRate;
	fs["cornerRefinementMethod"] >> params->cornerRefinementMethod;
	fs["cornerRefinementWinSize"] >> params->cornerRefinementWinSize;
	fs["cornerRefinementMaxIterations"] >> params->cornerRefinementMaxIterations;
	fs["cornerRefinementMinAccuracy"] >> params->cornerRefinementMinAccuracy;
	fs["markerBorderBits"] >> params->markerBorderBits;
	fs["perspectiveRemovePixelPerCell"] >> params->perspectiveRemovePixelPerCell;
	fs["perspectiveRemoveIgnoredMarginPerCell"] >> params->perspectiveRemoveIgnoredMarginPerCell;
	fs["maxErroneousBitsInBorderRate"] >> params->maxErroneousBitsInBorderRate;
	fs["minOtsuStdDev"] >> params->minOtsuStdDev;
	fs["errorCorrectionRate"] >> params->errorCorrectionRate;
	return true;
}

int main(int argc, char **argv)
{
	CommandLineParser parser(argc, argv, keys);
	parser.about(about);

	if (argc < 2 || !parser.has("v") || !parser.has("c")) {
		parser.printMessage();
		return 0;
	}

	int dictionaryId = parser.get<int>("d");
	bool showRejected = parser.has("r");
	float markerLength = parser.get<float>("l");

	Ptr<aruco::DetectorParameters> detectorParams = aruco::DetectorParameters::create();
	if (parser.has("dp")) {
		bool readOk = readDetectorParameters(parser.get<string>("dp"), detectorParams);
		if (!readOk) {
			cerr << "Invalid detector parameters file" << endl;
			return 0;
		}
	}

	if (parser.has("refine")) {
		//override cornerRefinementMethod read from config file
		detectorParams->cornerRefinementMethod = parser.get<int>("refine");
	}
	std::cout << "Corner refinement method (0: None, 1: Subpixel, 2:contour, 3: AprilTag 2): " << detectorParams->cornerRefinementMethod << std::endl;

	String video = parser.get<String>("v");

	if (!parser.check()) {
		parser.printErrors();
		return 0;
	}

	Ptr<aruco::Dictionary> dictionary =
		aruco::getPredefinedDictionary(aruco::PREDEFINED_DICTIONARY_NAME(dictionaryId));

	Mat camMatrix, distCoeffs;
	bool readOk = readCameraParameters(parser.get<string>("c"), camMatrix, distCoeffs);
	if (!readOk) {
		cerr << "Invalid camera file" << endl;
		return 0;
	}

	VideoCapture inputVideo;
	inputVideo.open(video);

	std::vector<MotionMatrixd> trjs;

	while (inputVideo.grab()) {
		Mat image, imageCopy;
		inputVideo.retrieve(image);

		double tick = (double)getTickCount();

		vector< int > ids;
		vector< vector< Point2f > > corners, rejected;
		vector< Vec3d > rvecs, tvecs;

		// detect markers and estimate pose
		aruco::detectMarkers(image, dictionary, corners, ids, detectorParams, rejected);
		if (ids.size() > 0) {

			aruco::estimatePoseSingleMarkers(corners, markerLength, camMatrix, distCoeffs, rvecs,
				tvecs);
			MotionMatrixd tmp;
			cv::Mat _r;
			cv::Rodrigues(rvecs[0], _r);
			cv::cv2eigen(_r, tmp.R());
			cv::cv2eigen(tvecs[0], tmp.T());
			trjs.push_back(tmp);
		}
		// draw results
		image.copyTo(imageCopy);
		if (ids.size() > 0) {
			aruco::drawDetectedMarkers(imageCopy, corners, ids);

			for (unsigned int i = 0; i < ids.size(); i++)
				aruco::drawAxis(imageCopy, camMatrix, distCoeffs, rvecs[i], tvecs[i],
					markerLength * 0.5f);
		}

		if (showRejected && rejected.size() > 0)
			aruco::drawDetectedMarkers(imageCopy, rejected, noArray(), Scalar(100, 0, 255));

		imshow("out", imageCopy);
		waitKey(10);
	}
	// solve
	RevolutePair r_pair;
	EstError error;
	std::vector<MotionMatrixd> tmp_trjs;
	r_pair.Estimation(tmp_trjs, trjs, error);
	std::vector<MotionMatrixd> cor_trjs;
	r_pair.EstimationConsideringFreeDOF(cor_trjs, tmp_trjs, trjs, error.transDiv);

	inputVideo.open(video);

	// draw result
	// video analysis
	Eigen::Vector3d axis = r_pair.AxisDirection(RevolutePair::BASE);
	Eigen::Vector3d center = r_pair.CenterOfRotation(RevolutePair::BASE);
	std::cout << axis.transpose() << " " << center.transpose() << std::endl;
	cv::Mat p3d(2, 3, CV_32F);
	const float length = 0.20f;
	p3d.at<float>(0, 0) = static_cast<float>((center + length * axis)[0]);
	p3d.at<float>(0, 1) = static_cast<float>((center + length * axis)[1]);
	p3d.at<float>(0, 2) = static_cast<float>((center + length * axis)[2]);
	p3d.at<float>(1, 0) = static_cast<float>((center - length * axis)[0]);
	p3d.at<float>(1, 1) = static_cast<float>((center - length * axis)[1]);
	p3d.at<float>(1, 2) = static_cast<float>((center - length * axis)[2]);
	cv::Mat p2d;
	cv::projectPoints(p3d, cv::Vec3d::zeros(), cv::Vec3d::zeros(), camMatrix, distCoeffs, p2d);
	cv::Point p1(p2d.at<float>(0, 0), p2d.at<float>(0, 1));
	cv::Point p2(p2d.at<float>(1, 0), p2d.at<float>(1, 1));
	while (inputVideo.grab()) {
		cv::Mat image;
		inputVideo.retrieve(image);
		cv::line(image, p1, p2, CV_RGB(255, 0, 0), 3);
		cv::imshow("out", image);
		cv::waitKey(10);
	}
	return 0;
}
