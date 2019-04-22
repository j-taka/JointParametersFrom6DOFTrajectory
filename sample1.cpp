#include <iostream>
#include <iomanip>
#include <cmath>
#include "RevolutePair.h"
#include <ctime>
#include "s2rand.h"
#include <boost/random.hpp>
#include "eigen_matrix_utility.h"

// pcl
#include <boost/thread/thread.hpp>
#include <pcl/common/common_headers.h>
#include <pcl/visualization/pcl_visualizer.h>

const size_t NUM_OF_DATA = 100;

// global variables
boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer;
int count(0);
bool correct(false); 
std::vector<MotionMatrixd> trjs;
std::vector<MotionMatrixd> cor_trjs;

boost::shared_ptr<pcl::visualization::PCLVisualizer> InitializeViewer()
{
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer(new pcl::visualization::PCLVisualizer("Renderer"));
	viewer->setBackgroundColor(0.0, 0.0, 0.0);
	viewer->addCoordinateSystem(1.0);
	viewer->initCameraParameters();

	return (viewer);
}

void updateViewer()
{
	viewer->removeShape("line1");
	viewer->removeShape("line2");
	viewer->removeShape("line3");
	const double length = 0.5;
	Eigen::Vector3d cen, xd, yd, zd;
	if (correct){
		cen = cor_trjs[count].T();
		xd = length * cor_trjs[count].R() * Eigen::Vector3d::UnitX();
		yd = length * cor_trjs[count].R() * Eigen::Vector3d::UnitY();
		zd = length * cor_trjs[count].R() * Eigen::Vector3d::UnitZ();
	}
	else{
		cen = trjs[count].T();
		xd = length * trjs[count].R() * Eigen::Vector3d::UnitX();
		yd = length * trjs[count].R() * Eigen::Vector3d::UnitY();
		zd = length * trjs[count].R() * Eigen::Vector3d::UnitZ();
	}
	viewer->addLine(pcl::PointXYZ(cen[0], cen[1], cen[2]), pcl::PointXYZ(cen[0] + xd[0], cen[1] + xd[1], cen[2] + xd[2]), 1.0, 0.0, 0.0, "line1");
	viewer->addLine(pcl::PointXYZ(cen[0], cen[1], cen[2]), pcl::PointXYZ(cen[0] + yd[0], cen[1] + yd[1], cen[2] + yd[2]), 0.0, 1.0, 0.0, "line2");
	viewer->addLine(pcl::PointXYZ(cen[0], cen[1], cen[2]), pcl::PointXYZ(cen[0] + zd[0], cen[1] + zd[1], cen[2] + zd[2]), 0.0, 0.0, 1.0, "line3");
}

static void Evaluation(double &tr_err, double &rot_err, const std::vector<MotionMatrixd> &src, const std::vector<MotionMatrixd> &gt)
{
	assert(src.size() == gt.size());
	rot_err = 0;
	tr_err = 0;
	for (size_t i(0); i < gt.size(); ++i){
		double angle;
		Eigen::Vector3d axis;
		Mat2AngleAxis<double>(angle, axis, src[i].R().transpose() * gt[i].R());
		rot_err += angle;
		tr_err += (src[i].T() - gt[i].T()).norm(); 
	}
	rot_err /= gt.size();
	tr_err /= gt.size();
}

int main(int argc, char **argv)
{
	MotionMatrixd trA(Eigen::Vector3d(2, 3, 4), Eigen::Matrix3d::Identity());
	MotionMatrixd trB(Eigen::Vector3d(5, 4, 3), Eigen::Matrix3d::Identity());
	trB.R() = Eigen::AngleAxisd(M_PI / 2.0, Eigen::Vector3d::UnitX());
	// data making
	std::vector<MotionMatrixd> gt_motion(NUM_OF_DATA); 
	for (size_t i(0); i < gt_motion.size(); ++i){
		MotionMatrixd rot(Eigen::Vector3d::Zero(), Eigen::AngleAxisd(i * M_PI / NUM_OF_DATA, Eigen::Vector3d::UnitZ()).matrix());
		gt_motion[i] = trA * rot * trB;
	}
	// add noise
	// init rand
	srand((unsigned int) time(0));
	boost::random::mt19937 gen(static_cast<unsigned long>(time(0)));
	boost::random::normal_distribution<> dst(0, 1);
	boost::random::variate_generator<boost::random::mt19937, boost::random::normal_distribution<> > mt(gen, dst);
	trjs.resize(gt_motion.size());
	for (size_t i(0); i < gt_motion.size(); ++i){
		trjs[i] = gt_motion[i];
#if 1
		trjs[i].T()[0] += 0.1 * mt();
		trjs[i].T()[1] += 0.1 * mt();
		trjs[i].T()[2] += 0.1 * mt();
		double axis[3];
		::s2rand(axis);
		trjs[i].R() = trjs[i].R() * Eigen::AngleAxisd(mt() * 0.1, Eigen::Vector3d(axis[0], axis[1], axis[2]));
#endif
	}
	RevolutePair r_pair;
	EstError error;
	std::vector<MotionMatrixd> tmp_trjs;
	r_pair.Estimation(tmp_trjs, trjs, error);
	r_pair.EstimationConsideringFreeDOF(cor_trjs, tmp_trjs, trjs, error.transDiv);

	std::cout << "Axis: " << r_pair.AxisDirection(KinematicPair::BASE).transpose() << " (0 0 1) in base coordinates" << std::endl;	
	std::cout << "Axis: " << r_pair.AxisDirection(KinematicPair::TARGET).transpose() << " (0 1 0) in target coordinates" << std::endl;
	std::cout << "Center: " << r_pair.CenterOfRotation(KinematicPair::BASE).transpose() << " (2 3 3.5) in base coordinates" << std::endl;
	std::cout << "Center: " << r_pair.CenterOfRotation(KinematicPair::TARGET).transpose() << " (-5 -3.5 4) in target coordinates" << std::endl;
	// evaluation
	double tr_err(0), rot_err(0);
	Evaluation(tr_err, rot_err, trjs, gt_motion);
	std::cout << "Error in original - Orientation: " << rot_err << " " << " Translation: " << tr_err << std::endl;
	Evaluation(tr_err, rot_err, tmp_trjs, gt_motion);
	std::cout << "Error in intermediate - Orientation: " << rot_err << " " << " Translation: " << tr_err << std::endl;
	Evaluation(tr_err, rot_err, cor_trjs, gt_motion);
	std::cout << "Error in correction - Orientation: " << rot_err << " " << " Translation: " << tr_err << std::endl;
	viewer = InitializeViewer();
	while (!viewer->wasStopped()){
		updateViewer();
		count++;
		if (count == NUM_OF_DATA){
			correct = !correct;
			if (correct){
				std::cout << "Correct" << std::endl;
			}
			else{
				std::cout << "Original" << std::endl;
			}
			count = 0;
		}
		viewer->spinOnce();
		boost::this_thread::sleep(boost::posix_time::microseconds(10));
	}
	return 0;
}