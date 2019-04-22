// FixedPair.cpp

#include "FixedPair.h"

/*! \brief Estimate transformaiton matrix between two configurations
 *  \retval dest  relative configuration without noises
 *  \param  src   original relative configuration
 *  \param  dest2 estimation accuracy
 */
void FixedPair::Estimation(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd> &src, EstError &dest2)
{
	// clear data
	dest.clear();
    param.clear();
    // estimate orientation
    std::vector<MotionMatrixd> interTrj;
	EstOri0Param(interTrj, src, dest2.rotDiv);
    // estimate location
    EstLoc0Param(dest, interTrj, dest2. transDiv);
	param[1] += param[2];
	param.pop_back();
}

/*! \brief return translation
 *  \param src not use
 */
Eigen::Vector3d FixedPair::Translation(int src) const
{
	return param[1];
}

/*! \brief return rotation
 *  \param src not use
 */
Eigen::Vector3d FixedPair::Rotation(int src) const
{
	return param[0];
}

/*! \brief print estimation result (cannot use)
 *  \param o stream	
 */
void FixedPair::Print(std::ostream &o) const
{
	o << "Rotation: " << param[0] << "  Translation: " << param[1] << std::endl; 
}

std::ostream& operator<<(std::ostream &o,const FixedPair &src)
{
    src.Print(o);
    return o;
}
