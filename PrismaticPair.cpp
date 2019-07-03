// PrismaticPair.cpp
#include "PrismaticPair.h"
#include "eigen_matrix_utility.h"
#include <cassert>
#include <iostream>
#include <gsl/gsl_multimin.h>
#include <cstdlib>

int PrismaticPair::Save(const std::string &filename) const
{
	std::ofstream ofs(filename);
	if (!ofs) {
		return _FILENOTFOUND;
	}
	ofs << "Prismatic" << std::endl;
	SaveBase(ofs);
	return 0;
}

int PrismaticPair::Load(const std::string &filename)
{
	std::ifstream ifs(filename);
	if (!ifs) {
		return _FILENOTFOUND;
	}
	try {
		std::string str;
		ifs >> str;
		if (str != "Prismatic") {
			return _FORMATERROR;
		}
		LoadBase(ifs);
		return 0;
	}
	catch (...) {
		return _FORMATERROR;
	}
}

void PrismaticPair::Estimation(std::vector<MotionMatrixd > &dest,
					  const std::vector<MotionMatrixd > &src, EstError &dest2)
{
    std::vector<MotionMatrixd > interTrj;
    /* clear data */
    dest.clear();
    param.clear();
    /* estimate parameter in orientation */
    EstOri0Param(interTrj, src, dest2.rotDiv);
    /* estimate parameter in translation */
    EstOri0Loc1Param(dest, interTrj, dest2.transDiv);
}

/* ’¼i•ûŒü
 * @param src : 1=Šî€•¨‘Ì 0=‘ÎÛ•¨‘Ì */
Eigen::Vector3d PrismaticPair::TranslationDirection(int src) const
{
    assert(0 <= src && src < 2);
    return param[src + 1];
}

/* print parameter */
void PrismaticPair::Print(std::ostream &o) const
{
    o << "Target Object: " << param[0] << std::endl;
    o << "Base Object: " << param[1] << std::endl;
}

std::ostream& operator<<(std::ostream &o, const PrismaticPair &src)
{
    src.Print(o);
    return o;
}
