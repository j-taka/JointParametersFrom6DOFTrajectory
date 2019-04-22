/*! \brief  Estimate parameters of various joints from observation
 *           (If you know details, please see Section 5 in doctor dissertation of Jun Takamatsu) 
 *  \file   KinematicPair.h
 *  \author Jun Takamatsu <j-taka@is.naist.jp>
 *  \date   21 May 2004
 *  \date   17 Feb 2016 (re-implement)
 */
#pragma once

#include <MotionMatrix.h>
#include <vector>

typedef MotionMatrix<double> MotionMatrixd;

template<typename T>
class Polynomial;

/*! \brief accuracy of the estimation */
struct EstError
{
    double rotDiv;
    double transDiv;
	friend std::ostream& operator<<(std::ostream &o, const EstError &src);
};

/*! 
 *  \brief Base class of joint-parameter estimetor
 */
class KinematicPair
{
protected:
    std::vector<Eigen::Vector3d> param;
	double _NEARLY_ZERO;

public:
	enum TypeOfJoint { _FIXED, _REVOLUTE, _SPHERICAL, _PRISMATIC, _DOUBLE_REVOLUTE };
	static const int TARGET = 0;
	static const int BASE = 1;
	/*! \brief constructor */
    KinematicPair() : _NEARLY_ZERO(1.0e-6) {}
    /*! \brief destroctor */
	virtual ~KinematicPair(){}

    virtual void Estimation(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd> &src, EstError &dest2) = 0;
	virtual TypeOfJoint GetType() const = 0;
	virtual void Print(std::ostream &dest) const = 0;
	friend std::ostream& operator<<(std::ostream& dest, const KinematicPair &src){ 
		src.Print(dest); 
		return dest; 
	}
	static Eigen::Matrix3d RotXY(double src1, double src2);
	static Eigen::Matrix3d RotXZ(double src1, double src2);
    
protected:
    // Functions for estimating joint parameters in general
	// there is no DOF in orientation    
	void EstOri0Param(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd> &src, double &);
	// there is one DOF in orientation    
    void EstOri1Param(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd> &src, double &);
	// there is two DOF in orientation    
    void EstOri2Param(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd> &src, double &);
	// there is no constraint in orientation    
    void EstOri3Param(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd> &src, double &);
    void EstLoc0Param(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd> &src, double &);
    void EstOri0Loc1Param(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd > &src,double &);
    void EstOri1Hel1Param(std::vector<MotionMatrixd> &dest, const std::vector<MotionMatrixd > &src,double &);
	void EstLoc3Param(std::vector<MotionMatrixd > &dest, const std::vector<MotionMatrixd > &src, double &);
    
	/* 適当な座標系を設定する */
    static Eigen::Vector3d SetCoordinate(const Eigen::Vector3d &src);
	static void ZAxis2Coordinate(Eigen::Matrix3d &dest, const Eigen::Vector3d &src);
    
private:
	void InitOri2Param(std::vector<double> &dest, const std::vector<MotionMatrixd> &src) const;
	void SetEquationInInitOri2(std::vector<Polynomial<double> > &I, const std::vector<size_t> &samp, const std::vector<MotionMatrixd> &src) const;
	bool isValidSolutionInOri2(const std::vector<std::complex<double> > &src) const;
	void SetSolutionInOri2(std::vector<double> &dest, const Eigen::Vector3d &la, const Eigen::Vector3d &lb) const;
	bool CalculateAnotherAxis(Eigen::Vector3d &dest, const std::vector<Eigen::Vector3d> &src) const;
    void RefineOri0(std::vector<MotionMatrixd > &dest,
                    const std::vector<MotionMatrixd > &src, double &) const;
    void RefineOri1(std::vector<MotionMatrixd > &dest,
                    const std::vector<MotionMatrixd > &src, double &) const;
    void RefineOri2(std::vector<MotionMatrixd > &dest,
                    const std::vector<MotionMatrixd > &src, double &) const;
    void RefineOri3(std::vector<MotionMatrixd > &dest,
                    const std::vector<MotionMatrixd > &src, double &) const;
    void RefineOri0Loc1(std::vector<MotionMatrixd > &dest,
                        const std::vector<MotionMatrixd > &src, double &) const;
    void RefineLoc0(std::vector<MotionMatrixd > &dest,
                    const std::vector<MotionMatrixd > &src, double &) const;
    void RefineOri1Hel1(std::vector<MotionMatrixd > &dest, const std::vector<MotionMatrixd > &src,
                        const std::vector<double> &src2, double &) const;
    void RefineLoc3(std::vector<MotionMatrixd> &dest,
					const std::vector<MotionMatrixd> &src, double &) const;
    
    void EstOriAmount(std::vector<double> &dest,const std::vector<MotionMatrixd > &src) const;
	static void ChooseRandom(std::vector<size_t> &dest, size_t max_val);
};

#if 0
/* ネジ対偶用 */
class ScrewKP : public KinematicPair
{
  public:
    ScrewKP();
    void Estimation(std::vector<MotionMatrixd > &dest, const std::vector<MotionMatrixd > &src,
                    EstError &dest2);
    void Print(std::ostream &dest) const;
    Eigen::Vector3d Dic(int src) const; /* ネジの方向 */
    Eigen::Vector3d Cen(int src) const; /* ネジ座標系の中心 */
    Eigen::Vector3d Ori(int src) const; /* ネジ座標系の姿勢(rpy) */
    double Ratio() const; /* ネジ率(?) */
    friend std::ostream& operator<<(std::ostream &dest, const ScrewKP &src);
};

/* 姿勢方向１自由度のみの対偶 */
class Ori1DLoc3DKP : public KinematicPair
{
public:
	Ori1DLoc3DKP();
	void Estimation(std::vector<MotionMatrixd > &dest, const std::vector<MotionMatrixd > &src,
                    EstError &dest2);
    void Print(std::ostream &dest) const;
    Eigen::Vector3d Dic(int src) const; /* 使えない */
    Eigen::Vector3d Cen(int src) const; /* 回転中心 */
    Eigen::Vector3d Ori(int src) const; /* 使えない */
    double Ratio() const; /* 使えない */
    friend std::ostream& operator<<(std::ostream &dest,const Ori1DLoc3DKP &src);
};

/* 姿勢のみ固定 */
class Ori0DLoc3DKP: public KinematicPair
{
public:
	Ori0DLoc3DKP();
	void Estimation(std::vector<MotionMatrixd > &dest, const std::vector<MotionMatrixd > &src,
					EstError &dest2);
	void Print(std::ostream &dest) const;
	Eigen::Vector3d Dic(int src) const; /* 使えない */
	Eigen::Vector3d Cen(int src) const; /* 使えない */
	Eigen::Vector3d Ori(int src) const; /* 使えない */
	double Ratio() const; /* 使えない */
	friend std::ostream& operator<<(std::ostream &dest, const Ori0DLoc3DKP &src);
};

#endif











