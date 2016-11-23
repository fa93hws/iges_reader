#pragma once
#include "..\shared\config.h"
#include "basis.h"
#include <vector>
#include "point4f.h"
#include "../geo/point3d.h"
#include "../geo/vector3d.h"
#include "../algebra/polynomials.h"
#include <assert.h>
#include <stdlib.h>
#include <armadillo>

NSN_BEG
class Curve {
protected:
	std::vector<arma::mat> _mRep;
	void curveDerivsAlg1(const std::vector<NS::Point3D> &P, const double u, const int d, std::vector<NS::Vector3D> &SKL) const;
	int calCombiation(const int a, const int b) const;
	void getComp(const int i, std::vector<double> &comp) const;
	void calculateMRep(const int v);
	void calMrepVal(const NS::Point3D &pt,arma::mat &Mp) const;
	void calMrepVal(const NS::Vector3D &pt, arma::mat &Mp) const;
public:
	Basis basis;
	std::vector<double> weightVector;
	std::vector<NS::Point3D> controlPoints;
	std::string name;
	int id;
	bool rational = false;
	double absTOL = 1e-6;
	struct CurveInfo {
	public:
		std::vector<std::vector<NS::Point3D>> pnts;
		std::vector<double> knots;
		CurveInfo(std::vector<std::vector<NS::Point3D>> p, std::vector<double> k) { pnts = p; knots = k; }
		CurveInfo() { }
	};

	void getPw(const std::vector<NS::Point3D> &cp, std::vector<Point4f> &pw) const;
	void getCpFromPw(const std::vector<Point4f> &pw, std::vector<NS::Point3D> &cp, std::vector<double> &w) const;
	int getOrder() const { return basis.order; };
	void getMRep(std::vector<arma::mat> &mRep) const { mRep = _mRep; };

	void getCoord(const double kont, NS::Point3D &P) const;
	void ratCurveDerives(const int d, const double u, std::vector<NS::Vector3D> &SKL) const;
	void _getScatters(const int numPoints, CurveInfo &c);
	double arc_length( const double k0, const double k1) const;
	Curve() { };
	Curve(Basis &b, std::vector<NS::Point3D> &_cP);
	Curve(Basis &b, std::vector<NS::Point3D> &_cP, std::vector<double> &w);
	Curve::Curve(const Curve &c);
	void reverse();
	//advance algorithm
	void curveKnotIns(double u, int r);
	double calculatePerImage(const NS::Point3D &pt) const;
	void calculateLineIntecs(const NS::Point3D &ptOnLine, const NS::Vector3D &slope, 
		std::vector<double> &intecUs, std::vector<NS::Point3D> &intecPts) const;

	//io
	void getScattors(const int numPoints, CurveInfo &c);
	void getScattors(CurveInfo &c);
};

class ParentCurve : public Curve {
private:
	std::vector<Curve> children;
public:
	ParentCurve() {};
	// getset
	void getChildren(std::vector<Curve> &ch) const { ch = children; }
	// constructor
	ParentCurve(const NURBS::Curve _c);
	ParentCurve(const NURBS::Basis &b, const std::vector<NS::Point3D> &_cP);
	ParentCurve(const NURBS::Basis &b, const std::vector<NS::Point3D> &_cP, const std::vector<double> &w);
	//method
	void subDivision();

};
NS_END