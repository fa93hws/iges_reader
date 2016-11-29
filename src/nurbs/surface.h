#pragma once
#include "basis.h"
#include <vector>
#include "../geo/point3d.h"
#include "../geo/vector3d.h"
#include "point4f.h"
#include <stdexcept>
#include <assert.h>
#include "../geo/geometry.h"
#include "../algebra/matrix.h"
#include "../algebra/linearAlgebra.h"
#include "../algebra/polynomials.h"
#include "convexHull.h"
#include <armadillo>

NSN_BEG
class Surface {
private:
	std::vector<arma::mat> _mRep;
	int calCombiation(const int a, const int b) const;
	void calculateMRep();
	void getComp(const int i, const int j, std::vector<double> &comp) const;
	void surfaceDerivsAlg1(const std::vector<std::vector<NS::Point3D>> &P, const double u,
						const double v, const int d, std::vector<std::vector<NS::Vector3D>> &SKL);
	void calMrepVal(const NS::Point3D &pt, arma::mat &Mp) const;
	void calMrepVal(const NS::Vector3D &pt, arma::mat &Mp) const;
	void ratSurfaceDerives(const int d, const double u, const double v, 
						std::vector<std::vector<NS::Vector3D>> &SKL);
		
		
public:
	std::vector<Basis> basis;			// basis0-> xdir, basis1 ->ydir
	std::vector<std::vector<NS::Point3D>> controlPoints; // nbasis0 rows, nbasis1 column controlPoints[row][column]   
	std::vector<std::vector<Point4f>> pw;			// controlPoints with weight (4D)
	std::vector<std::vector<double>> weightVector;  // nbasis0 rows, nbasis1 column
	std::vector<std::vector<NS::Point3D>> referPoints;// used to determine the initial guess points
	bool is_plane = false;
	std::vector<double> plane_coe;						// plane : plane_coe[0] * x + plane_coe[1] * y + plane_coe[2] * z + plane_coe[3] = 0
	int numReferPoints = 11;							// num of refPoints
	double absTOL = 1e-7;
	std::string name;
	int id;												// id and name are used for any usage that required
	std::vector<ConvexHull::face> convexHull;
	bool rational = false;
	struct SurfaceInfo {
	public:
		std::vector<std::vector<NS::Point3D>> pnts;
		std::vector<std::vector<double>> knots;
		SurfaceInfo(std::vector<std::vector<NS::Point3D>> p, std::vector<std::vector<double>> k) { pnts = p; knots = k; }
		SurfaceInfo() { }
	};

	void getPw();
	void getCpFromPw();
	void getCoord(const double knot0, const double knot1, NS::Point3D &pt) const;
	void calculatePerImage(const NS::Point3D &pt,double &u, double &v) const;
	void calculateLineIntecs(const NS::Point3D &ptOnLine, const NS::Vector3D &slope,
		std::vector<double> &intecUs, std::vector<NS::Point3D> &intecPts) const;
	void projPt_initialGuessByLS(const NS::Point3D &pt, double &u0, double &u1);
	//void projPt_initialGuessByLS(const NS::Point3D &pt, double &u0, double &u1,NS::);
	void _getGridPoints(const int numPoints, SurfaceInfo &s);

	Surface(std::vector<Basis> b, std::vector<std::vector<NS::Point3D>> cP);
	Surface(std::vector<Basis> b, std::vector<std::vector<NS::Point3D>> cP,
		std::vector<std::vector<double>> w);
	Surface() {};
		
	//Advance alogirthm
	void surfaceKnotIns(const int dir, const double u, const int r);
	void getConvexHull();
	void surfacePointInv(const NS::Point3D &pt, NS::Point3D &projPt, double &u, double &v);
	void calDist(const NS::Point3D &pt, NS::Point3D &projPt);
	bool lineIntersection(const NS::Point3D &nearPt, const NS::Point3D &seg0, const NS::Point3D &seg1, const double TOL, double &u0, double &u1, NS::Point3D &intPt);
	bool lineIntersection(const NS::Point3D &nearPt, const NS::Point3D &seg0, const NS::Point3D &seg1, const double TOL, NS::Point3D &intPt);
	//Export
	void getScattors(const int numPoints, SurfaceInfo &s);
	void getScattors(SurfaceInfo &s);
	void offset(double dx, double dy, double dz);
	void discrete(const int n,std::vector<std::vector<int>> &faces, std::vector<NS::Point3D> &pts) const;
};

class ParentSurface : public Surface {
private:
public:
	ParentSurface() {};
	ParentSurface(std::vector<NURBS::Basis> b, std::vector<std::vector<NS::Point3D>> cP);
	ParentSurface(std::vector<NURBS::Basis> b, std::vector<std::vector<NS::Point3D>> cP,
					std::vector<std::vector<double>> w);
	std::vector<Surface> children;
	void subDivision();
	void offset(double dx, double dy, double dz);
	int containsPoint(const NS::Point3D &pt, std::vector<int> &idx,double TOL_ratio) const;
};
NS_END