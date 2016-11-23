#pragma once
// general nurbs functor
#include <vector>
#include "surface.h"
#include "curve.h"
#include "../geo/point3d.h"

namespace NURBS {
	class NURBS_General {
	private:
		
	public:
		static int findContainingSurface(const std::vector<ParentSurface> &surfaces, const NS::Point3D &pt, std::vector<std::vector<int>> &index);
		static void findProjPt_Surface(const NS::Point3D& pt, std::vector<ParentSurface> &s, double &u, double &v , NS::Point3D& cpt);
		static void findProjPt_Surface(const NS::Point3D& pt, std::vector<ParentSurface> &s, NS::Point3D& cpt) {
			double u, v;
			return findProjPt_Surface(pt, s,u, v, cpt);
		}

		static bool findLineIntersection(const NS::Point3D &nearPt, const NS::Point3D &seg0, const NS::Point3D &seg1,
			std::vector<ParentSurface> &s, const double TOL, double &u0, double &v0, NS::Point3D &intPt);
		static bool findLineIntersection(const NS::Point3D &nearPt, const NS::Point3D &seg0, const NS::Point3D &seg1,
			std::vector<ParentSurface> &s, const double TOL, NS::Point3D &intPt) {
			double u, v;
			return findLineIntersection(nearPt, seg0, seg1, s, TOL, u, v, intPt);
		}
	};
}