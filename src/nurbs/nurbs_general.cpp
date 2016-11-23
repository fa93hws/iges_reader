#include "nurbs_general.h"

namespace NURBS {
	double TOL_ratio = 0.2; // use to find whether the point is in convex hull
	int NURBS_General::findContainingSurface(const std::vector<ParentSurface> &surfaces, const NS::Point3D &pt, std::vector<std::vector<int>> &index) {
		size_t nParentSurface = surfaces.size();
		int total = 0;
		for (size_t i = 0; i < nParentSurface; i++) {
			std::vector<int> temp;
			int nIdx = surfaces[i].containsPoint(pt,temp, TOL_ratio);
			if (nIdx > 0) {
				total += temp.size();
				index.push_back(temp);
				index.back().insert(index.back().begin(), i);
			}
		}
		if (total == 0) {
			std::vector<int> temp(2);
			total=1; double dist_cp = 9999.0;
			for (size_t i = 0; i < surfaces.size(); i++) {
				for (size_t j = 0; j < surfaces[i].children.size(); j++) {
					double average_dist =0.0;
					for (size_t k = 0; k < surfaces[i].children[j].controlPoints.size(); k++) {
						for (size_t w = 0; w < surfaces[i].children[j].controlPoints[k].size(); w++) {
							double dist_tmp = pt.distanceTo(surfaces[i].children[j].controlPoints[k][w]);
							average_dist += dist_tmp;
						}
					}
					average_dist = average_dist / (surfaces[i].children[j].controlPoints.size() * surfaces[i].children[j].controlPoints[0].size());
					if (average_dist < dist_cp) {
						temp[0] = i; temp[1] = j; dist_cp = average_dist;
					}
				}
			}
			index.push_back(temp);
		}
		return total;
	}

	void NURBS_General::findProjPt_Surface(const NS::Point3D& pt, std::vector<ParentSurface> &s, double &u, double &v, NS::Point3D& cpt) {
		std::vector<std::vector<int>> idx;
		int total = findContainingSurface(s, pt, idx);
		if (total == 1) {
			double u, v;
			s[idx[0][0]].children[idx[0][1]].surfacePointInv(pt, cpt, u, v);
		}
		else {
			int k = 0;
			std::vector<NS::Point3D> candidatePts(total);
			for (size_t i = 0; i < idx.size(); i++) {
				for (size_t j = 1; j < idx[i].size(); j++) {
					double u, v;
					s[idx[i][0]].children[idx[i][j]].surfacePointInv(pt, candidatePts[k], u, v);
					k++;
				}
			}
			double dist = pt.distanceTo(candidatePts[0]);
			int candidateIdx = 0;;
			for (int i = 1; i < total; i++) {
				if (pt.distanceTo(candidatePts[i]) < dist) {
					candidateIdx = i;
					dist = pt.distanceTo(candidatePts[i]);
				}
			}
			cpt = candidatePts[candidateIdx];
		}
	}
	
	bool NURBS_General::findLineIntersection(const NS::Point3D &nearPt, const NS::Point3D &seg0, const NS::Point3D &seg1,
		std::vector<ParentSurface> &s, const double TOL, double &u0, double &v0, NS::Point3D &intPt) {
		std::vector<std::vector<int>> idx;
		int total = findContainingSurface(s, nearPt, idx);
		if (total == 1) 
			return s[idx[0][0]].children[idx[0][1]].lineIntersection(nearPt, seg0, seg1, TOL, intPt);
		else {
			for (size_t i = 0; i < idx.size(); i++) {
				for (size_t j = 1; j < idx[i].size(); j++) {
					if (s[idx[i][0]].children[idx[i][j]].lineIntersection(nearPt, seg0, seg1, TOL, intPt))
						return true;
				}
			}
		}
		return false;
	}
}