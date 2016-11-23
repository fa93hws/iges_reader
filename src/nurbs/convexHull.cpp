#include "convexHull.h"
namespace NURBS {
	ConvexHull::ConvexHull(const std::vector<std::vector<NS::Point3D>> vertex) {
		int N = 0;
		for (size_t i = 0; i < vertex.size(); i++) {
			for (size_t j = 0; j < vertex[0].size(); j++) {
				A.push_back(vec3(vertex[i][j].x(), vertex[i][j].y(), vertex[i][j].z()));
				N++;
			}
		}
		//for (int iijj = 0 ; iijj<N ; iijj++) {
		//for (int i = 0; i < N; i++) {
		//	cin >> A[i].X[0] >> A[i].X[1] >> A[i].X[2];
		//}
		//A.resize(7);
		//N = A.size();
		// find four points that make a tetrahedron whoes volumn is not zero
		std::vector<int> initList;
		bool tetraFlag = false;
		for (int i = 0; i < N; i++) {
			for (int j = i + 1; j < N; j++) {
				for (int k = j + 1; k < N; k++) {
					// these three points mush not be collinear
					double x0 = A[i].X[0], y0 = A[i].X[1], z0 = A[i].X[2];
					double x1 = A[j].X[0], y1 = A[j].X[1], z1 = A[j].X[2];
					double x2 = A[k].X[0], y2 = A[k].X[1], z2 = A[k].X[2];
					std::vector<double> n(3);
					n[0] = (y1 - y0) * (z2 - z1) - (z1 - z0) * (y2 - y1);
					n[1] = -(x1 - x0) * (z2 - z1) + (z1 - z0) * (x2 - x1);
					n[2] = (x1 - x0) * (y2 - y1) - (y1 - y0) * (x2 - x1);
					double norm = 0.0;
					for (double d : n)	norm += d * d;
					norm = sqrt(norm);
					if (abs(norm) < absTOL)
						continue;
					for (double &d : n) d /= norm;
					for (int ii = k + 1; ii < N; ii++) {
						// fourth point should not on the plane from plane(i,j,k)
						double D = 0.0;
						double distance = 0.0;
						for (int jj = 0; jj<3; jj++) {
							distance += n[jj] * A[ii].X[jj];
							D -= A[i].X[jj] * n[jj];
						}
						if (abs(distance + D) > absTOL) {
							tetraFlag = true;
							initList.push_back(i); initList.push_back(j); initList.push_back(k); initList.push_back(ii);
							break;
						}
					}
					if (tetraFlag) break;
				}
				if (tetraFlag) break;
			}
			if (tetraFlag) break;
		}
		if (initList.empty())
			return;
			//assert(0);
		/* Initially construct the hull as containing only the first four points. */
		face f;
		memset(E, -1, sizeof(E));

		for (int i = 0; i < 4; i++)
			for (int j = i + 1; j < 4; j++)
				for (int k = j + 1; k < 4; k++) {
					int rest_idx = initList[6 - i - j - k];
					faces.push_back(make_face(initList[i], initList[j], initList[k], initList[6 - i - j - k]));
				}

		/* Now add a point into the hull one at a time. */
		for (int i = 0; i < N; i++) {
			bool multiFlag = false;
			for (int ii = 0; ii < 4; ii++) {		// same in index
				if (initList[ii] == i) {
					multiFlag = true;
					break;
				}
				for (int ii = 0; ii < 4; ii++) {		//remove duplicated pts
					vec3 tmp = A[i] - A[initList[ii]];
					double dist2 = tmp.dot(tmp);
					if (dist2 < absTOL) {
						multiFlag = true;
						break;
					}
				}
				if (multiFlag)
					break;
			}

			for (int ii = 0; ii < i; ii++) {		// same in coordinate
				vec3 tmp = A[i] - A[ii];
				double dist2 = tmp.dot(tmp);
				if (dist2 < absTOL) {
					multiFlag = true;
					break;
				}
			}
			if (multiFlag)
				continue;
			//if (i == 6)
			//	int p = 0;

			/* Find and delete all faces with their outside 'illuminated' by this
			* point. */
			for (int j = 0; j < faces.size(); j++) {
				f = faces[j];
				double dist = f.norm.dot(A[i]);
				if (dist - f.disc > absTOL) { // not illuminated face
					E[f.I[0]][f.I[1]].erase(f.I[2]);
					E[f.I[0]][f.I[2]].erase(f.I[1]);
					E[f.I[1]][f.I[2]].erase(f.I[0]);
					faces[j--] = faces.back();
					faces.resize(faces.size() - 1);
				}
			}
			//if (i == 6)
			//	int p = 0;
			/* Now for any edge still in the hull that is only part of one face
			* add another face containing the new point and that edge to the hull. */
			int nfaces = faces.size();
			for (int j = 0; j < nfaces; j++) {
				f = faces[j];
				for (int a = 0; a < 3; a++)
					for (int b = a + 1; b < 3; b++) {
						int c = 3 - a - b;
						if (E[f.I[a]][f.I[b]].size() == 2)
							continue;
						if (abs(f.norm.dot(A[f.I[c]]) - f.disc) < absTOL) { // f.I[c] must not be on the face of f.I[a], f.I[b] , i
							bool findFlag = false;
							for (size_t iFace = 0; iFace < faces.size(); iFace++) {
								for (int iPt = 0; iPt < 3; iPt++)
									if (abs(f.norm.dot(A[faces[iFace].I[iPt]]) - f.disc) > absTOL) {
										faces.push_back(make_face(f.I[a], f.I[b], i, faces[iFace].I[iPt]));
										findFlag = true;
										break;
									}
								if (findFlag)
									break;
							}
						}
						else
							faces.push_back(make_face(f.I[a], f.I[b], i, f.I[c]));
					}
			}
		}
	}

	ConvexHull::face ConvexHull::make_face(int i, int j, int k, int inside_i) {
		E[i][j].insert(k); E[i][k].insert(j); E[j][k].insert(i);
		face f;
		f.I[0] = i; f.I[1] = j; f.I[2] = k;
		f.norm = (A[j] - A[i]) * (A[k] - A[i]);
		f.disc = f.norm.dot(A[i]); //-D
		if (f.norm.dot(A[inside_i]) > f.disc) {
			f.norm = -f.norm;
			f.disc = -f.disc;
		}
		return f;
	}
}