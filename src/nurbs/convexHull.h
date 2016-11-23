#pragma once
#include <vector>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <algorithm>
#include "../geo/point3d.h"
#include "../geo/geometry.h"


namespace NURBS {
	class ConvexHull {
	private:
	public:
		double absTOL = 1e-7;
		struct vec3 {
			vec3() { X[0] = X[1] = X[2] = 0; }
			vec3(double x, double y, double z) { X[0] = x; X[1] = y; X[2] = z; }

			/* 3D cross product */
			vec3 operator*(const vec3& v) const {
				return vec3(X[1] * v.X[2] - X[2] * v.X[1],
					X[2] * v.X[0] - X[0] * v.X[2],
					X[0] * v.X[1] - X[1] * v.X[0]);
			}

			vec3 operator-(const vec3& v) const {
				return vec3(X[0] - v.X[0], X[1] - v.X[1], X[2] - v.X[2]);
			}

			vec3 operator-() const {
				return vec3(-X[0], -X[1], -X[2]);
			}

			double dot(const vec3& v) const {
				return X[0] * v.X[0] + X[1] * v.X[1] + X[2] * v.X[2];
			}

			double X[3];
		};
		struct twoset {
			void insert(int x) { (a == -1 ? a : b) = x; }
			bool contains(int x) { return a == x || b == x; }
			void erase(int x) { (a == x ? a : b) = -1; }
			int size() { return (a != -1) + (b != -1); }
			int a, b;
		};
		//E[1000][1000];

		struct face {
			vec3 norm;
			double disc;
			int I[3];
		};
		std::vector<face> faces;
		std::vector<vec3> A;
		twoset E[100][100];
		//std::vector<vector<twoset>> E;
		/* Compute the half plane {x : c^T norm < disc}
		* defined by the three points A[i], A[j], A[k] where
		* A[inside_i] is considered to be on the 'interior' side of the face. */
		face make_face(int i, int j, int k, int inside_i);

		ConvexHull(const std::vector<std::vector<NS::Point3D>> vertex);
		ConvexHull() {};
	};
}