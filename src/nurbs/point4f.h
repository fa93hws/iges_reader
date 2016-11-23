#pragma once
#include "../geo/point3d.h"

namespace NURBS {
	class Point4f {
	private:
	public:
		double x, y, z, w;
		Point4f();
		Point4f(double a, double b, double c, double d);
		Point4f(const Point4f &P);
		Point4f(const NS::Point3D &P, const double weight);
		Point4f operator + (const Point4f &P1);
		Point4f operator - (const Point4f &P1);
		Point4f operator * (const double d);
		Point4f operator / (const double d);
		Point4f& operator /= (const double d);
	};
}