#include "point4f.h"

namespace NURBS {
	Point4f::Point4f() { x = 0.0; y = 0.0; z = 0.0; w = 0.0; }
	Point4f::Point4f(double a, double b, double c, double d)
	{
		x = a;
		y = b;
		z = c;
		w = d;
	}
	Point4f::Point4f(const Point4f &P)
	{
		w = P.w;
		x = P.x;
		y = P.y;
		z = P.z;
	}

	Point4f::Point4f(const NS::Point3D &P, const double d) {
		x = P.x() * d;
		y = P.y() * d;
		z = P.z() * d;
		w = d;
	}
	Point4f Point4f::operator +(const Point4f &P1)
	{
		return Point4f(P1.x + this->x, P1.y + this->y, P1.z + this->z, P1.w + this->w);
	}
	Point4f Point4f::operator -(const Point4f &P2)
	{
		return Point4f(this->x - P2.x, this->y - P2.y, this->z - P2.z, this->w - P2.w);
	}
	Point4f Point4f::operator *(double a)
	{
		Point4f R(a*this->x, a*this->y, a*this->z, a*this->w);
		return R;
	}
	Point4f Point4f::operator /(double a)
	{
		Point4f R(this->x / a, this->y / a, this->z / a, this->w / a);
		return R;
	}
	Point4f& Point4f::operator /=(double a){
		x /= a;
		y /= a;
		z /= a;
		w /= a;
		return *this;
	}
}