//*************************************************************************
// File:     geometry.h
// Author:   Yan LIU
// Date:     09-10-2013
// Email:    liuyan0411@gmail.com
//*************************************************************************

#ifndef _GEOMETRY_H_09_10_2013__
#define _GEOMETRY_H_09_10_2013__



#include "..\shared\config.h"
#include "point3d.h"
#include "vector3d.h"
#include <cmath>


NS_BEGIN

class   Geometry
{
public:
	static const double PI;
	static const double EPS;
public:
	static double Volume6(const Point3D &a, const Point3D &b, const Point3D &c, const Point3D &d);

    //test if segment is intersected with triangle
    //if segment intersects triangle on one of triangle's edges
    //function also return true
    static int IsIntersected(const Point3D& tri0,const Point3D& tri1,const Point3D& tri2,
                              const Point3D& seg0,const Point3D& seg1,Point3D* intPt=NULL,double eps=1.0e-10);

	static bool TriSegIntersection(const Point3D& tri0, const Point3D& tri1, const Point3D& tri2,
                                   const Point3D& seg0, const Point3D& seg1, double eps,
							       Point3D& intPt, int& tl, int& sl);

	//return value: -1 no intersection, 0 intersected in the middle, 1 on seg0, 2 on seg1
	static int SegIntersectTriRelative(const Point3D& seg0, const Point3D& seg1,
		                               const Point3D& tri0, const Point3D& tri1, const Point3D& tri2, 
		                               double pEPS, double oEPS,
		                               Point3D& intPt);
	//return value: -1 no intersection, 0 intersected in the middle, 1 on seg0, 2 on seg1
	static int SegIntersectTriAbsolute(const Point3D& seg0, const Point3D& seg1,
		                               const Point3D& tri0, const Point3D& tri1, const Point3D& tri2, 
		                               /*double pEPS,*/ double segEPS, double triEPS,
									   Point3D& intPt);

	//return value: -1 no intersection, [0] inside, [1:4] node, [5:8] edge, 
	static int SquareIntersectSegAbsolute(const Point3D& sqr0, const Point3D& sqr1, const Point3D& sqr2, const Point3D& sqr3, 
		const Point3D& seg0, const Point3D& seg1,
		/*double parEPS,*/ double sqrEPS, double segEPS,
		Point3D& intPt);

	static Point3D ProjectPtOnLine(const Point3D& pt, const Point3D& ln0, const Point3D& ln1);
	static Point3D ProjectPtOnPlane(const Point3D& pt, const Point3D& pl0, const Point3D& pl1, const Point3D& pl2);
    
    static Vector3D Normal(const Point3D& a, const Point3D& b, const Point3D& c);
    static Vector3D UnitNormal(const Point3D& a, const Point3D& b, const Point3D& c);
	static Vector3D TriNormal(const Point3D& a, const Point3D& b, const Point3D& c);
	/*a and b must be normalized vector*/
	/*0-4, 0-2PI*/
	static double UnitVectorAngle(const Vector3D& a, const Vector3D& b, const Vector3D& axis);

	/*a and b must be normalized vector*/
	/*(-2,2], (-PI,PI]*/
	static double UnitVectorAngleBeta(const Vector3D& a, const Vector3D& b, const Vector3D& axis);

    static double Area(const Point3D& a, const Point3D& b, const Point3D& c);
    static double Area2(const Point3D& a, const Point3D& b, const Point3D& c);

	//segment
	enum SEG_INT{ExAB,ExBA,ExCD,ExDC,Parallel,Collinear,Intersect,IsPtA,IsPtB,IsPtC,IsPtD};
	static SEG_INT SegIntersection2D(const Point3D& pta, const Point3D& ptb, const Point3D& ptc, const Point3D& ptd,Point3D& interPoint);
	//ppos and qpos: 0 inside; 1 on end 0; 2 on end 1
	static bool    SegIntersection3D(const Point3D& p0, const Point3D& p1, 
			                         const Point3D& q0, const Point3D& q1,
									 double eps,        
									 Point3D& ipt, int& ppos, int& qpos);

	//the function is an approximation of closest points on two segments
	static void    SegClosestPoints(const Point3D& p0, const Point3D& p1, 
			                       const Point3D& q0, const Point3D& q1,
								   Point3D& ppt, Point3D& qpt);
	static bool    IsColinear(const Point3D& pt0, const Point3D& pt1, const Point3D& pt2, double cosang=0.9999); //default 1 degree

	//point
	static double  PointDistanceToSeg(const Point3D& pt, const Point3D& sa, const Point3D& sb);
	static void    PointClosestPtToSeg(const Point3D& pt, const Point3D& sa, const Point3D& sb, Point3D& cpt);
	static double  PointDistanceToLine(const Point3D& pt, const Point3D& la, const Point3D& lb);
	static void    PointClosestPtToTri(const Point3D& pt, const Point3D& t0, const Point3D& t1, const Point3D& t2, Point3D& cpt);
	static double  PointDistanceToTri(const Point3D& pt, const Point3D& t0, const Point3D& t1, const Point3D& t2);
	static double  PointDistanceToConvexPolygon(const Point3D& pt, const std::vector<Point3D>& polygon);
	static bool    IsPointInConvexPolygon(const Point3D& pt, const std::vector<Point3D>& polygon);

	//element quality
	//triangle quality r/max_edge
	static double TriQuality(const Point3D& a, const Point3D& b, const Point3D& c);

	//quadrialteral quality
	static double QuadQuality(const Point3D& a, const Point3D& b, const Point3D& c, const Point3D& d);

	//facet operations
	static Point3D   FacetCenter(const std::vector<Point3D>& nodes);
	static Vector3D  FacetNormal(const std::vector<Point3D>& nodes);
	static double    FacetSize(const std::vector<Point3D>& nodes);
	static double    FacetDPTriangulate(const std::vector<Point3D>& nodes, double cosang, const Point3D& sun, std::vector<int>& triangles);
};

NS_END

#endif
