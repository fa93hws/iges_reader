//*************************************************************************
// File:     geometry.cpp
// Author:   Yan LIU
// Date:     09-10-2013
// Email:    liuyan0411@gmail.com
//*************************************************************************


#include "geometry.h"
#include "vector3d.h"
#include <cmath>

NS_BEGIN

const double Geometry::PI = 3.141592653589793238462643383279502884197169399375105820974944592307816406286;
const double Geometry::EPS = 1.0e-10;

double Geometry::Volume6(const Point3D &a, const Point3D &b, const Point3D &c, const Point3D &d)
{
	double tet6vol = ( - a[0] * b[1] * c[2] - a[1] * b[2] * c[0] - a[2] * b[0] * c[1] 
	                + a[0] * b[2] * c[1] + a[1] * b[0] * c[2] + a[2] * b[1] * c[0] 
	                + a[0] * b[1] * d[2] + a[1] * b[2] * d[0] + a[2] * b[0] * d[1] 
	                - a[0] * b[2] * d[1] - a[1] * b[0] * d[2] - a[2] * b[1] * d[0] 
	                - a[0] * c[1] * d[2] - a[1] * c[2] * d[0] - a[2] * c[0] * d[1] 
	                + a[0] * c[2] * d[1] + a[1] * c[0] * d[2] + a[2] * c[1] * d[0] 
	                + b[0] * c[1] * d[2] + b[1] * c[2] * d[0] + b[2] * c[0] * d[1] 
	                - b[0] * c[2] * d[1] - b[1] * c[0] * d[2] - b[2] * c[1] * d[0] );
	return tet6vol;
}

//the eps is relative eps
int Geometry::IsIntersected(const Point3D& tri0,const Point3D& tri1,const Point3D& tri2,
                                 const Point3D& seg0,const Point3D& seg1,Point3D* intPt,double eps)
{
    Vector3D e1,e2,p,s,q;
    double t,u,v,tmp;
    e1 = tri1 - tri0;
    e2 = tri2 - tri0;
    Vector3D dir = seg1 - seg0;
    p = dir.crossProduct(e2);
    tmp = p.dotProduct(e1);
    //test wether edge and triangle are parallel
	//Vector3D pv=p,pu=e1;
	//pv.normalize();
	//pu.normalize();
	//double pcos = fabs(pv.dotProduct(pu));
    //if(pcos< EPS)
    //{
    //    return -1;
    //}
	if (std::fabs(tmp) < EPS)
	{
		return -1;
	}
    tmp = 1.0/tmp;
    
    s = seg0 - tri0;
    u = tmp * (p.dotProduct(s));
    if(u < 0.0 || u > 1.0)
    {
        return -1;
    }
    
    q = s.crossProduct(e1);
    v = tmp * (q.dotProduct(dir));
    if(v < 0.0 || v > 1.0)
    {
        return -1;
    }
    
    if(u+v > 1.0)
    {
        return -1;
    }
    
    t = tmp * (q.dotProduct(e2));
    
    if(t < 0.0 || t > 1.0)
    {
        return -1;
    }
    
    double w = 1.0 - u - v;
    int rslt;
    if(-eps < w-1.0 && w-1.0 < eps)
    {
        //on the vertex v0
        rslt = 4;
    }else if(-eps < u-1.0 && u-1.0 < eps)
    {
        //on the vertex v1
        rslt = 5;
    }else if(-eps < v-1.0 && v-1.0 < eps)
    {
        //on the vertex v2
        rslt = 6;
    }else if(-eps < w && w < eps)
    {
        //on the edge v1-v2
        rslt = 2;
    }else if(-eps < u && u < eps) 
    {
        //on the edge v2-v0
        rslt = 3;
    }else if(-eps < v && v < eps) 
    {
        //on the edge v0-v1
        rslt = 1;
    }else
    {
        rslt = 0;
    }
    if(intPt != NULL)
    {
        *intPt = w*tri0+u*tri1+v*tri2;
    }
    return rslt;
}

//the eps is absolute eps
bool Geometry::TriSegIntersection(const Point3D& tri0, const Point3D& tri1, const Point3D& tri2,
                                  const Point3D& seg0, const Point3D& seg1, double eps,
							      Point3D& intPt, int& tl, int& sl)
{
	//check valid seg and tri
	double v_area = Area(tri0, tri1, tri2);
	if(v_area < EPS) return false;
	double v_dist = seg0.distanceTo(seg1);
	if(v_dist < EPS) return false;
	Vector3D tnor = UnitNormal(tri0, tri1, tri2);
	Vector3D sdir = Vector3D(seg0, seg1).normalize();
	double cosang = tnor.dotProduct(sdir);
	if(fabs(cosang) < EPS) return false;

	//calculate intersection 
	double p_s_dist = Vector3D(tri0, seg0).dotProduct(tnor);
	double s_dist   = p_s_dist/cosang;
	Point3D thept   = seg0 - sdir*s_dist;

	//check the segment
	double seg_0_dist = thept.distanceTo(seg0);
	double seg_1_dist = thept.distanceTo(seg1);
	sl = 0;
	if(seg_0_dist < seg_1_dist)
	{
		if(seg_0_dist < eps) 
		{
			sl = 1; //on side seg0
		}
	}else
	{
		if(seg_1_dist < eps) 
		{
			sl = 2; //on side seg1
		}
	}
	double seg_dist = seg_0_dist+seg_1_dist;
	if(sl==0 && (seg_dist-v_dist>EPS)) return false;

	//check the triangle
	tl = 0;
	double min_dist, dist; 

	double tri_area = Area(tri0, tri1, thept) + Area(tri1, tri2, thept) + Area(tri2, tri0, thept);
	if(tri_area-v_area > EPS) 
	{
		return false;
	}

	//node 0
	dist     = thept.distanceTo(tri0);
	min_dist = dist;
	if(min_dist < eps) tl = 4;

	//node 1
	dist     = thept.distanceTo(tri1);
	if(dist<min_dist)  
	{
		min_dist = dist;
		if(min_dist < eps) tl = 5;
	}

	//node 2
	dist     = thept.distanceTo(tri2);
	if(dist<min_dist)  
	{
		min_dist = dist;
		if(min_dist < eps) tl = 6;
	}

	if(tl != 0) 
	{
		intPt = thept;
		return true;
	}

	//edge 0->1
	dist     = PointDistanceToLine(thept, tri0, tri1);
	if(dist<min_dist)  
	{
		min_dist = dist;
		if(min_dist < eps) 
		{
			dist = (thept.distanceTo(tri0) + thept.distanceTo(tri1) - tri0.distanceTo(tri1))/2.0;
			if(dist < eps)
			{
				tl = 1;
			}
		}
	}

	//edge 1->2
	dist     = PointDistanceToLine(thept, tri1, tri2);
	if(dist<min_dist)  
	{
		min_dist = dist;
		if(min_dist < eps) 
		{
			dist = (thept.distanceTo(tri1) + thept.distanceTo(tri2) - tri1.distanceTo(tri2))/2.0;
			if(dist < eps)
			{
				tl = 2;
			}
		}
	}

	//edge 2->0
	dist     = PointDistanceToLine(thept, tri2, tri0);
	if(dist<min_dist)  
	{
		min_dist = dist;
		if(min_dist < eps) 
		{
			dist = (thept.distanceTo(tri2) + thept.distanceTo(tri0) - tri2.distanceTo(tri0))/2.0;
			if(dist < eps)
			{
				tl = 3;
			}
		}
	}

	if(tl != 0) 
	{
		//double check the tl
		int rt = IsIntersected(tri0,tri1,tri2,seg0,seg1,NULL,1.0e-2);
		if(rt==0) 
		{
			tl=0;
		} 
	}

	intPt = thept;
	return true;
}

//return value: -1 no intersection, 0 intersected in the middle, 1 on seg0, 2 on seg1
int Geometry::SegIntersectTriRelative(const Point3D& seg0, const Point3D& seg1,
	const Point3D& tri0, const Point3D& tri1, const Point3D& tri2,
	double pEPS, double oEPS,
	Point3D& intPt)
{
	Vector3D e1, e2, p, s, q;
	double t, u, v, tmp;
	e1 = tri1 - tri0;
	e2 = tri2 - tri0;
	Vector3D dir = seg1 - seg0;
	p = dir.crossProduct(e2);
	tmp = p.dotProduct(e1);
	//test wether edge and triangle are parallel
	//Vector3D pv = p, pu = e1;
	//pv.normalize();
	//pu.normalize();
	//double pcos = fabs(pv.dotProduct(pu));
	//if (pcos< pEPS)
	//{
	//	return -1;
	//}
	if (std::fabs(tmp) < pEPS)
	{
		return -1;
	}
	tmp = 1.0 / tmp;
	double zero = -oEPS;
	double one = 1.0 + oEPS;
	s = seg0 - tri0;
	u = tmp * (p.dotProduct(s));
	if (u < zero || u > one)
	{
		return -1;
	}

	q = s.crossProduct(e1);
	v = tmp * (q.dotProduct(dir));
	if (v < zero || v > one)
	{
		return -1;
	}

	if (u + v > one)
	{
		return -1;
	}

	t = tmp * (q.dotProduct(e2));

	if (t < zero || t > one)
	{
		return -1;
	}
	intPt = seg0 + dir*t;
	if (t < -zero)
	{
		return 1;
	}
	else if (t>(1.0 - oEPS))
	{
		return 2;
	}
	else
	{
		return 0;
	}
}

//return value: -1 no intersection, 0 intersected in the middle, 1 on seg0, 2 on seg1
int Geometry::SegIntersectTriAbsolute(const Point3D& seg0, const Point3D& seg1,
	const Point3D& tri0, const Point3D& tri1, const Point3D& tri2,
    /*double pEPS,*/ double segEPS, double triEPS,
	Point3D& intPt)
{
	////test wether edge and triangle are parallel
	//Vector3D tnor = UnitNormal(tri0, tri1, tri2);
	//Vector3D sdir = Vector3D(seg0, seg1).normalize();
	//double cosang = tnor.dotProduct(sdir);
	//if (std::fabs(cosang) < pEPS)
	//{
	//	return -1;
	//}

	Vector3D e1, e2, p, s, q;
	double t, u, v, w, tmp;
	e1 = tri1 - tri0;
	e2 = tri2 - tri0;
	Vector3D dir = seg1 - seg0;
	p = dir.crossProduct(e2);
	tmp = p.dotProduct(e1);
	if (std::fabs(tmp) < 1.0e-5) // the intersection is at infinity 
	{
		return -1;
	}
	//ASSERT(std::fabs(tmp) > parEPS);
	tmp = 1.0 / tmp;
	s = seg0 - tri0;
	q = s.crossProduct(e1);
	t = tmp * (q.dotProduct(e2));
	intPt = seg0 + dir*t;

	//check the segment
	double seg_0_dist = intPt.distanceTo(seg0);
	double seg_1_dist = intPt.distanceTo(seg1);
	if (seg_0_dist + seg_1_dist - seg0.distanceTo(seg1) > segEPS + segEPS)
	{
		return -1;
	}
	int rt = 0;
	if (seg_0_dist < segEPS)
	{
		rt = 1;
	}
	else if (seg_1_dist < segEPS)
	{
		rt = 2;
	}

	//check the triangle
	u = tmp * (p.dotProduct(s));
	v = tmp * (q.dotProduct(dir));
	w = 1.0 - u - v;
	if (0.0 < u && u < 1.0 && 0.0 < v && v < 1.0 && 0.0 < w && w < 1.0)
	{
		return rt;
	}

	double dist = 0.0;
	dist = PointDistanceToSeg(intPt, tri0, tri1);
	if (dist < triEPS)
	{
		return rt;
	}
	dist = PointDistanceToSeg(intPt, tri1, tri2);
	if (dist < triEPS)
	{
		return rt;
	}
	dist = PointDistanceToSeg(intPt, tri2, tri0);
	if (dist < triEPS)
	{
		return rt;
	}

	return -1;
}

//return value: -1 no intersection, [0] inside, [1:4] node, [5:8] edge, 
int Geometry::SquareIntersectSegAbsolute(const Point3D& sqr0, const Point3D& sqr1, const Point3D& sqr2, const Point3D& sqr3,
	const Point3D& seg0, const Point3D& seg1,
	/*double parEPS,*/ double sqrEPS, double segEPS,
	Point3D& intPt)
{
	////test wether edge and triangle are parallel
	//Vector3D tnor = UnitNormal(sqr0, sqr1, sqr2);
	//Vector3D sdir = Vector3D(seg0, seg1).normalize();
	//double cosang = tnor.dotProduct(sdir);
	//if (std::fabs(cosang) < parEPS)
	//{
	//	return -1;
	//}

	Vector3D e1, e2, p, s, q;
	double t, u01, u12, u23, u30, tmp;
	e1 = sqr1 - sqr0;
	e2 = sqr3 - sqr0;
	Vector3D dir = seg1 - seg0;
	p = dir.crossProduct(e2);
	tmp = p.dotProduct(e1);
	if (std::fabs(tmp) < 1.0e-5) // the intersection is at infinity 
	{
		return -1;
	}
	//ASSERT(std::fabs(tmp) > parEPS);
	tmp = 1.0 / tmp;
	s = seg0 - sqr0;
	q = s.crossProduct(e1);
	t = tmp * (q.dotProduct(e2));
	intPt = seg0 + dir*t;

	//check the segment
	double seg_0_dist = intPt.distanceTo(seg0);
	double seg_1_dist = intPt.distanceTo(seg1);
	if (seg_0_dist + seg_1_dist - seg0.distanceTo(seg1) > segEPS + segEPS)
	{
		return -1;
	}

	//check the square
	u30 = tmp * (p.dotProduct(s)) / 2.0; //edge sqr0-sqr3
	u01 = tmp * (q.dotProduct(dir)) / 2.0; //edge sqr0-sqr1

	e1 = sqr1 - sqr2;
	e2 = sqr3 - sqr2;
	dir = seg1 - seg0;
	p = dir.crossProduct(e2);
	tmp = p.dotProduct(e1);
	tmp = 1.0 / tmp;
	s = seg0 - sqr2;
	q = s.crossProduct(e1);

	u23 = tmp * (p.dotProduct(s)) / 2.0; //edge sqr2-sqr3
	u12 = tmp * (q.dotProduct(dir)) / 2.0; //edge sqr2-sqr1

	double dist = 0.0;
	Point3D pts[4] = { sqr0, sqr1, sqr2, sqr3 };
	//check vertex
	for (int i = 0;i < 4;i++)
	{
		dist = intPt.distanceTo(pts[i]);
		if (dist < sqrEPS)
		{
			return i + 1;
		}
	}
	//check edge
	for (int i = 0;i < 4;i++)
	{
		dist = PointDistanceToSeg(intPt, pts[i], pts[(i+1)%4]);
		if (dist < sqrEPS)
		{
			return i + 5;
		}
	}
	//check inner
	if (0.0 < u01 && u01 < 1.0 && 
		0.0 < u12 && u12 < 1.0 &&
		0.0 < u23 && u23 < 1.0 &&
		0.0 < u30 && u30 < 1.0)
	{
		return 0;
	}

	return -1;
}

Point3D Geometry::ProjectPtOnLine(const Point3D& pt, const Point3D& ln0, const Point3D& ln1)
{
	Vector3D vec(ln0, ln1);
	vec.normalize();
	double proj = vec.dotProduct(Vector3D(ln0, pt));
	return (ln0 + vec*proj);
}

Point3D Geometry::ProjectPtOnPlane(const Point3D& pt, const Point3D& pl0, const Point3D& pl1, const Point3D& pl2)
{
	double dis = pt.distanceTo(pl0);
	if(dis < 1.0e-10)
	{
		return pl0;
	}
	Vector3D fnorm = UnitNormal(pl0,pl1,pl2);
	double len = Vector3D(pl0, pt).dotProduct(fnorm);
	Point3D ppt = pt - len*fnorm;
	return ppt;
}

Vector3D Geometry::UnitNormal(const Point3D& a, const Point3D& b, const Point3D& c)
{
	Vector3D nor=Normal(a,b,c);
	nor.normalize();
	return nor;
}

Vector3D Geometry::Normal(const Point3D& a, const Point3D& b, const Point3D& c)
{
    Vector3D ab(a,b);
    Vector3D bc(b,c);
	Vector3D nor;
    nor = ab.crossProduct(bc);
	return nor;
}

/*a and b must be normalized vector*/
/*0-4, 0-2PI*/
double Geometry::UnitVectorAngle(const Vector3D& a, const Vector3D& b, const Vector3D& axis)
{
	double EPS = 1.0e-10;

	double dis = (b-a).lengthSqrd();
    /*check if the angle is 0*/
    if(dis < EPS) return 0.0;
    /*check if the angle is PI*/
    if(fabs(dis-4.0) < EPS) return 2.0;
	Vector3D dir = a.crossProduct(b);
	double cos = dir.dotProduct(axis);
    if(cos > 0.0)
    {
		cos = a.dotProduct(b);
        return 1.0-cos;
    }else
    {
        /*2.0+(1.0-(-uu).dotProduct(vv)); */
		cos = a.dotProduct(b);
        return 3.0+cos;
    }
}

/*a and b must be normalized vector*/
/*(-2,2], (-PI,PI]*/
double Geometry::UnitVectorAngleBeta(const Vector3D& a, const Vector3D& b, const Vector3D& axis)
{
	double EPS = 1.0e-10;

	double dis = (b - a).lengthSqrd();
	/*check if the angle is 0*/
	if (dis < EPS) return 0.0;
	/*check if the angle is PI*/
	if (fabs(dis - 4.0) < EPS) return 2.0;
	Vector3D dir = a.crossProduct(b);
	double cos = dir.dotProduct(axis);
	double sign = cos > 0.0 ? 1.0 : -1.0;
	cos = (1.0 - a.dotProduct(b))*sign;
	return cos;
}

double Geometry::Area(const Point3D& a, const Point3D& b, const Point3D& c)
{
	return Area2(a,b,c)/2.0;
}

double Geometry::Area2(const Point3D& a, const Point3D& b, const Point3D& c)
{
    Vector3D vec1(a,b);
    Vector3D vec2(a,c);
    return vec1.crossProduct(vec2).length();
}

Geometry::SEG_INT Geometry::SegIntersection2D(const Point3D& ptA, const Point3D& ptB, const Point3D& ptC, const Point3D& ptD,Point3D& interPoint)
{
	Point3D P0 = ptA;
	Point3D P1 = ptC;
	Vector3D D0(ptA,ptB);
	Vector3D D1(ptC,ptD);

	// segments P0 + s * D0 for s in [0, 1], P1 + t * D1 for t in [0,1] 
	Point3D E = P1 - P0; 
	double kross = D0.x() * D1.y() - D0.y() * D1.x(); 
	double sqrKross = kross * kross; 
	double sqrLen0 = D0.x() * D0.x() + D0.y() * D0.y(); 
	double sqrLen1 = D1.x() * D1.x() + D1.y() * D1.y(); 
	double sqrEpsilon = EPS;
	if (sqrKross > sqrEpsilon * sqrLen0 * sqrLen1) 
	{ 
		// lines of the segments are not parallel 
		double s = (E.x() * D1.y() - E.y() * D1.x()) / kross; 
		if(s < 0.0) return ExBA;
		if(s > 1.0) return ExAB;
		double t = (E.x() * D0.y() - E.y() * D0.x()) / kross; 
		if(t < 0.0) return ExDC;
		if(t > 1.0) return ExCD;
		// intersection of lines is a point on each segment 
		interPoint = P0 + s * D0; 
		return Intersect; 
	} 
	// lines of the segments are parallel 
	return Parallel;
}

namespace GeometryAux
{
//****************************************************************************80

double r8vec_dot_product ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DOT_PRODUCT computes the dot product of a pair of R8VEC's.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], A2[N], the two vectors to be considered.
//
//    Output, double R8VEC_DOT_PRODUCT, the dot product of the vectors.
//
{
	int i;
	double value;

	value = 0.0;
	for ( i = 0; i < n; i++ )
	{
		value = value + a1[i] * a2[i];
	}
	return value;
}

//****************************************************************************80

void lines_exp_near_3d (const double p1[3],const double p2[3],const double q1[3],const double q2[3], double pn[3], double qn[3] )

//****************************************************************************80
//
//  Purpose:
//
//    LINES_EXP_NEAR_3D computes nearest points on two explicit lines in 3D.
//
//  Discussion:
//
//    The explicit form of a line in 3D is:
//
//      the line through the points P1 and P2.
//
//    This routine uses a method that is essentially independent of dimension.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double P1[3], P2[3], two points on the first line.
//
//    Input, double Q1[3], Q2[3], two points on the second line.
//
//    Output, double PN[3], QN[3], the nearest points on the lines.
//
{
	double a;
	double b;
	double c;
	double d;
	double det;
	double e;
	int i;
	double sn;
	double tn;
	double u[3];
	double v[3];
	double w0[3];
	//
	//  Let U = (P2-P1) and V = (Q2-Q1) be the direction vectors on
	//  the two lines.
	//
	for ( i = 0; i < 3; i++ )
	{
		u[i] = p2[i] - p1[i];
	}
	for ( i = 0; i < 3; i++ )
	{
		v[i] = q2[i] - q1[i];
	}
	//
	//  Let SN be the unknown coordinate of the nearest point PN on line 1,
	//  so that PN = P(SN) = P1 + SN * (P2-P1).
	//
	//  Let TN be the unknown coordinate of the nearest point QN on line 2,
	//  so that QN = Q(TN) = Q1 + TN * (Q2-Q1).
	//
	//  Let W0 = (P1-Q1).
	//
	for ( i = 0; i < 3; i++ )
	{
		w0[i] = p1[i] - q1[i];
	}
	//
	//  The vector direction WC = P(SN) - Q(TC) is unique (among directions)
	//  perpendicular to both U and V, so
	//
	//    U dot WC = 0
	//    V dot WC = 0
	//
	//  or, equivalently:
	//
	//    U dot ( P1 + SN * (P2 - P1) - Q1 - TN * (Q2 - Q1) ) = 0
	//    V dot ( P1 + SN * (P2 - P1) - Q1 - TN * (Q2 - Q1) ) = 0
	//
	//  or, equivalently:
	//
	//    (u dot u ) * sn - (u dot v ) tc = -u * w0
	//    (v dot u ) * sn - (v dot v ) tc = -v * w0
	//
	//  or, equivalently:
	//
	//   ( a  -b ) * ( sn ) = ( -d )
	//   ( b  -c )   ( tc )   ( -e )
	//
	a = r8vec_dot_product ( 3, u, u );
	b = r8vec_dot_product ( 3, u, v );
	c = r8vec_dot_product ( 3, v, v );
	d = r8vec_dot_product ( 3, u, w0 );
	e = r8vec_dot_product ( 3, v, w0 );
	//
	//  Check the determinant.
	//
	det = - a * c + b * b;

	if ( det == 0.0 )
	{
		sn = 0.0;
		if ( std::fabs ( b ) < std::fabs ( c ) )
		{
			tn = e / c;
		}
		else
		{
			tn = d / b;
		}
	}
	else
	{
		sn = ( c * d - b * e ) / det;
		tn = ( b * d - a * e ) / det;
	}

	for ( i = 0; i < 3; i++ )
	{
		pn[i] = p1[i] + sn * ( p2[i] - p1[i] );
	}
	for ( i = 0; i < 3; i++ )
	{
		qn[i] = q1[i] + tn * ( q2[i] - q1[i] );
	}
	return;
}

}//GeometryAux


//ppos and qpos: 0 inside; 1 on end 0; 2 on end 1
bool    Geometry::SegIntersection3D(const Point3D& p0, const Point3D& p1, 
						  const Point3D& q0, const Point3D& q1,
						  double eps,        
						  Point3D& ipt, int& ppos, int& qpos)
{
	//calculate closest point on the two lines
	Point3D pn, qn;
	GeometryAux::lines_exp_near_3d ( p0.memory(), p1.memory(), q0.memory(), q1.memory(), pn.memory(), qn.memory());
	double dist,d0,d1;
	dist = pn.distanceTo(qn);
	if(dist > eps) return false;

	//check p 
	d0 = p0.distanceTo(pn);
	if(d0 < eps)
	{
		ppos = 1;
		pn   = p0;
	}else
	{
		d1 = p1.distanceTo(pn);
		if(d1 < eps)
		{
			ppos = 2;
			pn   = p1;
		}else
		{
			dist = p0.distanceTo(p1);
			if(std::fabs(d0+d1-dist) < (eps+eps))
			{
				ppos = 0;
			}else
			{
				return false;
			}
		}
	}

	//check q 
	d0 = q0.distanceTo(qn);
	if(d0 < eps)
	{
		qpos = 1;
		qn   = q0;
	}else
	{
		d1 = q1.distanceTo(qn);
		if(d1 < eps)
		{
			qpos = 2;
			qn   = q1;
		}else
		{
			dist = q0.distanceTo(q1);
			if(std::fabs(d0+d1-dist) < (eps+eps))
			{
				qpos = 0;
			}else
			{
				return false;
			}
		}
	}

	ipt = (pn+qn)/2.0;
	return true;
}

//the function is an approximation of closest points on two segments
void    Geometry::SegClosestPoints(const Point3D& p0, const Point3D& p1, 
							   const Point3D& q0, const Point3D& q1,
							   Point3D& ppt, Point3D& qpt)
{
	//calculate closest point on the two lines
	Point3D pn, qn;
	GeometryAux::lines_exp_near_3d ( p0.memory(), p1.memory(), q0.memory(), q1.memory(), pn.memory(), qn.memory());
	double dist, d0, d1, eps;

	//check p 
	d0   = p0.distanceTo(pn);
	d1   = p1.distanceTo(pn);
	dist = p0.distanceTo(p1);
	eps  = dist*1.0e-6;
	if(std::fabs(d0+d1-dist) < (eps+eps))
	{
		ppt = pn;
	}else
	{
		ppt = d0<d1?p0:p1;
	}

	//check q 
	d0   = q0.distanceTo(qn);
	d1   = q1.distanceTo(qn);
	dist = q0.distanceTo(q1);
	eps  = dist*1.0e-6;
	if(std::fabs(d0+d1-dist) < (eps+eps))
	{
		qpt = qn;
	}else
	{
		qpt = d0<d1?q0:q1;
	}
}


bool    Geometry::IsColinear(const Point3D& pt0, const Point3D& pt1, const Point3D& pt2, double cosang)
{
	Vector3D u(pt0,pt1);
	u.normalize();
	Vector3D v(pt1,pt2);
	v.normalize();
	double ang = u.dotProduct(v);
	ang = fabs(ang);
	if(ang > cosang)
	{
		return true;
	}else
	{
		return false;
	}
}

void Geometry::PointClosestPtToSeg(const Point3D& pt, const Point3D& sa, const Point3D& sb, Point3D& cpt)
{
	//calculate closest point
	Vector3D vec(sa,sb);
	double len = vec.length();
	vec /= len; //vec.normalize();
	double proj = vec.dotProduct(Vector3D(sa,pt));
	if(proj>=0.0 && proj<=len)
	{
		cpt= sa + proj*vec;
	}else
	{
		if(proj<0.0)
		{
			cpt= sa;
		}else
		{
			cpt= sb;
		}
	}
}

double  Geometry::PointDistanceToSeg(const Point3D& pt, const Point3D& sa, const Point3D& sb)
{
	Point3D cpt;
	PointClosestPtToSeg(pt, sa, sb, cpt);
	double dis = pt.distanceTo(cpt);
	return dis;
}

double  Geometry::PointDistanceToLine(const Point3D& pt, const Point3D& sa, const Point3D& sb)
{
	//calculate closest point
	Vector3D vec(sa,sb);
	vec.normalize();
	double proj = vec.dotProduct(Vector3D(sa,pt));
	Point3D close_pt = sa + proj*vec;

	double dis = pt.distanceTo(close_pt);
	return dis;
}

void Geometry::PointClosestPtToTri(const Point3D& pt, const Point3D& t0, const Point3D& t1, const Point3D& t2, Point3D& cpt)
{
	Point3D ppt = ProjectPtOnPlane(pt, t0, t1, t2);
	double areasub = 0.0;
	areasub += std::fabs(Area2(t0, t1, ppt));
	areasub += std::fabs(Area2(t1, t2, ppt));
	areasub += std::fabs(Area2(t2, t0, ppt));
	double areastd = std::fabs(Area2(t0, t1, t2));
	double eps = areastd*1.0e-3;
	if(std::fabs(areastd - areasub) < eps) //inside
	{
		cpt = ppt;
	}else
	{
		Point3D tpts[3] = {t0,t1,t2};
		double d,dis = 1.0e30;
		Point3D ppt;
		for(int i=0;i<3;i++)
		{
			PointClosestPtToSeg(pt, tpts[i], tpts[(i+1)%3], ppt);
			double d = pt.distanceToSquare(ppt);
			if(d<dis)
			{
				dis = d;
				cpt = ppt;
			}
		}
	}
}

double  Geometry::PointDistanceToTri(const Point3D& pt, const Point3D& t0, const Point3D& t1, const Point3D& t2)
{
	Point3D ppt = ProjectPtOnPlane(pt, t0, t1, t2);
	double areasub = 0.0;
	areasub += std::fabs(Area2(t0, t1, ppt));
	areasub += std::fabs(Area2(t1, t2, ppt));
	areasub += std::fabs(Area2(t2, t0, ppt));
	double areastd = std::fabs(Area2(t0, t1, t2));
	double eps = areastd*1.0e-3;
	if(std::fabs(areastd - areasub) < eps) //inside
	{
		double dis = pt.distanceTo(ppt);
		return dis;
	}else
	{
		double d,dis = 1.0e30;
		d = PointDistanceToSeg(pt, t0, t1);
		dis = d<dis?d:dis;
		d = PointDistanceToSeg(pt, t1, t2);
		dis = d<dis?d:dis;
		d = PointDistanceToSeg(pt, t2, t0);
		dis = d<dis?d:dis;
		return dis;
	}
}

double  Geometry::PointDistanceToConvexPolygon(const Point3D& pt, const std::vector<Point3D>& polygon)
{
	//assume the polygon is planar
	int psize = polygon.size();
	ASSERT(psize >= 3);
	Point3D ppt = ProjectPtOnPlane(pt, polygon[0], polygon[1], polygon[2]);
	bool inside = true;
	Vector3D u = Normal(polygon.back(), polygon.front(), ppt);
	for(int i=1;i<psize;i++)
	{
		Vector3D v = Normal(polygon[i-1], polygon[i], ppt); 
		double ang = u.dotProduct(v);
		if(ang < 0.0) 
		{
			inside = false;
			break;
		}
	}
	if(inside) //inside
	{
		double dis = pt.distanceTo(ppt);
		return dis;
	}else
	{
		double dis = PointDistanceToSeg(pt, polygon.back(), polygon.front());
		for(int i=1;i<psize;i++)
		{
			double d = PointDistanceToSeg(pt, polygon[i-1], polygon[i]);
			dis = d<dis?d:dis;
		}
		return dis;
	}
}

bool    Geometry::IsPointInConvexPolygon(const Point3D& pt, const std::vector<Point3D>& polygon)
{
	//assume the polygon is planar
	int psize = polygon.size();
	Point3D ppt = ProjectPtOnPlane(pt, polygon[0], polygon[1], polygon[2]);
	ASSERT(psize >= 3);
	bool inside = true;
	Vector3D u = Normal(polygon.back(), polygon.front(), ppt);
	for(int i=1;i<psize;i++)
	{
		Vector3D v = Normal(polygon[i-1], polygon[i], ppt); 
		double ang = u.dotProduct(v);
		if(ang < 0.0) 
		{
			inside = false;
			break;
		}
	}
	return inside;
}

//element quality r/max_edge
double Geometry::TriQuality(const Point3D& n0, const Point3D& n1, const Point3D& n2)
{
	double lens[3];
	lens[0] = n0.distanceTo(n1);
	lens[1] = n1.distanceTo(n2);
	lens[2] = n2.distanceTo(n0);
	double s = (lens[0] + lens[1] + lens[2]) / 2.0;
	double a = s - lens[0];
	if (a < 0.0) return 0.0;
	double b = s - lens[1];
	if (b < 0.0) return 0.0;
	double c = s - lens[2];
	if (c < 0.0) return 0.0;
	double r = sqrt((a*b*c) / s);
	ASSERT(ISNAN(r) == 0);
	double maxLen = Maximum3(lens[0], lens[1], lens[2]);
	ASSERT(ISNAN(maxLen) == 0);
	double q = r / maxLen;
	return q;
}

//quadrialteral quality
double Geometry::QuadQuality(const Point3D& a, const Point3D& b, const Point3D& c, const Point3D& d)
{
	double area2 = Area2(a, b, c);
	if (area2 < EPS) return 0.0;
	Vector3D norm =  UnitNormal(a, b, c);
	area2 = Area2(c, d, a);
	if (area2 < EPS) return 0.0;
	Vector3D normII =  UnitNormal(c, d, a);
	double dihedral = std::fabs(norm.dotProduct(normII));
//	if(dihedral < 0.999390827) //coplane check more than 2 degree
//	if(dihedral < 0.996194698) //coplane check more than 5 degree
	if(dihedral < 0.984807753012208) //coplane check more than 10 degree
	{
		return 0.0;
	}
	
	Point3D pnt[4] = { a,b,c,d };
	double maxang;
	double minlen, maxlen;
	Vector3D vec[4];
	vec[0] = Vector3D(a, b).normalize();
	vec[1] = Vector3D(b, c).normalize();
	vec[2] = Vector3D(c, d).normalize();
	vec[3] = Vector3D(d, a).normalize();
	maxang = UnitVectorAngle(vec[0], -vec[3], norm);
	minlen = d.distanceTo(a);
	maxlen = minlen;
	for (int i = 1;i < 4;i++)
	{
		double ang = UnitVectorAngle(vec[i], -vec[i-1], norm);
		maxang = ang > maxang ? ang : maxang;
		double len = pnt[i].distanceTo(pnt[i - 1]);
		minlen = len < minlen ? len : minlen;
		maxlen = len > maxlen ? len : maxlen;
	}
	if (maxang > 1.996) //>175
	{
		return 0.0;
	}
	double q = minlen / maxlen * (1.0 - std::fabs(1.0 - maxang));
	return q;
}


//facet operations
Point3D   Geometry::FacetCenter(const std::vector<Point3D>& nodes)
{
	int fnsize = nodes.size();
	ASSERT(fnsize > 2);
	Point3D center = Point3D(0.0, 0.0, 0.0);
	for (int i = 0;i < fnsize;i++)
	{
		center += nodes[i];
	}
	center = center / fnsize;
	return center;
}

Vector3D  Geometry::FacetNormal(const std::vector<Point3D>& nodes)
{
	double eps = 1.0e-5;
	int fnsize = nodes.size();
	ASSERT(fnsize > 2);


	std::vector<Vector3D> normals;
	normals.reserve(32);
	for (int i = 0;i < fnsize - 2;i++)
	{
		const Point3D& v0 = nodes[i];
		for (int j = i + 1;j < fnsize - 1;j++)
		{
			const Point3D& v1 = nodes[j];
			for (int k = j + 1;k < fnsize;k++)
			{
				const Point3D& v2 = nodes[k];
				double area = Geometry::Area2(v0, v1, v2);
				if (area > eps)
				{
					normals.push_back(Geometry::UnitNormal(v0, v1, v2));
				}
			}
		}
	}
	int nsize = normals.size();
	ASSERT(nsize > 0);
	Vector3D normal = Vector3D(0.0, 0.0, 0.0);
	for (int i = 0;i < nsize;i++)
	{
		normal += normals[i];
	}
	normal = normal / nsize;
	normal.normalize();

	return normal;
}

double    Geometry::FacetSize(const std::vector<Point3D>& nodes)
{
	int fnsize = nodes.size();
	ASSERT(fnsize > 2);
	double size = 0.0;
	for (int i = 0;i < fnsize;i++)
	{
		size += nodes[i].distanceTo(nodes[(i + 1) % fnsize]);
	}
	fnsize = fnsize / fnsize;
	return fnsize;
}

double    Geometry::FacetDPTriangulate(const std::vector<Point3D>& nodes, double cosang, const Point3D& sun, std::vector<int>& triangles)
{
	int size = nodes.size();
	ASSERT(size >= 3);
	std::vector<double> Q;
	std::vector<int>    K;
	double bad = -1.0e10;
	double flag = -1.0e9;//used to judge bad
	Q.resize(size*size, bad);
	K.resize(size*size, -1);
	for (int i = size - 3;i >= 0;i--)
	{
		const Point3D& vi = nodes[i];
		for (int j = i + 2;j<size;j++)
		{
			const Point3D& vj = nodes[j];
			for (int k = i + 1;k<j;k++)
			{
				const Point3D& vk = nodes[k];
				double q;
				Vector3D curnorm;
				double volume = Geometry::Volume6(vi, vk, vj, sun);
				if (volume<0.0)
				{
					q = bad;
				}
				else
				{
					q = Geometry::TriQuality(vi, vk, vj);
					curnorm = Geometry::UnitNormal(vi, vk, vj);
				}
				if (k<j - 1)
				{
					double qkj = Q[k*size + j];
					if (qkj<q) q = qkj;
					if (q > flag && qkj > flag)
					{
						int l = K[k*size + j];
						ASSERT(l >= 0);
						const Point3D& vl = nodes[l];
						Vector3D adjnorm = Geometry::UnitNormal(vk, vl, vj);
						double ang = Geometry::UnitVectorAngleBeta(adjnorm, curnorm, Vector3D(vk, vj));
						if (ang < cosang)
						{
							q = bad;
						}
					}
				}
				if (k>i + 1)
				{
					double qik = Q[i*size + k];
					if (qik<q) q = qik;
					if (q > flag && qik > flag)
					{
						int l = K[i*size + k];
						ASSERT(l >= 0);
						const Point3D& vl = nodes[l];
						Vector3D adjnorm = Geometry::UnitNormal(vi, vl, vk);
						double ang = Geometry::UnitVectorAngleBeta(adjnorm, curnorm, Vector3D(vi, vk));
						if (ang < cosang)
						{
							q = bad;
						}
					}
				}
				if (k == i + 1 || q>Q[i*size + j]) //k==i+1 just for initialization
				{
					Q[i*size + j] = q;
					K[i*size + j] = k;
				}
			}
		}
	}

	double q = Q[size - 1];
	if (q < flag)
	{
		return bad;
	}

	typedef std::pair<int, int> TEdge;
	std::stack<TEdge> edgestack;
	edgestack.push(TEdge(0, size - 1));
	triangles.reserve(size * 3);
	while (!edgestack.empty())
	{
		TEdge e = edgestack.top();
		edgestack.pop();
		if (e.second < e.first + 2) continue;
		ASSERT(e.first >= 0);
		ASSERT(e.second >= 0);
		int k = K[e.first*size + e.second];
		ASSERT(k >= 0);

		triangles.push_back(e.first);
		triangles.push_back(k);
		triangles.push_back(e.second);

		edgestack.push(TEdge(e.first, k));
		edgestack.push(TEdge(k, e.second));
	}
//	triangles.shrink_to_fit();
	triangles = triangles;

	return q;
}

NS_END
