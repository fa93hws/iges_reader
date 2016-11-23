//*************************************************************************
// File:     point3d.cpp
// Author:   Yan LIU
// Date:     04-10-2013
// Email:    liuyan0411@gmail.com
//*************************************************************************



#include "point3d.h"
#include "vector3d.h"
#include "geometry.h"
#include <math.h>
#include <stdlib.h>
#include <memory.h>


NS_BEGIN

/**
@returns hypotenuse of Real (non-complex) scalars a and b by 
avoiding underflow/overflow
using (a * sqrt( 1 + (b/a) * (b/a))), rather than
sqrt(a*a + b*b).
*/
template <typename Real>
Real Hypot(const Real &a, const Real &b)
{
    if (a== 0)
        return fabs(b);
    else
    {
        Real c = b/a;
        return fabs(a) * sqrt(1 + c*c);
    }
}

template <typename Real>
Real Hypot(const Real &a, const Real &b,const Real &c)
{
    if (a== 0)
        return Hypot(b,c);
    else
    {
        Real d = b/a;
        Real e = c/a;
        return fabs(a) * sqrt(1 + d*d+e*e);
    }
}

inline 
bool IsEqual(double a,double b,double eps){return  (fabs(a-b)<eps);}

Point3D::operator Vector3D()const{return Vector3D(_xyz[0],_xyz[1],_xyz[2]);}


Point3D&    Point3D::operator += (const Vector3D& vec)
{
#ifdef   USING_MTGARRAY
    _xyz.detach();
#endif
    _xyz[0]+=vec[0];
    _xyz[1]+=vec[1];
    _xyz[2]+=vec[2];
    return *this;
}

int   Point3D::compare (const Point3D& rkV) const
{
    return memcmp(memory(),rkV.memory(),3*sizeof(double));
}


Point3D&    Point3D::operator -= (const Vector3D& vec)
{
#ifdef   USING_MTGARRAY
    _xyz.detach();
#endif
    _xyz[0]-=vec[0];
    _xyz[1]-=vec[1];
    _xyz[2]-=vec[2];
    return *this;
}



double  Point3D::distanceTo  (const Point3D& pnt) const
{
    double dx=_xyz[0]-pnt[0];
    double dy=_xyz[1]-pnt[1];
    double dz=_xyz[2]-pnt[2];
    return Hypot(dx,dy,dz);
}


double  Point3D::distanceToSquare(const Point3D& pnt)const
{
    double dx=_xyz[0]-pnt[0];
    double dy=_xyz[1]-pnt[1];
    double dz=_xyz[2]-pnt[2];
    return dx*dx+dy*dy+dz*dz;
}


bool  Point3D::isEqualTo(const Point3D& pnt,double tol) const
{
    if(!IsEqual(_xyz[0],pnt[0],tol)) return false;//
    if(!IsEqual(_xyz[1],pnt[1],tol)) return false;
    if(!IsEqual(_xyz[2],pnt[2],tol)) return false;
    return true;
}

double         Point3D::operator [] (int idx) const
{
    return  _xyz[idx];
}

double&       Point3D:: operator [] (int idx)
{
    return  _xyz[idx];
}

Point3D   Point3D::offset(const Vector3D& vec,double dis)const
{
    double len=vec.module();

    double dx=dis*vec[0]/len;
    double dy=dis*vec[1]/len;
    double dz=dis*vec[2]/len;
    return Point3D(_xyz[0]+dx,_xyz[1]+dy,_xyz[2]+dz);
}

Point3D&   Point3D::offsetilize(const Vector3D& vec,double dis)
{
    double len=vec.module();

    _xyz[0]+=dis*vec[0]/len;
    _xyz[1]+=dis*vec[1]/len;
    _xyz[2]+=dis*vec[2]/len;
    return *this;
}

Point3D   Point3D::rotate(double ang , const Vector3D& axis,
                              const Point3D& center)const

{
    Vector3D vec(center,*this);
    vec.rotateBy(ang,axis);
    return center+vec;
}



double TetVolume6  (const Point3D& a,
                    const Point3D& b,
                    const Point3D& c,
                    const Point3D& d)
{
    //return 
    //    ( - a[0] * b[1] * c[2] - a[1] * b[2] * c[0] - a[2] * b[0] * c[1] 
    //+ a[0] * b[2] * c[1] + a[1] * b[0] * c[2] + a[2] * b[1] * c[0] 
    //+ a[0] * b[1] * d[2] + a[1] * b[2] * d[0] + a[2] * b[0] * d[1] 
    //- a[0] * b[2] * d[1] - a[1] * b[0] * d[2] - a[2] * b[1] * d[0] 
    //- a[0] * c[1] * d[2] - a[1] * c[2] * d[0] - a[2] * c[0] * d[1] 
    //+ a[0] * c[2] * d[1] + a[1] * c[0] * d[2] + a[2] * c[1] * d[0] 
    //+ b[0] * c[1] * d[2] + b[1] * c[2] * d[0] + b[2] * c[0] * d[1] 
    //- b[0] * c[2] * d[1] - b[1] * c[0] * d[2] - b[2] * c[1] * d[0] );

     
    double rt = ( + a[0]*b[1]*c[2] - a[0]*b[1]*d[2] - a[0]*c[1]*b[2]
               + a[0]*c[1]*d[2] + a[0]*d[1]*b[2] - a[0]*d[1]*c[2]
               - b[0]*a[1]*c[2] + b[0]*a[1]*d[2] + b[0]*c[1]*a[2]
               - b[0]*c[1]*d[2] - b[0]*d[1]*a[2] + b[0]*d[1]*c[2]
               + c[0]*a[1]*b[2] - c[0]*a[1]*d[2] - c[0]*b[1]*a[2]
               + c[0]*b[1]*d[2] + c[0]*d[1]*a[2] - c[0]*d[1]*b[2]
               - d[0]*a[1]*b[2] + d[0]*a[1]*c[2] + d[0]*b[1]*a[2]
               - d[0]*b[1]*c[2] - d[0]*c[1]*a[2] + d[0]*c[1]*b[2]);
    return -rt;
        
}

Point3D Centroid(const Point3D& a,
                  const Point3D& b,
                  const Point3D& c,
                  const Point3D& d)
{
//	double vol6 = Geometry::Orient3D(a,b,c,d);
////	ASSERT(vol6 > 0.0);
	Vector3D av(d,a);
	Vector3D bv(d,b);
	Vector3D cv(d,c);
	Vector3D ccv = (av.dotProduct(av))*(bv.crossProduct(cv)) + 
				   (bv.dotProduct(bv))*(cv.crossProduct(av)) + 
				   (cv.dotProduct(cv))*(av.crossProduct(bv));
	ccv = ccv/(2.0*av.dotProduct(bv.crossProduct(cv)));
//	ccv = ccv/(2.0*vol6);
//	Point3D circumcenter = -1.0*(d+ccv);
	Point3D circumcenter = d+ccv;
#ifdef _SELFCHECK
	ASSERT(ISNAN(circumcenter[0]) == 0);
	ASSERT(ISNAN(circumcenter[1]) == 0);
	ASSERT(ISNAN(circumcenter[2]) == 0);
//	ASSERT(ISNAN(vol6) == 0);
#endif
	return circumcenter;
}


Vector3D  NormalTri(const Point3D& a,
                      const Point3D& b,
                      const Point3D& c)
{
    Vector3D vector1(a,b); 
    Vector3D vector2(a,c); 
    return vector1.crossProduct(vector2) ;
}


double  Area2(const Point3D& pta,const Point3D& ptb,const Point3D& ptc)
{
    Vector3D vec1(pta,ptb);
    Vector3D vec2(pta,ptc);
    return vec1.crossProduct(vec2).length();
}


bool   Coplanarity(const Point3D& a, 
                   const Point3D& b, 
                   const Point3D& c, 
                   const Point3D& d,
                   double eps)
{
    return  fabs(TetVolume6(a, b, c, d))<eps ;
}


NS_END
