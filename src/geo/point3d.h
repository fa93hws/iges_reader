//*************************************************************************
// File:     point3d.h
// Author:   Yan LIU
// Date:     04-10-2013
// Email:    liuyan0411@gmail.com
//*************************************************************************

#ifndef _POINT3D_H_04_10_2013__
#define _POINT3D_H_04_10_2013__


#include "..\shared\config.h"
#include "..\shared\utility.h"
#include <iostream>

NS_BEGIN

class Vector3D;
class   Point3D
{
public:
    Point3D();
    Point3D(double xx, double yy, double zz);
    Point3D(const double xyz[3]);
    Point3D&    set(double xx,double yy,double zz);
    Point3D&    operator += (const Vector3D& vec);
    Point3D&    operator += (const Point3D& pt);
    Point3D&    operator -= (const Vector3D& vec);
    Point3D&    operator -= (const Point3D& pt);
    Point3D&    operator /= (double);
    Point3D&    operator *= (double);
    Point3D     operator -  ( ) const;
    Point3D     operator - (const Vector3D& vec)const;
    Point3D&    setToSum    (const Point3D& pnt, const Vector3D& vec);
    double        distanceTo  (const Point3D& pnt) const;
    double        distanceTo  (double pnt[3]) const;
    double        distanceToSquare(const Point3D& pnt)const;
    bool          isEqualTo      (const Point3D& pnt,double tol) const;
    double        operator [] (int i) const;
    double&       operator [] (int idx);
    const double* memory()const;
    double*       memory();
    Point3D     offset(const Vector3D& vec,double dis)const;
    Point3D&    offsetilize(const Vector3D& vec,double dis);
    Point3D&    negate();
    Point3D&    setToMin (const Point3D & p2);
    Point3D&    setToMax (const Point3D & p2);
    Point3D     rotate( double ang , const Vector3D& axis,const Point3D& center=Point3D(0,0,0))const;
    double        x() const{return _xyz[0];}
    double        y() const {return _xyz[1];}
    double        z() const {return _xyz[2];}
    operator      Vector3D()const;
    // comparison
    bool operator== (const Point3D& rkV) const;
    bool operator!= (const Point3D& rkV) const;
    bool operator<  (const Point3D& rkV) const;
    bool operator<= (const Point3D& rkV) const;
    bool operator>  (const Point3D& rkV) const;
    bool operator>= (const Point3D& rkV) const;
protected:
    int   compare (const Point3D& rkV) const;
#ifdef   USING_MTGARRAY
    MTMemArray<double> _xyz;
#else
    double     _xyz[3];
#endif
};

inline
std::ostream& operator<<(std::ostream& out,const Point3D& pt)
{
    out<<"("<<pt[0]<<"\t"<<pt[1]<<"\t"<<pt[2]<<")";
    return out;
}
 
double     Area2(const Point3D& pta,const Point3D& ptb,const Point3D& ptc);

 
double TetVolume6(const Point3D& a,
                  const Point3D& b,
                  const Point3D& c,
                  const Point3D& d);

Point3D Centroid(const Point3D& a,
                  const Point3D& b,
                  const Point3D& c,
                  const Point3D& d);

 
bool   Coplanarity(const Point3D& a, 
                   const Point3D& b, 
                   const Point3D& c, 
                   const Point3D& d,
                   double eps);

 
Vector3D  NormalTri(const Point3D& a,
                      const Point3D& b,
                      const Point3D& c);


inline
Point3D  operator + (const Point3D& pt,const Vector3D& vec)
{
    Point3D rt(pt);
    rt+=vec;
    return rt;
}

inline
Point3D  operator + (const Point3D& pt,const Point3D& vec)
{
    Point3D rt(pt);
    rt+=vec;
    return rt;
}

inline
Point3D  operator - (const Point3D& pt,const Point3D& vec)
{
    Point3D rt(pt);
    rt-=vec;
    return rt;
}

inline
Point3D::Point3D()
#ifdef   USING_MTGARRAY
:_xyz(3)
#endif
{
    _xyz[0]=_xyz[1]=_xyz[2]=0.0;
}

inline
Point3D::Point3D(const double xyz[3])
#ifdef   USING_MTGARRAY
:_xyz(3)
#endif
{
    _xyz[0]=xyz[0];
    _xyz[1]=xyz[1];
    _xyz[2]=xyz[2];
}

inline
Point3D::Point3D(double xx, double yy, double zz)
#ifdef   USING_MTGARRAY
:_xyz(3)
#endif
{
    _xyz[0]=xx;
    _xyz[1]=yy;
    _xyz[2]=zz;
}



inline 
Point3D& Point3D::setToMin (const Point3D & p2)
{
#ifdef   USING_MTGARRAY
    _xyz.detach();
#endif
    if (p2[0] < _xyz[0]) _xyz[0] = p2[0];
    if (p2[1] < _xyz[1]) _xyz[1] = p2[1];
    if (p2[2] < _xyz[2]) _xyz[2] = p2[2];
    return *this;
}

inline
Point3D       Point3D::operator -  ( ) const
{
    return  Point3D(-_xyz[0],-_xyz[1],-_xyz[2]);
}

inline  
Point3D & Point3D::setToMax (const Point3D & p2)
{
#ifdef   USING_MTGARRAY
    _xyz.detach();
#endif
    if (p2[0] > _xyz[0]) _xyz[0] = p2[0];
    if (p2[1] > _xyz[1]) _xyz[1] = p2[1];
    if (p2[2] > _xyz[2]) _xyz[2] = p2[2];
    return *this;
}

inline
const double*  Point3D::memory()const
{
#ifdef   USING_MTGARRAY
    return _xyz.data();
#else
    return _xyz;
#endif
}

inline
double*  Point3D::memory()
{
#ifdef   USING_MTGARRAY
    _xyz.detach();
    return _xyz.data();
#else
    return _xyz;
#endif
}


inline Point3D&   
Point3D::set (double xx , double yy, double zz)
{
#ifdef   USING_MTGARRAY
    _xyz.detach();
#endif
    _xyz[0]=xx;
    _xyz[1]=yy;
    _xyz[2]=zz;
    return *this;
}

inline
Point3D&    Point3D::operator += (const Point3D& pt)
{
#ifdef   USING_MTGARRAY
    _xyz.detach();
#endif
    _xyz[0]+=pt[0];
    _xyz[1]+=pt[1];
    _xyz[2]+=pt[2];
    return *this;
}

inline
Point3D&    Point3D::operator -= (const Point3D& pt)
{
#ifdef   USING_MTGARRAY
    _xyz.detach();
#endif
    _xyz[0]-=pt[0];
    _xyz[1]-=pt[1];
    _xyz[2]-=pt[2];
    return *this;
}

inline
Point3D&    Point3D::operator *= (double sc)
{
#ifdef   USING_MTGARRAY
    _xyz.detach();
#endif
    _xyz[0]*=sc;
    _xyz[1]*=sc;
    _xyz[2]*=sc;
    return *this;
}

inline
Point3D&    Point3D::operator /= (double sc)
{
    sc=1.0/sc;
    *this*=sc;
    return *this;
}



inline
Point3D      operator * (double sc, const Point3D& vec)
{
    return Point3D(sc*vec[0],sc*vec[1],sc*vec[2]);
}


inline
Point3D      operator * (const Point3D& vec,double sc)
{
    return Point3D(sc*vec[0],sc*vec[1],sc*vec[2]);
}

inline
Point3D      operator / (const Point3D& vec,double sc)
{
    sc=1.0/sc;
    return vec*sc;
}

inline
Point3D   Point3D::operator - (const Vector3D& vec)const
{
    Point3D rt(*this);
    rt-=vec;
    return rt;
}

inline
bool Point3D::operator== (const Point3D& rkV) const
{
    return compare(rkV) == 0;
}
//----------------------------------------------------------------------------
inline
bool Point3D::operator!= (const Point3D& rkV) const
{
    return compare(rkV) != 0;
}
//----------------------------------------------------------------------------
inline
bool Point3D::operator< (const Point3D& rkV) const
{
    return compare(rkV) < 0;
}
//----------------------------------------------------------------------------
inline
bool Point3D::operator<= (const Point3D& rkV) const
{
    return compare(rkV) <= 0;
}
//----------------------------------------------------------------------------
inline
bool Point3D::operator> (const Point3D& rkV) const
{
    return compare(rkV) > 0;
}
//----------------------------------------------------------------------------
inline
bool Point3D::operator>= (const Point3D& rkV) const
{
    return compare(rkV) >= 0;
}
//---------

NS_END
#endif
