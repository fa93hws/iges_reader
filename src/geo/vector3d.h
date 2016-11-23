//*************************************************************************
// File:     vector3d.h
// Author:   Yan LIU
// Date:     04-10-2013
// Email:    liuyan0411@gmail.com
//*************************************************************************

#ifndef _VECTOR3D_H_04_10_2013__
#define _VECTOR3D_H_04_10_2013__


#include "..\shared\\config.h"
#include "point3d.h"


NS_BEGIN

class   Vector3D
{
public:
    Vector3D();
    Vector3D(const Point3D& pnt1,const Point3D& pnt2);

    Vector3D(double xx, double yy, double zz);
    Vector3D&  set(double xx,double yy,double zz);
    Vector3D&  set(const Point3D& pnt1,const Point3D& pnt2);

    Vector3D         operator /  (double scl) const;
    Vector3D&        operator /= (double scl);
    Vector3D&        operator *= (double scl);
    Vector3D&        operator += (const Vector3D& vec);
    Vector3D&        operator -= (const Vector3D& vec);
    Vector3D         operator -  ( ) const;
    Vector3D         operator - (const Vector3D& vec)const;
    Vector3D         operator + (const Vector3D& vec)const;
    double             operator [] (int i) const;
    double&            operator [] (int idx);


    //bool                 isEqualToVec   (const Vector3D& vec,double Angletol) const;
    //bool                 isEqualToVecEx (const Vector3D& vec,double lineartol, double angleTol) const;
    Vector3D&          negate();
    bool                 isParallelTo(const Vector3D& vec,double tol ) const;

    bool                 isPerpendicularTo(const Vector3D& vec,double tol) const;

    bool                 isCodirectionalTo  (const Vector3D& vec,double tol) const;

    Vector3D&          rotateBy(double ang , const Vector3D& axis);
    Vector3D           rotate(double ang , const Vector3D& axis)const;

    Vector3D&          rotateByXaxis(double ang);
    Vector3D&          rotateByYaxis(double ang);
    Vector3D&          rotateByZaxis(double ang);  

    double               module() const{return length();};
    double               length() const;
    double               len()const{return length(); }
    double               lengthSqrd()const;
    bool                 isUnitLength(double tol) const;
    bool                 isZeroLength(double tol) const;
    double               dotProduct(const Vector3D& vec) const;
    double               dotProduct(const Point3D& vec) const;
    Vector3D           crossProduct(const Vector3D& vector)const;
    Vector3D           unitCrossProduct(const Vector3D& vector)const;
    double               mixProduct(const Vector3D& veca,const Vector3D& vecb)const;

    //double               angleTo(const Vector3D& vec) const;
    //double               angleToCos(const Vector3D& vec) const;
    //double               angleToCosEx(const Vector3D& vec,double& lenV1,double& lenV2) const;

    double               cosine() const;


    Vector3D&         normalize();
    double              normalizeEx();
    Vector3D          norm()const;

    double              x() const{return _xyz[0];}
    double              y() const {return _xyz[1];}
    double              z() const {return _xyz[2];}
    const double*        memory()const{return _xyz.memory();}
    double*              memory(){return _xyz.memory();}
private:
    Point3D   _xyz;
};


 
Vector3D triangle_normal(const Point3D&,const Point3D&, const Point3D&);


//======================global function for 3-D vector =====================================
inline
std::ostream& operator<<(std::ostream& out,const Vector3D& vec)
{
    out<<"("<<vec[0]<<"\t"<<vec[1]<<"\t"<<vec[2]<<")";
    return out;
}

inline
double operator*(const Vector3D& vec1,const Vector3D& vec2)
{
    return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
}

inline
double operator*(const Point3D& vec1,const Vector3D& vec2)
{
    return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
}

inline
double operator*(const Vector3D& vec1,const Point3D& vec2)
{
    return vec2* vec1;
}

inline
double operator*(const Point3D& vec1,const Point3D& vec2)
{
    return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
}

inline
Vector3D crossProduct(const Vector3D& a, const Vector3D& b)
{
    return a.crossProduct(b);
}

inline
double dotProduct(const Vector3D& a, const Vector3D& b)
{
    return a*b;
}

inline 
Vector3D operator*(double val,const Vector3D& pt)
{
    return Vector3D(pt[0]*val, pt[1]*val, pt[2]*val);
}

inline 
Vector3D operator*(const Vector3D& pt,double val)
{
    return Vector3D(pt[0]*val, pt[1]*val, pt[2]*val);
}

inline
double   Vector3D::dotProduct(const Point3D& pt) const{ return *this*pt;}

inline
double  mixProduct(const Vector3D& veca,
                   const Vector3D& vecb,
                   const Vector3D& vecc)
{
    return vecc.dotProduct(veca.crossProduct(vecb));
}

//======================inline function for 3-D vector =====================================
inline
Vector3D::Vector3D()
{
    _xyz[0]=_xyz[1]=_xyz[2]=0.0;
}

inline
Vector3D::Vector3D(double xx, double yy, double zz)
{
    _xyz[0]=xx;
    _xyz[1]=yy;
    _xyz[2]=zz;
}

inline
Vector3D::Vector3D(const Point3D& pt1,const Point3D& pt2)
{
    _xyz[0]=pt2[0]-pt1[0];
    _xyz[1]=pt2[1]-pt1[1];
    _xyz[2]=pt2[2]-pt1[2];
}

inline Vector3D& 
Vector3D::set(const Point3D& pt1,const Point3D& pt2)
{
    return set(pt2[0]-pt1[0],pt2[1]-pt1[1],pt2[2]-pt1[2]);
}

inline double 
Vector3D::dotProduct(const Vector3D& vec) const
{
    return vec*(*this);
}


inline Vector3D&  
Vector3D::operator /= (double scale)
{
    *this*=(1.0/scale);
    return *this;
}

inline Vector3D
Vector3D::operator / (double val) const
{
    return Vector3D (_xyz[0]/val, _xyz[1]/val, _xyz[2]/val);
}


inline Vector3D&
Vector3D::operator *= (double val)
{  
    set(_xyz[0]*val,_xyz[1]*val,_xyz[2]*val);
    return *this;
}

inline Vector3D& 
Vector3D::operator += (const Vector3D& vec)
{
    set(_xyz[0]+vec[0],_xyz[1]+vec[1],_xyz[2]+vec[2]);
    return *this;
}

inline Vector3D& 
Vector3D::operator -= (const Vector3D& vec) 
{
    set(_xyz[0]-vec[0],_xyz[1]-vec[1],_xyz[2]-vec[2]);
    return *this;
}


inline double     
Vector3D::operator [] (int i) const
{
    return _xyz[i];
}

inline double&    
Vector3D::operator [] (int i)
{
    return _xyz[i];
}



inline Vector3D&   
Vector3D::set (double xx , double yy, double zz)
{
    _xyz.set(xx,yy,zz);
    return *this;
}

inline Vector3D&    Vector3D::negate()
{
    _xyz.set(-_xyz[0],-_xyz[1],-_xyz[2]);
    return *this;
}

inline double 
Vector3D::lengthSqrd() const
{
    return _xyz[0]*_xyz[0]+_xyz[1]*_xyz[1]+_xyz[2]*_xyz[2];
}

inline Vector3D
Vector3D::operator - () const
{
    return Vector3D (-_xyz[0], -_xyz[1], -_xyz[2]);
}

inline
Vector3D     Vector3D::operator - (const Vector3D& vec)const
{
    return Vector3D (_xyz[0]-vec[0], _xyz[1]-vec[1],_xyz[2]-vec[2]);
}

inline
Vector3D     Vector3D::operator + (const Vector3D& vec)const
{
    return Vector3D (_xyz[0]+vec[0], _xyz[1]+vec[1],_xyz[2]+vec[2]);    
}

inline
double  Vector3D::mixProduct(const Vector3D& veca,
                               const Vector3D& vecb)const
{
    return dotProduct(veca.crossProduct(vecb));
}

inline
Vector3D&  Vector3D::normalize()
{
    double d=cosine();
    return set(_xyz[0]*d,_xyz[1]*d,_xyz[2]*d); 
}


inline
Vector3D Vector3D::norm()const
{
    double d=cosine();
    return Vector3D(_xyz[0]*d,_xyz[1]*d,_xyz[2]*d);		
}


NS_END
#endif
