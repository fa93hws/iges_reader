//*************************************************************************
// File:     vector3d.cpp
// Author:   Yan LIU
// Date:     04-10-2013
// Email:    liuyan0411@gmail.com
//*************************************************************************

#include "vector3d.h"
#include <math.h>

NS_BEGIN

double             Vector3D:: length() const{return sqrt(lengthSqrd());};


//bool  Vector3D::isEqualToVec   (const Vector3D& vec,double tol) const
//{
//    return angleTo(vec)<=tol;
//}

//bool Vector3D::isEqualToVecEx (const Vector3D& vec,double lineartol, double angleTol) const
//{
//    double len=this->length();
//    double lenVec=vec.length();
//    double diff=len-lenVec;
//    if(diff<0) diff=-diff;
//    if(len<=lineartol || lenVec<=lineartol)
//    {
//        return diff<=lineartol;
//    }else
//    {
//        return diff<=lineartol && angleTo(vec)<=angleTol;
//    }
//}


Vector3D Vector3D::crossProduct(const Vector3D& vec)const
{
    return Vector3D(_xyz[1]*vec[2]-_xyz[2]*vec[1],
        _xyz[2]*vec[0]-_xyz[0]*vec[2],
        _xyz[0]*vec[1]-_xyz[1]*vec[0]);
}

Vector3D  Vector3D::unitCrossProduct(const Vector3D& vector)const
{
    return crossProduct(vector).normalize();
}

//  return the angle of two vector, [0,PI]

//double Vector3D::angleTo(const Vector3D& vec) const
//{
//    double value=angleToCos(vec);
//    double angle=0.0;
//
//    if(fabs(fabs(value)-1.0)<=MTTolerance::EPS)
//    {
//        angle=value>0?0:M_PI;
//    }else
//    {		
//        angle=acos(value);
//    }
//    return angle;
//}

//double  Vector3D::angleToCos(const Vector3D& vec) const
//{
//    double dominator=module()*vec.module();
//
//    return dotProduct(vec)/dominator;
//}

//double  Vector3D::angleToCosEx(const Vector3D& vec,double& lenV1,double& lenV2) const
//{
//    lenV1=module();lenV2=vec.module();
//    double dominator=lenV1*lenV2;
//    if(dominator==0) dominator=MTTolerance::EPS;
//    return dotProduct(vec)/dominator;
//}

double Vector3D::cosine() const
{
    double l=module();
	if(l>0.0) return 1.0/l;
	else return 0.0;
}


//两个向量是否平行
bool Vector3D::isParallelTo(const Vector3D& vec,double tol ) const
{
    Vector3D vcx=norm().crossProduct(vec.norm());
    double dot=vcx.module();
    if(dot<tol) return true;
    return false;
}


//两个向量是否垂直
bool Vector3D::isPerpendicularTo(const Vector3D& vec,double tol) const
{
    Vector3D vcx=crossProduct(vec);
    double dot=vcx.module();
    double dot1=module();
    double dot2=vec.module();	
    if(fabs(dot-dot1*dot2)<tol) return true;
    return false;
}

//将this向量进行旋转
//绕轴x旋转, fTheta为弧度单位
Vector3D&  Vector3D::rotateByXaxis(double ang)
{	
    double cosin=cos(ang);
    double sinin=sin(ang);
    double yy = _xyz[1] * cosin - _xyz[2] * sinin;
    double zz = _xyz[1] * sinin + _xyz[2] * cosin;
    return set(_xyz[0],yy,zz);
}

//将this向量进行旋转
//绕轴y旋转, fTheta为弧度单位
Vector3D&  Vector3D::rotateByYaxis(double ang)
{
    double cosin=cos(ang);
    double sinin=sin(ang);
    double xx = _xyz[0] * cosin + _xyz[2] * sinin;
    double zz = -_xyz[0] * sinin + _xyz[2] * cosin;
    return set(xx,_xyz[1],zz);
}

//将this向量进行旋转
//绕轴z旋转, fTheta为弧度单位
Vector3D&  Vector3D::rotateByZaxis(double ang)
{
    double cosin=cos(ang);
    double sinin=sin(ang);
    double xx = _xyz[0] * cosin - _xyz[1] * sinin;
    double yy = _xyz[0] * sinin + _xyz[1] * cosin;
    return set(xx,yy,_xyz[2]);
}

//向量绕某轴旋转一个角度，可以计算某点绕一个轴旋转一个角度后的结果
Vector3D&   Vector3D::rotateBy(double ang , const Vector3D& axis)
{
    *this=rotate(ang,axis);
    return *this;
}

Vector3D   Vector3D::rotate(double ang , const Vector3D& axis)const
{
    double fx = axis[0],  fy = axis[1],  fz = axis[2];
    double c = cos(ang);
    double s = sin(ang);
    double xx = (fx * fx * (1.0f - c) + c)      * _xyz[0] +
        (fx * fy * (1.0f - c) - fz * s) * _xyz[1] +
        (fx * fz * (1.0f - c) + fy * s) * _xyz[2];
    double yy = (fy * fx * (1.0f - c) + fz * s) * _xyz[0] +
        (fy * fy * (1.0f - c) + c)      * _xyz[1] +
        (fy * fz * (1.0f - c) - fx * s) * _xyz[2];
    double zz = (fz * fx * (1.0f - c) - fy * s) * _xyz[0] + 
        (fz * fy * (1.0f - c) + fx * s) * _xyz[1] +
        (fz * fz * (1.0f - c) + c)      * _xyz[2];
    return Vector3D(xx,yy,zz);
}

//double  Vector3D::normalizeEx()
//{
//    double fLength=length();
//    if(fLength>MTTolerance::EPS)
//    {
//        double fInvLength = ((double)1.0)/fLength;
//        set(_xyz[0]*fInvLength,_xyz[1]*fInvLength,_xyz[2]*fInvLength);
//    }else
//    {
//        fLength = 0.0;
//        set(0.0,0.0,0.0);
//    }
//    return fLength;
//}


bool  Vector3D::isUnitLength(double tol) const
{
    double l=lengthSqrd();
    return fabs(l-1.0)<tol;
}

bool   Vector3D::isZeroLength(double tol) const
{
    if(fabs(_xyz[0])>tol) return false;
    if(fabs(_xyz[1])>tol) return false;
    if(fabs(_xyz[2])>tol) return false;
    return true;
}

bool   Vector3D::isCodirectionalTo  (const Vector3D& vec,double tol) const
{
    double eps=tol;
    double d0=cosine();
    double i0=_xyz[0]/d0;
    double j0=_xyz[1]/d0;
    double k0=_xyz[2]/d0;

    double d=vec.cosine();
    double ii=vec[0]/d;
    double jj=vec[1]/d;
    double kk=vec[2]/d;

    double dx=ii-i0;
    double dy=jj-j0;
    double dz=kk-k0;

    return (dx*dx+dy*dy+dz*dz)<eps;
}


//=========================================
Vector3D triangle_normal(const Point3D& pt1,
                           const Point3D& pt2,
                           const Point3D& pt3)
{
    //右手规则
    Vector3D vector1(pt1,pt2); 
    Vector3D vector2(pt1,pt3); 
    return vector1.crossProduct(vector2) ;
}




NS_END
