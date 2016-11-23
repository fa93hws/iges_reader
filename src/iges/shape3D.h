#pragma once
#include <vector>
#include "../geo/point3d.h"
#include "../geo/vector3d.h"
#include "../nurbs/curve.h"
#include "../nurbs/surface.h"
#include <math.h>		//sin,cos,atan2
#include <algorithm>    // std::reverse
#include <armadillo>
#include <memory> //std::shared_ptr

NSI_BEG
/**************************************** Shape3D **************************************************/
class Shape3D {
public:
	arma::vec _T;
	arma::mat _R;
	/* The Transformation Matrix Entity transforms three - row column vectors by means of a matrix multiplication
	and then a vector addition.The notation for this transformation is :
	[R_11	R_12	R_13]	[X_in]	+	[T1]	=	[X_out]
	|R_21	R_22	R_23|	|Y_in|	+	|T2|	=	|Y_out|
	|R_31	R_32	R_33|	|Z_in|	+	|T3|	=	|Z_out|
	*/
	double _color[3];
	int _direction = 1;
	std::string _layer;
	int _level;
	std::string _name;

	

	Shape3D();
	void matrixTran(NS::Point3D &pt) const;
	void setMatTran(const arma::mat &R, const arma::vec T);
	void setDir(int dir) { _direction = dir; }
};
/**************************************** Surface3D **************************************************/
class Surface3D : public Shape3D {
protected:
	int _subClass;
	/*
		0 = Plane;
		1 = NURBS;
	*/
public:
	int getSurfaceTyp() const { return _subClass; }
};
/**************************************** Curve3D **************************************************/
class Curve3D : public Shape3D {
protected:
	int _subClass;
	/*
		0 = Line;
		1 = Circular Arc;
		2 = NURBS;
	*/
public:
};
/**************************************** NurbsSurf 128 ************************************************/
class Surf128 : public Surface3D {
	// p.156
private:
	NURBS::ParentSurface _surf;
public:
	void setSurf(const NURBS::ParentSurface &surf);
	void getSurf(NURBS::ParentSurface &ss) { ss = _surf; }
};
/**************************************** Edge 504 **************************************************/
class Edge504 {
	//p.546
private:
	int _curve;
	std::vector<int> _ver = std::vector<int>(2);
public:
	void setCurve(const int curve) { _curve = curve; }
	void setPts(const int pt0, const int pt1) { _ver[0] = pt0; _ver[1] = pt1; }
};
/**************************************** Loop 508 **************************************************/
class Loop508 {
	// p.548
private:
	std::vector<Edge504> _edge;
	std::vector<int> _direction;
public:
	size_t count() const { return _edge.size() - 1; }
	void addEdge(Edge504 const &edge, int const dir);
};
/**************************************** NurbsCurve 126 ************************************************/
class Curve126 : public Curve3D {
private:
	NURBS::Curve _curve;
public:
	void setCurve(const NURBS::ParentCurve &cc);
};
/**************************************** Face 510 **************************************************/
class Face510 {
	// p.550
private:
	int _surf;
	int _outerLoopFlag;
	std::vector<Loop508> _loops;
public:
	void setSurf(const int surf) { _surf = surf; }
	void setDir(const int dir) { _outerLoopFlag = dir; }
	void addLoop(const Loop508 &loop) { _loops.push_back(loop); }
};
/**************************************** Shell 514 **************************************************/
class Shell514 {
	// p.551
private:
	std::vector<Face510> _faces;
	std::vector<int> _direction;
public:
	void setLastDir(int dir) { _direction.back() = dir; };
	void setDir(int idx, int dir) { _direction[idx] = dir; };
	void addFace(const Face510 &face, const int dir);
};
/**************************************** Solid 186 **************************************************/
class Solid186 : public Shape3D {
	// p.229
private:
	Shell514 _shell;
	std::vector<Shell514> _voidShell;
public:

	void setShell(const Shell514 &shell) { _shell = shell; }
};


/**************************************** Group **************************************************/
template<class COMP,class DIR>
class ShapesGroup {
protected:
	std::vector<DIR> _directoryIdx;
	std::vector<COMP> _comp;
public:
	int count() const;
	void addComp(const COMP &comp, const DIR &idx) { _comp.push_back(comp); _directoryIdx.push_back(idx); }
	int findDirIdx(const DIR &dir);

	void getLast(COMP &comp) const { comp = _comp.back(); }
	void getComp(const int i, COMP &comp) const { comp = _comp[i]; }
	void getComps(std::vector<COMP> &comps) const { comps = _comp; }

};
/**************************************** Vertex 502 **************************************************/
class Vertex502 : public ShapesGroup<NS::Point3D, std::pair<int, int>> {
public:
	
};
/**************************************** Surfaces Group **************************************************/
class SurfacesGroup :public ShapesGroup<std::shared_ptr<Surface3D>,int>{
public:
	void addComp(const std::shared_ptr<Surface3D> &comp, const int &idx);
};
/**************************************** Curves Group **************************************************/
class CurvesGroup : public ShapesGroup<std::shared_ptr<Curve3D>,int> {
public:

};
/**************************************** Solids Group **************************************************/
class SolidsGroup : public ShapesGroup<Solid186,int>{
public:

};





NS_END