#pragma once
#include <vector>
#include "../geo/point3d.h"
#include "../geo/vector3d.h"
#include "../nurbs/curve.h"
#include "../nurbs/surface.h"

NSI_BEG
	//class ShapeManger {
	//private:
	//	std::vector<Line> lines;									// type_id = 0
	//	std::vector<CircularArc> arcs;								// type_id = 1
	//	std::vector<IGES_NURBSCurve> n_curves;						// type_id = 2
	//	std::vector<IGES_NURBSSurface> n_surface;					// type_id = 3
	//public:
	//	void newShape(const Shape &s, int type_id);

	//};

	/***************************************** All Shapes Goes Here *******************************************/
class Shape {
protected:
	bool _propHierarchy=false;																	// use parent's prop
	std::string _name;
	std::string _layer;
	int _level;
	double _color[3];
	int _direction = 1;																		// -1 = reversed
	std::vector<std::vector<double>> R = { { 1,0,0 },{ 0,1,0 },{ 0,0,1 } };
	std::vector<double> T = { 0,0,0 };
	/* The Transformation Matrix Entity transforms three - row column vectors by means of a matrix multiplication
	and then a vector addition.The notation for this transformation is :
	[R_11	R_12	R_13]	[X_in]	+	[T1]	=	[X_out]
	|R_21	R_22	R_23|	|Y_in|	+	|T2|	=	|Y_out|
	|R_31	R_32	R_33|	|Z_in|	+	|T3|	=	|Z_out|
	*/
	void MatrixTransform(const std::vector<std::vector<double>> &R, const std::vector<double> &T, NS::Point3D &pt) const;
	double linear_interpolation(double seg0, double seg1, double pt) const;

public:																

	Shape() {};
	void set_tran_matrix(std::vector<std::vector<double>> &_R, std::vector<double> &_T);
	std::vector<std::vector<double>> get_R() const { return R; };
	std::vector<double> get_T() const { return T; };
	void set_name(const std::string &name) { _name = name; };
	std::string get_name() const { return _name; };
	void set_layer(const std::string &layer) { _layer = layer; };
	std::string get_layer()  const { return _layer; };
	void set_level(const int level) { _level = level; };
	int get_level() const { return _level; };
	void set_color(const double r, const double g, const double b) { _color[0] = r; _color[1] = g; _color[2] = b; }
	std::string get_color() const  {return (std::to_string(_color[0]) + "," + std::to_string(_color[1]) + "," + std::to_string(_color[2]));}
	void get_color(std::vector<double> &color) const { color.clear(); for (double d : _color) color.push_back(d); };
	void set_dir(const int dir)   { _direction = dir; }
	int get_dir() const { return _direction; }
	void copyAttrFrom(const Shape &s) { 
		_layer = s.get_layer(); 
		std::vector<double> color;
		s.get_color(color);
		set_color(color[0], color[1], color[2]);
	}
	void setpropHierarchy(bool hi) { _propHierarchy = hi; }
	bool getpropHierarchy() const { return _propHierarchy; }

	virtual void discrete(const int num_pt, std::vector<NS::Point3D> &pts) =0;
	virtual void discrete(const double TOL, std::vector<NS::Point3D> &pts) =0;
	virtual void get_pt_from_para(const double t, NS::Point3D &pt) const =0;			// parameter space from -1 to 1
};

/***************************************** Circlular arc-100 *******************************************/
class CircularArc :public Shape {
	// Anticlockwise
private:
	std::vector<NS::Point3D> discretePts;
	NS::Point3D cnt;						// [cnt_x,cnt_y,z_t] where ZT is the coordinate of a point along the ZT axis
	std::vector<double> theta;				// [theta1,theta2],arc starts from theta1, end theta2
	double radius;
	bool _isClosed = false;
public:
	CircularArc() {};
	CircularArc(NS::Point3D &pt1, NS::Point3D &pt2, NS::Point3D &pt3);
	CircularArc(const NS::Point3D &_cnt, const double radius, const double theta0, const double theta1);
	void discrete(const int num_pt, std::vector<NS::Point3D> &pts) ;
	void discrete(const double num_pt, std::vector<NS::Point3D> &pts);
	void get_pt_from_para(const double t, NS::Point3D &pt) const;
	NS::Point3D getCnt() const { return cnt; }
	double getRadius() const { return radius; }
	void getTheta(std::vector<double> &t) const { t = theta; }
	bool isClosed() const { return _isClosed; }
	void set_tran_matrix(std::vector<std::vector<double>> &_R, std::vector<double> &_T);

	bool isConincideWithArc(const CircularArc &arc, const double TOL) const;
};


/***************************************** Line-110 *******************************************/
class Line :public Shape {
private:
	NS::Point3D pt0;
	NS::Point3D pt1;
public:
	int form = 0;							// 0 = line segement, 1=semi-bounded line, 2=unbounded line
	Line(NS::Point3D &_pt0, NS::Point3D &_pt1) { pt0 = _pt0; pt1 = _pt1; }
	Line() {}
	void discrete(const int num_pt, std::vector<NS::Point3D> &pts) ;
	void discrete(const double TOL, std::vector<NS::Point3D> &pts)  { discrete(2, pts); };
	void get_pt_from_para(const double t, NS::Point3D &pt) const ;
	void set_tran_matrix(std::vector<std::vector<double>> &_R, std::vector<double> &_T);

	bool isConincideWithLine(const Line &l,const double TOL) const;
};



/***************************************** NURBS Curve - 126 *******************************************/
class IGES_NURBSCurve : public Shape{
private:
	NURBS::Curve c;
	std::vector<NS::Point3D> discretePts;
	double polyline_length(const NS::Point3D &pt0, const NS::Point3D &pt1, const NS::Point3D &pt2) const;
	double arc_len_pt_loc_convex(const double u0, const double u1) const;
	double arc_len_pt_loc_concave1(const double u0, const double u1,const double ui) const;
	double arc_len_pt_loc_general(const double u0, const double u1, const double arc_len) const;
	double arc_len_pt_loc(const double u0, const double u1, const double arc_len) const;
		
public:
	IGES_NURBSCurve() { };
	IGES_NURBSCurve(NURBS::Basis &b, std::vector<NS::Point3D> &_cP) { c = NURBS::Curve(b, _cP); }
	IGES_NURBSCurve(NURBS::Basis &b, std::vector<NS::Point3D> &_cP, std::vector<double> &w) { c = NURBS::Curve(b, _cP, w); };
	IGES_NURBSCurve( NURBS::Curve &_c ){c = _c; };
	void discrete(const int num_pt, std::vector<NS::Point3D> &pts) ;
	void discrete(const double TOL, std::vector<NS::Point3D> &pts) ;
	void get_pt_from_para(const double t, NS::Point3D &pt) const;
	void get_curve(NURBS::Curve &_c) const { _c = c; }
	int isConincideWith(const IGES_NURBSCurve &_c, const double TOL) const;
		
};
/***************************************** NURBS Surface - ??? *******************************************/
class IGES_NURBSSurface : public NURBS::Surface{

};
/***************************************** Shape Group *******************************************/
class ShapeGroup : public Shape {
protected:
	std::vector<Line> lines;									// type_id = 0
	std::vector<CircularArc> arcs;								// type_id = 1
	std::vector<IGES_NURBSCurve> curves;						// type_id = 2
	std::vector<int> _idx;										// [0,1,1,2] = first ->lines; second -> arcs; third -> arcs; fourth -> nurbs curve
public:
	ShapeGroup() {}
	void add_line(const Line &line);
	void add_arc(const CircularArc &arc);
	void add_nurbs_curves(const IGES_NURBSCurve &curve);
	void discrete(const int num_pt, std::vector<NS::Point3D> &pts) ;
	void discrete(const double TOL, std::vector<NS::Point3D> &pts) ;
	void discrete(const double TOL, std::vector<std::vector<NS::Point3D>> &pts) ;
	void discrete(const double TOL, std::vector<std::list<int>> &sequence) ;
	void discreteByIdx(const int idx, const double TOL, std::vector<NS::Point3D> &pts,  std::string &layer, std::vector<double> &color) ;
	void get_pt_from_para(const double t, NS::Point3D &pt) const;
	void find_dir(const NS::Point3D &pt);

	void getLine(const int idx, Line &l) const { l = lines[idx]; };
	void getArc(const int idx, CircularArc &l) const { l = arcs[idx]; };
	void getNURBSCurve(const int idx, IGES_NURBSCurve &l) const { l = curves[idx]; };
	std::vector<int> getIdx() { return _idx; }
		

	int getLineSize() const { return lines.size(); };
	int getArcSize() const { return arcs.size(); };
	int getNURBSCurveSize() const { return curves.size(); };
};
/***************************************** Surface-2D *******************************************/
class Surface2D : public ShapeGroup {
private:
	std::vector<ShapeGroup> _sg;								// type_id = 3
public:
	Surface2D() {}
	void add_sg(const ShapeGroup &_sg);
	void discrete(const double TOL, std::vector<std::vector<NS::Point3D>> &pts) ;
};



/***************************************** Shape-Manger *******************************************/
class ShapeManger { // not in use
protected:
	std::vector<Line> _lines;									// type_id = 0
	std::vector<CircularArc> _arcs;								// type_id = 1
	std::vector<IGES_NURBSCurve> _curves;						// type_id = 2
	std::vector<IGES_NURBSSurface> _surface;					// type_id = 3
	std::vector<Surface2D> _surface2d;							// type_id = 4
	std::vector<ShapeGroup> _group;								// type_id = 5
public:
	void addLine(Line &l) { _lines.push_back(l); };
	void addArc(CircularArc &a) { _arcs.push_back(a); };
	void addNURBSCurve(IGES_NURBSCurve &c) { _curves.push_back(c); }
	void addNURBSSurface(IGES_NURBSSurface &s) { _surface.push_back(s); }
	void addSurface(Surface2D &s) { _surface2d.push_back(s); }
	void addShapeGroup(ShapeGroup &sg) { _group.push_back(sg); }

	void getLine(const int idx, Line &l) const { l = _lines[idx]; };
	void getArc(const int idx, CircularArc &l) const { l = _arcs[idx]; };
	void getNURBSCurve(const int idx, IGES_NURBSCurve &l) const { l = _curves[idx]; };
	void getNURBSSurface(const int idx, IGES_NURBSSurface &l) const { l = _surface[idx]; };
	void getSurface(const int idx, Surface2D &s) const { s = _surface2d[idx]; };
	void getShapeGroup(const int idx, ShapeGroup &sg) const { sg = _group[idx]; }

	int getLineSize() const { return _lines.size(); };
	int getArcSize() const { return _arcs.size(); };
	int getNURBSCurveSize() const { return _curves.size(); };
	int getNURBSSurfaceSize() const { return _surface.size(); };
	int getSurfaceSize() const { return _surface2d.size(); };
	int getShapeGroupSize() const { return _group.size(); }
};


/***************************************** Shape-Manger2DSurfs *******************************************/
class ShapeManager2DSurfs {
protected:
	std::vector<Line> _lines;									// type_id = 0
	std::vector<CircularArc> _arcs;								// type_id = 1
	std::vector<IGES_NURBSCurve> _curves;						// type_id = 2
	struct Member {
	public :
		int _idx;
		int _typ;
		int _dir = 1;
		Member(const int idx, const int typ, const int dir) { _idx = idx; _typ = typ; _dir = dir; }
	};
	std::vector<std::vector<Member>> _surfs;
	std::vector<int> _level;
	std::vector<std::string> _layerName;

public:
	size_t getSurfSize() const { return _surfs.size(); }
	size_t getSurfSize(const int i) const { return _surfs[i].size(); }
	int getShapeType(const int i, const int j) const { return _surfs[i][j]._typ; }
	int getShapeDir(const int i, const int j) const { return _surfs[i][j]._dir; }
	void getLines(std::vector<Line> &lines) const { lines = _lines; }
	void getArcs(std::vector<CircularArc> &arcs) const { arcs = _arcs; }
	void getNurbsCurves(std::vector<IGES_NURBSCurve> &nurbs) const { nurbs = _curves; }
	void set_layer(const std::string &l, const int i) { _layerName[i] = l; }
	void set_level(const int l, const int i) { _level[i] = l; }
	std::string get_layer(const int i) const { return _layerName[i]; };
	int get_level(const int i) const { return _level[i]; }
	int getDir(const int i, const int j) const { return _surfs[i][j]._dir; }

	void addSurf() { _surfs.push_back(std::vector<Member>()); _level.push_back(-1); _layerName.push_back(""); };
	void addLine(const int sid, const Line &l,const int dir);
	void addArc(const int sid, CircularArc &a, const int dir);
	void addNURBSCurve(const int sid, const IGES_NURBSCurve &c, const int dir);

	void discrete(const int i, const int j,const double TOL, std::vector<NS::Point3D> &pts) ;
	void discrete(const int i ,const double TOL, std::vector<std::vector<NS::Point3D>> &pts) ;

	void removeConincideLines();
	void removeConincideArcs();
	void removeConcincideNURBSCrrves();

};
NS_END