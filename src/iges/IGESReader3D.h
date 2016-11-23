#pragma once
#include <string>
#include <iostream>
#include <stdlib.h> 
#include <fstream>
#include <math.h>
#include <vector>
#include <string>
#include "../nurbs/curve.h"
#include "../nurbs/surface.h"
#include "Shape3D.h"
#include "../shared/utility.h"
#include <armadillo>
#include <memory> //std::shared_ptr
#define _PARSEINPUT const std::vector<std::vector<std::string>> &dt_all, const int dt_idx, const std::vector<std::string> &pt


NSI_BEG
class IGESReader3D {
private:
	SolidsGroup _solidsGroup;
	SurfacesGroup _surfaceGroup;
	CurvesGroup _curvesGroup;
	Vertex502 _vertex;

	std::string _empty = "        ";
	struct StatusNumber {
	public:
		int blank;
		int subordinate;
		int useFlag;
		int hierarchy;
		StatusNumber(const std::string &status_number) {
			blank = stoi(status_number.substr(0, 2));
			subordinate = stoi(status_number.substr(2, 2));
			useFlag = stoi(status_number.substr(4, 2));
			hierarchy = stoi(status_number.substr(6, 2));
		};
	};
	double _stod(const std::string &str);
	
	// D,P part
	void readPPart(const std::string &line, std::vector<std::string> &p_text) const;
	void getPPart(const int pid_begin, const int pid_end, const std::vector<std::string> &p_text, std::vector<std::string> &pt) const;
	void readDPart(const std::string &line, std::vector<std::vector<std::string>> &d_text) const;
	// Shapes
	void parseAllShapes(std::vector<std::vector<std::string>> &dt, const std::vector<std::string> &pt, SolidsGroup &sg);
	void parseSolid186(_PARSEINPUT, Solid186 &solid);
	void parseShell514(_PARSEINPUT, Shell514 &shell);
	void parseFace510(_PARSEINPUT, Face510 &face);
	void parseLoop508(_PARSEINPUT, Loop508 &loop);
	int parseSurf(_PARSEINPUT, SurfacesGroup &sg);
	void parseNURBSSurf128(_PARSEINPUT, Surf128 &ss); // nurbs surf
	int parseCurve(_PARSEINPUT,CurvesGroup &cg);
	void parseCurve126(_PARSEINPUT, Curve126 &cc); // nurbs curve
	int parseVer502(_PARSEINPUT, const int &ptidx, Vertex502 &ver);
	void parseEdge504(_PARSEINPUT,const int edgeIdx, Edge504 &edge);
	void parseMatrixTran(_PARSEINPUT, std::shared_ptr<Shape3D> &s3d);

public:
	IGESReader3D() {};
	IGESReader3D(std::string &path);
	void getCurvesGroup(CurvesGroup &cg) const { cg = _curvesGroup; };
};
NS_END