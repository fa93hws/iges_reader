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
#include "Shape.h"
#include "../shared/utility.h"

namespace IGES {
	class IGESReader {
	private:
		std::vector<int> entities_typ;
		std::vector<std::string> entities_name;
		std::vector<std::vector<std::string>> s_text;
		std::vector<std::vector<std::string>> g_text;
		std::vector<std::vector<std::string>> d_text;
		std::vector<std::string>			  p_text;
		std::vector<std::vector<std::string>> t_text;
		std::vector<std::string> layer_name;
		ShapeManger  shape_collection;
		ShapeManager2DSurfs _surfs2D;
		char subClass[5] = { 'S','G','D','P','T' };
		std::ifstream myfile;
		struct StatusNumber {
		public:
			int blank;
			int subordinate;
			int useFlag;
			int hierarchy;
			StatusNumber(const std::string &status_number){
				blank = stoi(status_number.substr(0, 2));
				subordinate = stoi(status_number.substr(2, 2));
				useFlag = stoi(status_number.substr(4, 2));
				hierarchy = stoi(status_number.substr(6, 2));
			};
		};

		double _stod(const std::string &str);
		// Subrutine

		template <typename T>
		std::vector<std::vector<T>> transpose(std::vector<std::vector<T>> in) {
			int nRow = in.size();
			int nCol = in[0].size();
			std::vector<std::vector<T>> out(nCol);
			for (int i = 0; i < nCol; i++)
				for (int j = 0; j < nRow; j++)
					out[i].push_back(in[j][i]);
			return out;
		}

		//auto parseEntityByDT(const std::vector<std::string> &dt);
		//Method
		void readSPart(const std::string &line);
		void savePPart(const std::string &line);
		void readPPart(const int pid_begin, const int pid_end, std::vector<std::string> &pt);
		void readDPart(const std::string &line);
		void readGPart(const std::string &line);
		void readTPart(const std::string &line);
		void reOrganizePart(const std::string &line);

		//parse entities


		void parseNurbSurface(const std::vector<std::string> &pt, NURBS::ParentSurface &ss);
		
		void parse_matrix_tran(const std::vector<std::string> &dt, std::vector<std::vector<double>> &R, std::vector<double> &T);
		void set_attribute(const std::vector<std::string> &dt,Shape &s);
		void set_attribute(const std::vector<std::string> &dt, ShapeManager2DSurfs &s, const int surf_id);
		//void set_shape_attributes(const std::string &shapeStringIdx, const int parentShapeIdx,StatusNumber &status, const StatusNumber &parentStatus, Shape &shape);
		/*
			h_flag = 0 ====> use parent's attributes
			h_flag = 1 ====> use own attributes
			h_flag = 2 ====> use parents' attributes when own attributes is empty.
		*/

		void parseCircularArc(const std::vector<std::string> &dt, CircularArc &cc);																				// 100
		void parseLine(const std::vector<std::string> &dt, Line &line);																							// 110		
		void parseNurbCurve(const std::vector<std::string> &dt,  IGES_NURBSCurve &curve);																		// 126
		void parseSurface2D143(const std::vector<std::string> &dt, ShapeManager2DSurfs &surface2D);																			// 143
		void parseSubFigure408(const std::vector<std::string> &dt, ShapeGroup &sg);																				// 408
		void parseProp406_layer(const std::vector<std::string> &dt, std::vector<std::string> &layer_name);														// 406


	public:		
		//Cons
		IGESReader() {}
		IGESReader(std::string &path);

		//io
		//SubSection getSubSection() { return p; };

		//method		
		void parseAllShapes();
		void parseEntities(std::vector<NURBS::ParentSurface> &s3D, std::vector<NURBS::Curve> &c2D, std::vector<IGES::Shape> &shp);
		void parseEntities(std::vector<NURBS::ParentSurface> &s3D, std::vector<NURBS::Curve> &c2D) { std::vector<IGES::Shape> _t; parseEntities(s3D, c2D, _t); }
		void getShapes(ShapeManger &s) const { s= shape_collection; }
		void getsurfaces2D(ShapeManager2DSurfs &s2d) const {s2d = _surfs2D; }

		
	};
}