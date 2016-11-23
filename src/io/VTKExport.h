#pragma once
#include <string>
#include "../geo/point3d.h"
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <Windows.h>
#include <stdio.h>
#include "../nurbs/curve.h"
#include "../nurbs/surface.h"
#include "../IGES/shape.h"
#include "../IGES/IGESReader.h"

namespace IGES {
	class VTKExport {
	private:
		std::string outPutPath;
		std::string exampleName;
		struct PlotInfo {
		public:
			std::vector<int> idx;
			std::vector<int> cpIdx;
			std::vector<int> cvIdx;

			std::vector<NS::Point3D> pts;
			std::vector<std::vector<int>> connev;
			std::vector<NS::Point3D> controlPoints;
			std::vector<std::vector<int>> cpConnev;
			std::vector<NS::Point3D> convexHullPoints;
			std::vector<std::vector<int>> chpConnec;
			PlotInfo() { };
			void contGeo(std::vector<NS::Point3D> &points, std::vector<std::vector<int>> &curveConnev) {
				int oriNumPts = pts.size();
				for (NS::Point3D p : points)
					pts.push_back(p);
				if (curveConnev.size() > 0)
					for (size_t  i = 0; i < curveConnev.size(); i++)
						for (size_t  j = 0; j < curveConnev[i].size(); j++)
							curveConnev[i][j] += oriNumPts;
				for (std::vector<int> cc : curveConnev)
					connev.push_back(cc);
				idx.push_back(connev.size());
			}
			void contControlPoints(std::vector<NS::Point3D> &cp, std::vector<std::vector<int>> &_cpConnev) {
				int oriNumPts = controlPoints.size();
				for (NS::Point3D p : cp)
					controlPoints.push_back(p);
				for (size_t  i = 0; i < _cpConnev.size(); i++)
					for (size_t  j = 0; j < _cpConnev[i].size(); j++)
						_cpConnev[i][j] += oriNumPts;
				for (std::vector<int> cc : _cpConnev)
					cpConnev.push_back(cc);
				cpIdx.push_back(cpConnev.size());
			}
			void contConvexHull(std::vector<NS::Point3D> &cp, std::vector<std::vector<int>> &_cpConnev) {
				int oriNumPts = convexHullPoints.size();
				for (NS::Point3D p : cp)
					convexHullPoints.push_back(p);
				for (size_t  i = 0; i < _cpConnev.size(); i++)
					for (size_t  j = 0; j < _cpConnev[i].size(); j++)
						_cpConnev[i][j] += oriNumPts;
				for (std::vector<int> cc : _cpConnev)
					chpConnec.push_back(cc);
				cvIdx.push_back(chpConnec.size());
			}
		} curvePlotInfo, surfacePlotInfo;
		int findDulp(NS::Point3D &p, std::vector<NS::Point3D> &pts);
		void writeHeader(std::ofstream &myfile) const;
		void writeLines(const std::vector<std::vector<int>> &line, std::ofstream &myfile) const;
		void writeCell(const std::vector<std::vector<int>> &poly, std::ofstream &myfile) const;
		void writePointsCoord(const std::vector<NS::Point3D> &pts, std::ofstream &myfile) const;
		void writeCellData(const std::vector<int> &idx, std::ofstream &myfile) const;
		

	public:
		VTKExport() {};
		VTKExport(std::string pathOut, std::string name) {
			outPutPath = pathOut;
			exampleName = name;
		}
		void setPath(const std::string path) { outPutPath = path; };
		void setName(const std::string name) { exampleName = name; };
		void addCurve(NURBS::Curve &c);
		void addSurface(NURBS::Surface &s);
		void generateVTKFile();
	};

}