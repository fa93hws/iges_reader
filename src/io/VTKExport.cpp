#include "VTKExport.h"

namespace IGES {
	int VTKExport::findDulp(NS::Point3D &p, std::vector<NS::Point3D> &pts) {
		for (size_t  i = 0; i < pts.size(); i++) {
			if (p == pts[i])
				return i;
		}
		return -1;
	}

	void VTKExport::writeLines(const std::vector<std::vector<int>> &line, std::ofstream &myfile) const {
		if (line.size() == 0) return;
		int lineSize = 0;
		for (std::vector<int> t : line) lineSize += t.size();
		myfile << "LINES\t\t" << line.size() << "\t\t" << lineSize + line.size() << std::endl;
		for (std::vector<int> t : line) {
			myfile << "\t\t" << t.size();
			for (int tt : t)
				myfile << "\t\t" << tt;
			myfile << std::endl;
		}
		myfile << std::endl;
	}
	void VTKExport::writeHeader(std::ofstream &myfile) const {
		myfile << "# vtk DataFile Version 3.1" << std::endl;
		myfile << "Something" << std::endl;
		myfile << "ASCII" << std::endl;
		myfile << "DATASET POLYDATA" << std::endl;
	}
	void VTKExport::writeCell(const std::vector<std::vector<int>> &poly, std::ofstream &myfile) const {
		if (poly.size() == 0)	return;
		int cellSize = 0;
		for (std::vector<int> t : poly) cellSize += t.size();
		myfile << "POLYGONS\t\t" << poly.size() << "\t\t" << cellSize + poly.size() << std::endl;
		for (std::vector<int> t : poly) {
			myfile << "\t\t" << t.size();
			for (int tt : t)
				myfile << "\t\t" << tt;
			myfile << std::endl;
		}
		myfile << std::endl;

		//myfile << "CELL_TYPE\t\t" << poly.size() << endl << "\t";
		//for (int i = 0; i < poly.size(); i++)
		//	myfile << "5\t\t";
	}
	void VTKExport::writePointsCoord(const std::vector<NS::Point3D> &pts, std::ofstream &myfile) const {
		myfile << "POINTS\t" << pts.size() << "\tFLOAT" << std::endl;
		for (NS::Point3D pt : pts) {
			myfile << std::fixed << std::setprecision(5);
			pt.x() > -1e-6 ? myfile << "\t " << pt.x() : myfile << "\t" << pt.x();
			pt.y() > -1e-6 ? myfile << "\t\t " << pt.y() : myfile << "\t\t" << pt.y();
			pt.z() > -1e-6 ? myfile << "\t\t " << pt.z() : myfile << "\t\t" << pt.z();
			myfile << std::endl;
		}
		myfile << std::endl;
	}
	void VTKExport::writeCellData(const std::vector<int> &idx, std::ofstream &myfile) const {
		if (idx.empty())
			return;
		myfile << std::endl;
		myfile << "CELL_DATA\t" << idx[idx.size() - 1] << std::endl;
		myfile << "SCALARS\tLine\tint" << std::endl;
		myfile << "LOOKUP_TABLE default" << std::endl;
		int j = 0;
		for (int i = 0; i < idx[idx.size() - 1]; i++) {
			if (i >= idx[j])
				j++;
			myfile << "\t" << j << std::endl;
		}
		myfile << std::endl;
	}
	void VTKExport::generateVTKFile() {
		std::ofstream myfile;
		std::string path = outPutPath + exampleName + "\\";
		CreateDirectory(path.c_str(), NULL);
		myfile.open(path + "line.vtk");
		writeHeader(myfile);
		writePointsCoord(curvePlotInfo.pts, myfile);
		writeLines(curvePlotInfo.connev, myfile);
		writeCellData(curvePlotInfo.idx, myfile);
		myfile.close();
		//
		myfile.open(path + "line_control.vtk");
		writeHeader(myfile);
		writePointsCoord(curvePlotInfo.controlPoints, myfile);
		writeLines(curvePlotInfo.cpConnev, myfile);
		writeCellData(curvePlotInfo.cpIdx, myfile);
		myfile.close();
		//
		myfile.open(path + "surface.vtk");
		writeHeader(myfile);
		writePointsCoord(surfacePlotInfo.pts, myfile);
		writeCell(surfacePlotInfo.connev, myfile);
		writeCellData(surfacePlotInfo.idx, myfile);
		myfile.close();
		//
		myfile.open(path + "surface_control.vtk");
		writeHeader(myfile);
		writePointsCoord(surfacePlotInfo.controlPoints, myfile);
		writeCell(surfacePlotInfo.cpConnev, myfile);
		writeCellData(surfacePlotInfo.cpIdx, myfile);
		myfile.close();
		//
		myfile.open(path + "surface_control_convexHull.vtk");
		writeHeader(myfile);
		writePointsCoord(surfacePlotInfo.convexHullPoints, myfile);
		writeCell(surfacePlotInfo.chpConnec, myfile);
		writeCellData(surfacePlotInfo.cvIdx, myfile);
		myfile.close();
	}

	void VTKExport::addCurve(NURBS::Curve &c) {
		//Plot Points
		NURBS::Curve::CurveInfo curveInfo;
		c.getScattors(200, curveInfo);
		std::vector<std::vector<NS::Point3D>> scatter = curveInfo.pnts;
		std::vector<NS::Point3D> pts;
		std::vector<std::vector<int>> controlNet;
		std::vector<std::vector<int>> lines;
		int index = 0;
		lines.push_back({});
		for (size_t j = 0; j < scatter[0].size(); j++) {
			pts.push_back(scatter[0][j]);
			lines[0].push_back(index);
			index++;
		}
		//Control points
		controlNet.push_back({});
		index = 0;
		for (NS::Point3D p : c.controlPoints) {
			controlNet[0].push_back(index);
			index++;
		}
		// finish
		curvePlotInfo.contGeo(pts, lines);
		curvePlotInfo.contControlPoints(c.controlPoints, controlNet);
	}

	void VTKExport::addSurface(NURBS::Surface &s) {
		//Plot Points
		NURBS::Surface::SurfaceInfo surfaceInfo;
		s.getScattors(20, surfaceInfo);
		std::vector<std::vector<NS::Point3D>> scatter = surfaceInfo.pnts;
		std::vector<NS::Point3D> pts;
		std::vector<std::vector<int>> poly;
		int nRow = scatter.size();
		int nCol = scatter[0].size();
		for (std::vector<NS::Point3D> P3 : scatter)
			for (NS::Point3D P33 : P3)
				pts.push_back(P33);
		for (int i = 0; i < nRow - 1; i++)
			for (int j = 0; j < nCol - 1; j++) {
				int _poly[] = { i * nCol + j, i * nCol + j + 1 , (i + 1) * nCol + j + 1 , (i + 1) * nCol + j };
				std::vector<int> v_poly(_poly, _poly + 4);
				poly.push_back(v_poly);
			}
		// Control Polygons
		std::vector<std::vector<NS::Point3D>> cp = s.controlPoints;
		nRow = cp.size();
		nCol = cp[0].size();
		std::vector<std::vector<int>> cpConnev;
		for (int i = 0; i < nRow - 1; i++)
			for (int j = 0; j < nCol - 1; j++) {
				int _poly[] = { i * nCol + j, i * nCol + j + 1 , (i + 1) * nCol + j + 1 , (i + 1) * nCol + j };
				std::vector<int> v_poly(_poly, _poly + 4);
				cpConnev.push_back(v_poly);
			}
		std::vector<NS::Point3D> _cps; // _cps: convert 2d matrix into 1d for plot
		for (std::vector<NS::Point3D> vp : s.controlPoints)
			for (NS::Point3D p : vp)
				_cps.push_back(p);
		// ConvexHull Polygons
		if (!s.convexHull.empty()) {
			size_t nFace = s.convexHull.size();
			std::vector<NS::Point3D> _cvs;
			std::vector<std::vector<int>> cvConnev;
			for (size_t i = 0; i < s.controlPoints.size(); i++)
				for (size_t j = 0; j < s.controlPoints[i].size(); j++)
					_cvs.push_back(s.controlPoints[i][j]);
			for (size_t i = 0; i < nFace; i++)
				cvConnev.push_back({ s.convexHull[i].I[0],  s.convexHull[i].I[1] , s.convexHull[i].I[2] });
			surfacePlotInfo.contConvexHull(_cvs, cvConnev);
		}
		//vector<vector<Point3D>> cv = s.convexHull;
		//if (!cv.empty()) {
		//	nRow = cv.size();
		//	nCol = cv[0].size();
		//	vector<vector<int>> cvConnev;
		//	for (int i = 0; i < nRow - 1; i++)
		//		for (int j = 0; j < nCol - 1; j++) {
		//			int _poly[] = { i * nCol + j, i * nCol + j + 1 , (i + 1) * nCol + j + 1 , (i + 1) * nCol + j };
		//			vector<int> v_poly(_poly, _poly + 4);
		//			cvConnev.push_back(v_poly);
		//		}
		//	vector<Point3D> _cvs; // _cps: convert 2d matrix into 1d for plot
		//	for (vector<Point3D> vp : s.convexHull)
		//		for (Point3D p : vp)
		//			_cvs.push_back(p);
		//	surfacePlotInfo.contConvexHull(_cvs, cvConnev);
		//}
		// finish
		surfacePlotInfo.contGeo(pts, poly);
		surfacePlotInfo.contControlPoints(_cps, cpConnev);

	}
}