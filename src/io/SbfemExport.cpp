#include "SbfemExport.h"
#include "iostream"
#include <math.h>
#include <Windows.h>
#include <iomanip>      // std::setprecision
#include "../geo/geometry.h"



namespace IGES {
	// IO
	void SBFEMExport::getSurfaces(std::vector<std::vector<int>> &surf) const {
		surf.clear();
		for (auto i = _surfaces.begin(); i != _surfaces.end(); i++) {
			std::vector<int> tmp;
			for (auto j = (*i).begin(); j != (*i).end(); j++) {
				tmp.push_back((*j)._idx);
			}
			surf.push_back(tmp);
		}
	}
	void SBFEMExport::getReversed(std::vector<std::vector<bool>> &rev) const {
		rev.clear();
		for (auto i = _surfaces.begin(); i != _surfaces.end(); i++) {
			std::vector<bool> tmp;
			for (auto j = (*i).begin(); j != (*i).end(); j++) {
				tmp.push_back((*j)._reversed);
			}
			rev.push_back(tmp);
		}
	}
	//
	double SBFEMExport::calculateAngle3pt(const NS::Point3D &pt0, const NS::Point3D &pt1, const NS::Point3D &pt2) {
		double l01 = pt0.distanceTo(pt1);
		double l02 = pt0.distanceTo(pt2);
		double l12 = pt1.distanceTo(pt2);
		assert(l01 > 1e-10 && l02 > 1e-10 && l12 > 1e-10);
		double para = l01*l01 + l12*l12 - l02 * l02;
		para = para / 2.0 / l01 / l12;
		return acos( para) * 180.0 / 3.141592654;
	}

	int SBFEMExport::buildPolyLines(const std::vector<NS::Point3D> &member_pts, int pt_idx) {
		std::vector<int> connev;
		for (NS::Point3D pt : member_pts) {
			_keyPoints.push_back(pt);
			connev.push_back(pt_idx);
			pt_idx++;
		}
		_polyLines.push_back(connev);
		return connev.size();
	}
	//void SBFEMExport::parseShapeGrouptoPolylines(const ShapeGroup &sg, const double TOL) {
	//	int pt_idx = _keyPoints.size();
	//	for (int i = 0; i < sg.getNURBSCurveSize(); i++) {
	//		IGES_NURBSCurve curve;
	//		sg.getNURBSCurve(i, curve);
	//		std::vector<NS::Point3D> member_pts;
	//		curve.discrete(TOL, member_pts);
	//		pt_idx += buildPolyLines(member_pts, pt_idx);
	//		std::vector<double> c;
	//		if (sg.getpropHierarchy()) {
	//			sg.get_color(c);
	//			_color.push_back(c);
	//			_layerName.push_back(sg.get_layer());
	//		}
	//		else {
	//			curve.get_color(c);
	//			_color.push_back(c);
	//			_layerName.push_back(curve.get_layer());
	//		}
	//	}
	//	for (int i = 0; i < sg.getArcSize(); i++) {
	//		CircularArc arc;
	//		sg.getArc(i, arc);
	//		std::vector<NS::Point3D> member_pts;
	//		arc.discrete(TOL, member_pts);
	//		pt_idx += buildPolyLines(member_pts, pt_idx);
	//		std::vector<double> c;
	//		if (sg.getpropHierarchy()) {
	//			sg.get_color(c);
	//			_color.push_back(c);
	//			_layerName.push_back(sg.get_layer());
	//		}
	//		else {
	//			arc.get_color(c);
	//			_color.push_back(c);
	//			_layerName.push_back(arc.get_layer());
	//		}
	//	}
	//	for (int i = 0; i < sg.getLineSize(); i++) {
	//		Line line;
	//		sg.getLine(i, line);
	//		std::vector<NS::Point3D> member_pts;
	//		line.discrete(TOL, member_pts);
	//		pt_idx += buildPolyLines(member_pts, pt_idx);
	//		std::vector<double> c;
	//		if (sg.getpropHierarchy()) {
	//			sg.get_color(c);
	//			_color.push_back(c);
	//			_layerName.push_back(sg.get_layer());
	//		}
	//		else {
	//			line.get_color(c);
	//			_color.push_back(c);
	//			_layerName.push_back(line.get_layer());
	//		}
	//	}
	//}
	void SBFEMExport::findHardPt(const double angle) {
		_isHardpt.clear();
		_isHardpt.resize(_keyPoints.size(), false);
		for (auto i = _surfaces.cbegin(); i != _surfaces.cend();i++) {
			for (auto j = (*i).cbegin(); j != (*i).cend();j++) {
				auto line = _polyLines[(*j)._idx];
				bool isThisReversed = (*j)._reversed;
				// check is there any hard points within the lines
				if (_polyLines[(*j)._idx].size() > 2) {					
					for (auto k = line.cbegin()+1; k != line.cend()-1; k++) {
						auto pt0 = _keyPoints[*(k - 1)];
						auto pt1 = _keyPoints[*(k)];
						auto pt2 = _keyPoints[*(k + 1)];
						double ang = calculateAngle3pt(pt0, pt1, pt2);
						if (ang < angle)
							_isHardpt[*(k)] = true;
					}
				}
				// check is the first point in the poline hard point
				if (isThisReversed) std::reverse(line.begin(), line.end());
				std::vector<int> preLine;
				bool isPreReversed;
				if ((*i).size() == 1) { // closed;
					preLine = line;
					isPreReversed = false;
				}
				else if (j == (*i).cbegin()) {
					preLine = _polyLines[(*i).back()._idx];
					isPreReversed = (*i).back()._reversed;
				}
				else {
					preLine = _polyLines[(*(j - 1))._idx];
					isPreReversed = (*(j - 1))._reversed;
				}
				if (isPreReversed) std::reverse(preLine.begin(),preLine.end());
				auto pt0 = _keyPoints[ *(preLine.end() - 2) ];
				auto pt1 = _keyPoints[line[0]];
				auto pt2 = _keyPoints[line[1]];
				if (calculateAngle3pt(pt0, pt1, pt2) < angle)
					_isHardpt[line[0]] = true;
			}
		}
	}
	void SBFEMExport::parseSurf2D(ShapeManager2DSurfs &s2d, const double TOL) {
		_surfaces = std::vector<std::vector<Member>>(s2d.getSurfSize(), std::vector<Member>());
		// find all points
		int pt_idx = 0;
		for (size_t i = 0; i < s2d.getSurfSize(); i++) {
			std::vector<std::vector<NS::Point3D>> member_pts;
			s2d.discrete(i,TOL, member_pts);
			std::vector<int> surfconnev;
			std::vector<bool> surfreverse;
			for (auto j = member_pts.begin(); j != member_pts.end(); j++) {
				pt_idx += buildPolyLines(*j, pt_idx);
				int dir = s2d.getDir(i, j - member_pts.begin());
				_surfaces[i].push_back(Member(_polyLines.size() - 1, dir == -1));
			}
			_level.push_back(s2d.get_level(i));
			_layerName.push_back(s2d.get_layer(i));			
		}

		double xmax = _keyPoints[0].x(), xmin = _keyPoints[0].x();
		double ymax = _keyPoints[0].y(), ymin = _keyPoints[0].y();
		double zmax = _keyPoints[0].y(), zmin = _keyPoints[0].z();
		for (NS::Point3D pt : _keyPoints) {
			if (pt.x() > xmax) xmax = pt.x();
			if (pt.x() < xmin) xmin = pt.x();
			if (pt.y() > ymax) ymax = pt.y();
			if (pt.y() < ymin) ymin = pt.y();
			if (pt.z() > zmax) zmax = pt.z();
			if (pt.z() < zmin) zmin = pt.z();
		}
		double diffmax = NS::Maximum3(xmax - xmin, ymax - ymin, zmax - zmin);
		double scale = diffmax / 20.0;
		NS::Point3D offset;
		offset.set((xmax + xmin) / 2.0, (ymax + ymin) / 2.0, (zmax + zmin) / 2.0);
		std::vector<NS::Point3D> uniformalized_kpt;
		for (NS::Point3D pt : _keyPoints)
			uniformalized_kpt.push_back((pt - offset) / scale);
		//  rem coincide pts
		std::vector<int> multi_idx(_keyPoints.size(), -1);
		for (size_t  i = 0; i < uniformalized_kpt.size() - 1; i++) {
			if (multi_idx[i]>-1) continue;
			for (size_t  j = i + 1; j < uniformalized_kpt.size(); j++)
				if (uniformalized_kpt[i].distanceTo(uniformalized_kpt[j]) < 1e-8)
					multi_idx[j] = i;
		}
		std::vector<NS::Point3D> newKpt;
		for (size_t  i = 0; i < _keyPoints.size(); i++)
			if (multi_idx[i] < 0)
				newKpt.push_back(_keyPoints[i]);
		_keyPoints = newKpt;
		// rebuild connev
		int shift = 0;
		for (size_t  i = 0; i < multi_idx.size();i++) {
			if (multi_idx[i] > -1) shift++;
			if (shift == 0) continue;
			if (multi_idx[i] == -1)
				multi_idx[i] = i - shift;
			else
				multi_idx[i] = multi_idx[multi_idx[i]] == -1 ? multi_idx[i] : multi_idx[multi_idx[i]];
		}
		for (std::vector<int> &con : _polyLines)
			for (int &idx : con)
				if (multi_idx[idx]>-1)
					idx = multi_idx[idx];

		// remove coincidence in _polylines
		std::vector<int> lineMultiIdx = std::vector<int>(_polyLines.size(), -1);
		std::vector<int> dir = std::vector<int>(_polyLines.size(), 1);
		for (size_t i = 0; i < _polyLines.size()-1; i++) {
			if (lineMultiIdx[i] > -1) continue;
			for (size_t j = i + 1; j < _polyLines.size(); j++) {
				if (lineMultiIdx[i] > -1) continue;
				if (_polyLines[i] == _polyLines[j]) {
					lineMultiIdx[j] = i;
					continue;
				}
				else if (_polyLines[i].back() != _polyLines[j].front()) continue;
				else {
					std::vector<int> temp = _polyLines[j];
					std::reverse(temp.begin(), temp.end());
					if (temp == _polyLines[i]) {
						lineMultiIdx[j] = i;
						dir[j] = -1;
					}
				}
			}
		}
		int count = 0;
		for (size_t i = 0; i < lineMultiIdx.size(); i ++) {
			if (lineMultiIdx[i] == -1) {
				lineMultiIdx[i] = i - count;
				_polyLines[i - count] = _polyLines[i];
			}
			else {
				lineMultiIdx[i] = lineMultiIdx[lineMultiIdx[i]];
				count++;
			}
		}
		for (int i = 0; i < count; i++) _polyLines.pop_back();
		for (std::vector<Member> &m : _surfaces) {
			for (Member &mm : m) {
				if (dir[mm._idx] == -1) mm._reversed = !mm._reversed;
				mm._idx = lineMultiIdx[mm._idx];				
			}
		}
	}
	void SBFEMExport::filterByLayer(const std::vector<std::string> &layer, const char mode) {
		bool include;
		if (mode == 'i' || mode == 'I')
			include = true;
		else if (mode == 'o' || mode == 'O')
			include = false;
		else
			ASSERT(0);
		std::vector<bool> inRecord(_polyLines.size(),false);
		for (size_t  i = 0; i < _polyLines.size(); i++) {
			for (size_t  j = 0; j < layer.size(); j++) {
				if (_layerName[i].compare(layer[j]) == 0) {
					inRecord[i] = true;
					break;
				}
			}
		}
		std::vector<std::vector<int>> polyLines;
		for (size_t  i=0; i < _polyLines.size(); i++)
			if ((include && inRecord[i]) || (!include && !inRecord[i]))
				polyLines.push_back(_polyLines[i]);
		_polyLines = polyLines;
	}
	void SBFEMExport::setLayerColor(const std::string &layer, const int color) {
		auto it = _layerName.cbegin();
		do {
			it = std::find(it, _layerName.cend(), layer);			
			if (it == _layerName.cend()) break;
			_level[it++ - _layerName.cbegin()] = color;
		} while (true);
	}
	void SBFEMExport::setLayerColor(const std::vector<int> &level, const std::vector<int> &color) {
		std::vector<std::vector<int>> existIdx;
		for (size_t i = 0; i < level.size(); i++) {
			existIdx.push_back(std::vector<int>());
			for (size_t j = 0; j < _level.size(); j++) {
				if (level[i] == _level[j]) {
					existIdx[i].push_back(j);
				}
			}
		}
		for (size_t i = 0; i < existIdx.size(); i++) {
			for (size_t j = 0; j < existIdx[i].size(); j++) {
				_level[existIdx[i][j]] = color[i];
			}
		}
	}
	SBFEMExport::SBFEMExport(ShapeManager2DSurfs &in,const double TOL, const double angle) {
		parseSurf2D(in,TOL);
		findHardPt(angle);
	}

	void SBFEMExport::export2MFile2D(const std::string &pathOut, const std::string &name) const {
		std::ofstream myfile;
		std::string path = pathOut + "\\";
		CreateDirectory(path.c_str(), NULL);
		myfile.open(path + name + "_input.m");
		myfile << "% Auto Generated" << std::endl;
		myfile << "function[objnodes, objlines, objsurfs, objregions, regioncolors,ishardPt] = ";
		myfile << name << "_input()" << std::endl << std::endl;
		// objnodes
		myfile << "objnodes = [" << std::endl;
		int i = 0;
		for (NS::Point3D pt : _keyPoints) {
			myfile << std::fixed << std::setprecision(5);
			pt.x() > -1e-6 ? myfile << "\t " << pt.x() : myfile << "\t" << pt.x();
			myfile << "\t\t,\t";
			pt.y() > -1e-6 ? myfile << "\t\t " << pt.y() : myfile << "\t\t" << pt.y();
			myfile << ";\t\t%\tpt"<< ++i << std::endl;
		}
		myfile << "]; % " <<i << "pts in total" << std::endl << std::endl;
		// isHardPt
		myfile << "ishardPt = [";
		for (bool b : _isHardpt)
			myfile << b << " ";
		myfile << "];" << std::endl << "ishardPt = ishardPt==1;" << std::endl;
		// objlines
		i = 0;
		myfile << "objlines = cell(" << _polyLines.size() <<",1) "  ";" << std::endl;
		for (std::vector<int> p : _polyLines) {
			myfile << "objlines{" << ++i << "} \t= \t[";
			for (auto pp = p.cbegin(); pp != p.cend() - 1; pp++)
				myfile << *pp + 1<< ",";
			myfile << p.back() + 1<< "];"<<std::endl;
		}
		myfile << std::endl;
		// objsurfs
		i = 0;
		myfile << "objsurfs = cell(" << _surfaces.size() << ",1) " ";" << std::endl;
		for (size_t j = 0; j < _surfaces.size(); j++) {
			myfile << "objsurfs{" << ++i << "} \t = \t[\t";
			for (auto pp = _surfaces[j].cbegin(); pp != _surfaces[j].cend() - 1; pp++)
				myfile << (*pp)._idx + 1 << ",";
			myfile << _surfaces[j].back()._idx + 1 << ";" << std::endl << "\t\t\t\t\t\t";
			for (auto pp = _surfaces[j].cbegin(); pp != _surfaces[j].cend() - 1; pp++)
				myfile << 1 - (*pp)._reversed * 2 <<"," ;
			myfile << (1- _surfaces[j].back()._reversed *2) << "]';" << std::endl;
		}
		// objregions
		i = 0;
		myfile << std::endl << "objregions = cell(" << _surfaces.size() << ",1) " ";" << std::endl;
		for (size_t j = 0; j < _surfaces.size(); j++) {
			myfile << "objregions{" << ++i << "} \t = \t[" << j + 1 << "];" << std::endl;
		}
		// color
		myfile << std::endl << "regioncolors=[";
		for (auto j = _level.cbegin(); j != _level.cend() - 1; j++)
			myfile << (*j) + 1 << ",";
		myfile << _level.back() + 1 << "];" << std::endl << std::endl;
		// layer		
		std::vector<int> layer_table;
		for (int j : _level)
			if (std::find(layer_table.begin(), layer_table.end(), j) == layer_table.end())
				layer_table.push_back(j);
		std::sort(layer_table.begin(), layer_table.end(), [](int a, int b) {return a<b; });

		i = 0;
		myfile << "%%\t\t\t regioncolors - layers" << std::endl;
		myfile << "%color\tlayer_name"<<std::endl;
		for (int j: layer_table)
			myfile << "% " << j+1 << "\t\t" << _layerName[std::find(_level.begin(), _level.end(), j) - _level.begin()] << std::endl;
		// finish
		myfile << "end" << std::endl;
		myfile.close();
		std::cout << "exported to m file" << std::endl;

	}
	void SBFEMExport::export2vtk(const std::string &pathOut, const std::string &name) const  {
		std::ofstream myfile;
		std::string path = pathOut + "\\";
		CreateDirectory(path.c_str(), NULL);
		myfile.open(path + name + ".vtk");
		myfile << "# vtk DataFile Version 3.1" << std::endl;
		myfile << "Something" << std::endl;
		myfile << "ASCII" << std::endl;
		myfile << "DATASET POLYDATA" << std::endl;
		std::vector<std::vector<int>> _poly;
		std::vector<int> _lineValue;
		std::vector<int> _lineId;
		for (size_t i = 0; i < _surfaces.size(); i++) {
			for (size_t j = 0; j < _surfaces[i].size(); j++) {
				_lineId.push_back(_surfaces[i][j]._idx);
				if (_surfaces[i][j]._reversed) {
					std::vector<int> temp;
					temp = _polyLines[_surfaces[i][j]._idx];
					std::reverse(temp.begin(), temp.end());
					_poly.push_back(temp);
				}
				else
					_poly.push_back(_polyLines[_surfaces[i][j]._idx]);
				_lineValue.push_back(_level[i]);
			}
		}
		myfile << "POINTS\t" << _keyPoints.size() << "\tFLOAT" << std::endl;
		for (NS::Point3D pt : _keyPoints) {
			myfile << std::fixed << std::setprecision(5);
			pt.x() > -1e-6 ? myfile << "\t " << pt.x() : myfile << "\t" << pt.x();
			pt.y() > -1e-6 ? myfile << "\t\t " << pt.y() : myfile << "\t\t" << pt.y();
			pt.z() > -1e-6 ? myfile << "\t\t " << pt.z() : myfile << "\t\t" << pt.z();
			myfile << std::endl;
		}
		myfile << std::endl;
		if (_poly.size() == 0) return;
		int lineSize = 0;
		for (std::vector<int> t : _poly) lineSize += t.size();
		myfile << "LINES\t\t" << _poly.size() << "\t\t" << lineSize + _poly.size() << std::endl;
		for (std::vector<int> t : _poly) {
			myfile << "\t\t" << t.size();
			for (int tt : t)
				myfile << "\t\t" << tt;
			myfile << std::endl;
		}
		myfile << std::endl;
		myfile << std::endl;
		myfile << "CELL_DATA\t" << _lineValue.size() << std::endl;
		myfile << "SCALARS\tLayer_Id\tint" << std::endl;
		myfile << "LOOKUP_TABLE default" << std::endl;
		for (int j : _lineValue)
			myfile << "\t\t" << j << std::endl;

		myfile << "SCALARS\tLine_ID\tint" << std::endl;
		myfile << "LOOKUP_TABLE default" << std::endl;
		for (auto j : _lineId)
			myfile << "\t\t" << j << std::endl;
		myfile << std::endl;

		myfile << "POINT_DATA\t" << _isHardpt.size() << std::endl;
		myfile << "SCALARS\tIsHardPt\tint" << std::endl;
		myfile << "LOOKUP_TABLE default" << std::endl;
		for (bool b : _isHardpt)
			myfile << "\t\t" << (int)b << std::endl;
		myfile << std::endl;
		myfile.close();
		std::cout << "exported to VTK file" << std::endl;
	}
};