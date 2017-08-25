
#include "ReadIGES.h"
#include "JSON.hpp"
#include "../iges/IGESReader.h"
#include "../io/SbfemExport.h"

JNIEXPORT jstring JNICALL Java_ReadIGES_readFile
(JNIEnv * env , jclass jboj, jstring __path) {
	const char *_path = env->GetStringUTFChars(__path, 0);
	std::string path(_path);

	IGES::IGESReader gen = IGES::IGESReader(path);
	IGES::ShapeManager2DSurfs s2d; gen.getsurfaces2D(s2d);
	IGES::SBFEMExport se(s2d, 0.001, 150);
	//pts
	std::vector<NS::Point3D> _pts;
	se.getKPT(_pts);
	nlohmann::json coor;
	for (NS::Point3D pt : _pts) {
		nlohmann::json _pt;
		_pt["x"] = pt.x();
		_pt["y"] = pt.y();
		coor.push_back(_pt);
	}
	nlohmann::json out;
	out["coor"] = coor;
	// polylines
	std::vector<std::vector<int>> polyLines;
	se.getPolyLines(polyLines);
	out["polylines"] = polyLines;
	// surfaces
	std::vector<std::vector<int>> _surfaces;
	std::vector<std::vector<bool>> _reversed;
	se.getSurfaces(_surfaces);
	se.getReversed(_reversed);
	nlohmann::json surfaces;
	for (size_t i = 0; i < _surfaces.size(); i++) {
		nlohmann::json temp;
		for (size_t j = 0; j < _surfaces[i].size(); j++) {
			temp["surf"].push_back(_surfaces[i][j]);
			temp["dir"].push_back(_reversed[i][j] ? -1 : 1);
		}
		surfaces.push_back(temp);
	}
	out["surfaces"] = surfaces;
	// isHardpt
	std::vector<bool> isHardpt;
	se.getIsHardPt(isHardpt);
	out["isHardPt"] = isHardpt;
	// level
	std::vector<int> level;
	se.getLevel(level);
	out["level"] = level;
	// layerName
	std::vector<std::string> layerNames;
	se.getLayerName(layerNames);
	out["layerNames"] = layerNames;

	const char* outString = out.dump().c_str();
	jstring outJString = env->NewStringUTF(outString);
	return outJString;
}