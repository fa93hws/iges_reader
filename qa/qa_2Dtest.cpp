#include "qa_2Dtest.h"
#include "../src/iges/IGESReader.h"
#include "../src/io/SbfemExport.h"


void QA_IGES2D::dam_part_surface() {
	std::string fileName = "dam_part_surface";
	IGES::IGESReader dam = IGES::IGESReader(input_file(fileName + ".iges"));
	IGES::ShapeManager2DSurfs s2d; dam.getsurfaces2D(s2d);
	IGES::SBFEMExport se(s2d, 0.005,150);
	//se.setLayerColor("8Hvoid mat", -2);
	se.export2vtk(_outputDir, fileName);
	se.export2MFile2D(_outputDir, fileName);
	return;
}
void QA_IGES2D::arcs_conincide_test() {
	std::string fileName = "arcs_test";
	IGES::IGESReader arcs = IGES::IGESReader(input_file(fileName + ".iges"));
	IGES::ShapeManager2DSurfs s2d; arcs.getsurfaces2D(s2d);
	IGES::SBFEMExport se(s2d, 0.0001, 150);
	se.export2vtk(_outputDir, fileName);
	se.export2MFile2D(_outputDir, fileName);
	return;
}
void QA_IGES2D::splines_test() {
	std::string fileName = "splines_test";
	IGES::IGESReader spls = IGES::IGESReader(input_file(fileName + ".iges"));
	IGES::ShapeManager2DSurfs s2d; spls.getsurfaces2D(s2d);
	// debug
	//spls.getsurfaces2D(s2d);
	//std::vector<IGES::IGES_NURBSCurve> cs;
	//s2d.getNurbsCurves(cs);
	//int i = 0;
	//for (IGES::IGES_NURBSCurve _c : cs) {
	//	NURBS::Curve c;
	//	_c.get_curve(c);
	//	std::cout << i++ << "\t" << c.controlPoints[0].x() << ",  " << c.controlPoints[0].y();
	//	std::cout << "\t\t" << c.controlPoints.back().x() << ",  " << c.controlPoints.back().y() << std::endl;
	//}
	// finish debug
	IGES::SBFEMExport se(s2d, 0.001, 150);
	se.export2vtk(_outputDir, fileName);
	se.export2MFile2D(_outputDir, fileName);
	return;
}
void QA_IGES2D::general_test() {
	std::string fileName = "Drawing2";
	IGES::IGESReader gen = IGES::IGESReader(input_file(fileName + ".iges"));
	IGES::ShapeManager2DSurfs s2d; gen.getsurfaces2D(s2d);
	IGES::SBFEMExport se(s2d, 0.001, 150);
	se.export2vtk(_outputDir, fileName);
	se.export2MFile2D(_outputDir, fileName);
	return;
}
void QA_IGES2D::buildings() {
	std::string fileName = "buildings";
	IGES::IGESReader bd = IGES::IGESReader(input_file(fileName + ".iges"));
	IGES::ShapeManager2DSurfs s2d; bd.getsurfaces2D(s2d);
	IGES::SBFEMExport se(s2d, 0.001, 150);
	se.export2vtk(_outputDir, fileName);
	se.export2MFile2D(_outputDir, fileName);
	return;
}
void QA_IGES2D::others() {
	std::string fileName = "L_shape";
	IGES::IGESReader bd = IGES::IGESReader(input_file(fileName + ".iges"));
	IGES::ShapeManager2DSurfs s2d; bd.getsurfaces2D(s2d);
	IGES::SBFEMExport se(s2d, 0.001, 150);
	se.export2vtk(_outputDir, fileName);
	se.export2MFile2D(_outputDir, fileName);
	return;
}