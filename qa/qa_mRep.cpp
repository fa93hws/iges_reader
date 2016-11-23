#include "qa_mRep.h"


void QA_MRep::calMRep_Test() {
	std::vector<NURBS::Basis> _b(2);
	std::vector<std::vector<NS::Point3D>> _pts;
	std::vector<std::vector<double>> _w;

	std::vector<double> _btemp = { -1,-1,1,1 };
	_b[0] = NURBS::Basis(_btemp);
	_btemp = { -1,-1,-1,1,1,1 };
	_b[1] = NURBS::Basis(_btemp);

	std::vector<NS::Point3D> _pttemp = { NS::Point3D(1,0,0), NS::Point3D(1,0,1), NS::Point3D(0,0,1) };
	_pts.push_back(_pttemp);
	_pttemp = { NS::Point3D(1,1,0), NS::Point3D(1,1,1), NS::Point3D(0,1,0) };
	_pts.push_back(_pttemp);

	_w = { {1,1,2},{1,1,2} };

	NURBS::Surface surf = NURBS::Surface(_b, _pts, _w);
	//IGES::VTKExport vtk;
	//vtk.setPath("mRep");
	//vtk.setName("Surf01");
	//vtk.addSurface(surf);
	//vtk.generateVTKFile();
	NS::Point3D image;
	surf.getCoord(0.1, 0.5,image);
	double u, v;
	surf.calculatePerImage(image, u, v);

	NS::Point3D pt(0, 0, 0);
	NS::Vector3D V(1, 1, 1);
	std::vector<double> uv;
	std::vector<NS::Point3D> p;
	//surf.calculateLineIntecs(pt, V, uv, p);
	auto s = 1;
}

void QA_MRep::calMRep_CurveTest() {
	
}
