#include <vector>
#include "../geo/point3d.h"
#include "../iges/shape3D.h"
#include "../iges/IGESReader3D.h"
#include "../nurbs/curve.h"
#include <assert.h> 
#include "iostream"
#include <math.h>
#include <Windows.h>
#include <iomanip>      // std::setprecision
#include "../geo/geometry.h"
#include <thread>

NSI_BEG
class IGESVTKExport3D {
private:

	SolidsGroup _solidsGroup;
	SurfacesGroup _surfacesGroup;
	CurvesGroup _curvesGroup;
	Vertex502 _vertex;

	void discreteSurf(const int idx, std::vector<std::vector<int>>& faces, std::vector<NS::Point3D>& _pts);
	void surfsFileO(const std::string &pathOut, const std::string &fileName,
		const std::vector<std::vector<std::vector<int>>>& faces, std::vector<std::vector<NS::Point3D>>& _pts) const;
public:
	void foo() {};
	void bar(int i) {};
	IGESVTKExport3D(const IGESReader3D &in);
	void exportSurfsToVTK(const std::string &pathOut, std::string &fileName);
};

NS_END