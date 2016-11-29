#include "IGESVTKExport3D.h"
#define __numberThread 10

NSI_BEG
void IGESVTKExport3D::discreteSurf(int idx, std::vector<std::vector<int>> & faces, std::vector<NS::Point3D> & pts) {
	//auto temp = &faces;
	faces.clear();
	pts.clear();
	int numPoints = 51;
	std::shared_ptr<Surface3D> surf;
	_surfacesGroup.getComp(idx, surf);
	surf->discrete(numPoints, faces, pts);
}
void IGESVTKExport3D::surfsFileO(const std::string &pathOut, const std::string &fileName,
	const std::vector<std::vector<std::vector<int>>>& faces, std::vector<std::vector<NS::Point3D>>& _pts) const{
	std::ofstream myfile;
	std::string path = pathOut + "\\";
	CreateDirectory(path.c_str(), NULL);
	myfile.open(path + fileName + ".vtk");
	myfile << "# vtk DataFile Version 3.1" << std::endl; 
	myfile << "Something" << std::endl;
	myfile << "ASCII" << std::endl;
	myfile << "DATASET POLYDATA" << std::endl;

	int numPts = 0;
	for (std::vector<NS::Point3D> pts : _pts) numPts += pts.size();
	myfile << "POINTS\t" << numPts << "\tFLOAT" << std::endl;
	for (std::vector<NS::Point3D> pts : _pts)
		for (NS::Point3D pt : pts) {
			myfile << std::fixed << std::setprecision(5);
			pt.x() > -1e-6 ? myfile << "\t " << pt.x() : myfile << "\t" << pt.x();
			pt.y() > -1e-6 ? myfile << "\t\t " << pt.y() : myfile << "\t\t" << pt.y();
			pt.z() > -1e-6 ? myfile << "\t\t " << pt.z() : myfile << "\t\t" << pt.z();
			myfile << std::endl;
		}
	myfile << std::endl;

	int numCells = 0;
	for (std::vector<std::vector<int>> _faces : faces) numCells += _faces.size();
	myfile << "POLYGONS\t\t" << numCells << "\t\t" << numCells*4 << std::endl;
	int startIdx = 0;
	int currentColor = 0;
	std::vector<int> cellColor;
	for (size_t i = 0; i < faces.size(); i++) {
		std::vector<std::vector<int>> _faces = faces[i];
		for (std::vector<int> _face : _faces) {
			myfile << "\t\t3";
			for (int idx : _face)
				myfile << "\t\t" << idx + startIdx;
			myfile << std::endl;
			cellColor.push_back(currentColor);
		}
		startIdx += _pts[i].size();
		currentColor++;
	}
	myfile << std::endl;
	myfile << std::endl;

	myfile << "CELL_DATA\t" << cellColor.size() << std::endl;
	myfile << "SCALARS\tSurf_Id\tint" << std::endl;
	myfile << "LOOKUP_TABLE default" << std::endl;
	for (int j : cellColor)
		myfile << "\t\t" << j << std::endl;
	myfile.close();
	std::cout << "exported to VTK file" << std::endl;
}

void IGESVTKExport3D::exportSurfsToVTK(const std::string &pathOut, std::string &fileName) {
	static const int nSurf = _surfacesGroup.count();
	std::vector<std::vector<std::vector<int>>> faces;
	faces.resize(nSurf);
	std::vector<std::vector<NS::Point3D>> pts;
	pts.resize(nSurf);

	int idxStart = 0;
	std::thread th[__numberThread];
	do {
		int idxEnd = idxStart + __numberThread;
		bool exitFlag = false;
		if (idxEnd >= nSurf - 1) {
			idxEnd = nSurf-1;
			exitFlag = true;
		}
		for (int i = idxStart; i <= idxEnd; i++)
			th[i] = std::thread([this, i,&faces, &pts] {discreteSurf(i,faces[i],pts[i]); });
		for (int i = 0; i < idxEnd - idxStart + 1; i++)
			th[i].join();
		idxStart += __numberThread;
	} while (idxStart < nSurf - 1);
	surfsFileO(pathOut, fileName, faces, pts);

}
IGESVTKExport3D::IGESVTKExport3D(const IGESReader3D &in) {
	in.getCurvesGroup(_curvesGroup);
	in.getSolidsGroup(_solidsGroup);
	in.getVertexes(_vertex);
	in.getSurfacesGroup(_surfacesGroup);
}

NS_END