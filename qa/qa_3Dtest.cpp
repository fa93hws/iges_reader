#include "qa_3Dtest.h"
#include <memory>

void QA_IGES3D::bullet() {
	std::string fileName = "bullet";
	IGES::IGESReader3D bullet = IGES::IGESReader3D(input_file(fileName + ".iges"));

	IGES::IGESVTKExport3D ve3d = IGES::IGESVTKExport3D(bullet);
	ve3d.exportSurfsToVTK(_outputDir, fileName);
	return;
}

void QA_IGES3D::spinner() {
	std::string fileName = "spinner";
	IGES::IGESReader3D spinner = IGES::IGESReader3D(input_file(fileName + ".iges"));

	IGES::IGESVTKExport3D ve3d = IGES::IGESVTKExport3D(spinner);
	ve3d.exportSurfsToVTK(_outputDir, fileName);

}