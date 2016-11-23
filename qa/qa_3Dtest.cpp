#include "qa_3Dtest.h"
#include <memory>

void QA_IGES3D::bullet() {
	std::string fileName = "bullet";
	IGES::IGESReader3D bullet = IGES::IGESReader3D(input_file(fileName + ".iges"));	
	return;
}