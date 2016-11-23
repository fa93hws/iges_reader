#pragma once
#include "..\src\nurbs\basis.h"
#include "..\src\nurbs\curve.h"
#include "..\src\nurbs\surface.h"
#include "..\src\io\VTKexport.h"
#include "..\src\iges\IGESReader.h"
#include "qa_template.h"
#include "../src/iges/IGESReader.h"
#include "../src/io/SbfemExport.h"
#include "../src/iges/IGESReader3D.h"

class QA_IGES3D :public Test {

public:
	void bullet();
};