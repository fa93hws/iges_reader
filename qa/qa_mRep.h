#pragma once
#include "..\src\shared\config.h"
#include "..\src\nurbs\basis.h"
#include "..\src\nurbs\surface.h"
#include "qa_template.h"
#include "..\src\io\VTKExport.h"
#include <vector>
#include <armadillo>


class QA_MRep:public Test{
public:
	void calMRep_Test();
	void calMRep_CurveTest();
};
