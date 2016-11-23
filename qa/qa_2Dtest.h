#pragma once
#include "..\src\nurbs\basis.h"
#include "..\src\nurbs\curve.h"
#include "..\src\nurbs\surface.h"
#include "..\src\io\VTKexport.h"
#include "..\src\iges\IGESReader.h"
#include "qa_template.h"


class QA_IGES2D:public Test {
	
public:
	void dam_part_surface();
	void arcs_conincide_test();
	void splines_test();
	void general_test();
	void buildings();
	void others();
};