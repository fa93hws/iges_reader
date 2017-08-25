#include <iostream>
#include "qa_2Dtest.h"
#include "qa_mRep.h"
#include "qa_3Dtest.h"

std::string case_directory = "../case/";

int main() {
	QA_IGES2D qa_iges2D;
	//qa_iges2D.set_directory(case_directory ,"dam_part");
	//qa_iges2D.dam_part_surface();
	//qa_iges2D.set_directory(case_directory ,"arcs_test");
	//qa_iges2D.arcs_conincide_test();
	//qa_iges2D.set_directory(case_directory ,"splines_test");
	//qa_iges2D.splines_test();
	//qa_iges2D.set_directory(case_directory ,"general_test");
	//qa_iges2D.general_test();
	//qa_iges2D.set_directory(case_directory, "buildings");
	//qa_iges2D.buildings();
	//qa_iges2D.set_directory(case_directory, "others");
	//qa_iges2D.others();

	//
	//QA_MRep qa_mRep;
	//qa_mRep.calMRep_Test();

	QA_IGES3D qa_iges3D;
	//qa_iges3D.set_directory(case_directory, "bullet");
	//qa_iges3D.bullet();
	qa_iges3D.set_directory(case_directory, "spinner");
	qa_iges3D.spinner();

	getchar();
	return 0;
}