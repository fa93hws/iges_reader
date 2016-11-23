#include "IGESReader.h"

NSI_BEG
double IGESReader::_stod(const std::string &str) {
	int idx = str.find("D");
	if (idx > -1) {
		double d1 = stod(str.substr(0, idx));
		double d2 = stod(str.substr(idx + 1));
		double result = d1 * pow(10.0, d2);
		if (result < 1e-10)
			result = 0;
		return result;
	}
	else
		return stod(str);
}

//void IGESReader::parseEntityByDT(const std::vector<std::string> &dt, Shape shape) {

//}
void  IGESReader::readSPart(const std::string &line) {
	//Start Section
}
void  IGESReader::readGPart(const std::string &line) {
	//Global Section
}
void  IGESReader::readDPart(const std::string &line) {
	//Directory Entry Section
	int idx = d_text.size();
	int line_idx = stoi(line.substr(73,7));
	if (line_idx % 2 == 1)
		d_text.push_back({});
	else
		idx--;
	for (int i = 0; i < 9; i++) {
		d_text[idx].push_back(line.substr(i * 8, 8));
	}
}
void  IGESReader::readPPart(const int pid_begin, const int pid_end, std::vector<std::string> &pt) {
	std::string str;
	for (int i = pid_begin; i < pid_end; i++) {
		str += p_text[i - 1].substr(0, 64);
	}
	//int idx = p_text.size() - 1;
	//if (idx == -1) {
	//	p_text.push_back({});
	//	idx = 0;
	//}
	std::string temp = "";
	for (size_t  i = 0; i < str.size(); i++) {
		if (str[i] == ',') {
			pt.push_back(temp);
			temp = "";
			continue;
		}
		else if (str[i] == ';') {
			pt.push_back(temp);
			return;
		}
		else {
			temp += str[i];
		}
	}
}
void  IGESReader::savePPart(const std::string &line) {
	//Parameter Data Section
	p_text.push_back(line);
}
void  IGESReader::readTPart(const std::string &line) {
	//Terminate Section
}
void IGESReader::reOrganizePart(const std::string &line) {
	//reOrganize Each Part and put in corresponding variable
	char f = line[72]; //flag char
	switch (f) {
	case 'S':
		readSPart(line); break;
	case 'G':
		readGPart(line); break;
	case 'D':
		readDPart(line); break;
	case 'P':
		savePPart(line); break;
	case 'T':
		readTPart(line); break;
	default:
		std::cout << "f not right" << std::endl;
	}

	//int idx;
	//for (int i = 0; i < 5; i++)
	//	if (f == subClass[i]) { idx = i; break; }
	//int endIdx = idx == 3 ? 64 : 72;
	//for (int i = 0; i < endIdx; i++) {
	//	if (i == 0) {
	//		subSection[idx] += line[i]; continue;
	//	}
	//	bool f2 = line[i - 1] == ',' || line[i - 1] == ';' || line[i-1] == ' ';
	//	bool f3 = line[i] == ' ';
	//	if (f2 && f3)
	//		continue;
	//	else
	//		subSection[idx] += line[i];
	//	if (line[i] == ';')
	//		subSection[idx] += "\r\n";
	//}

}
void IGESReader::parseNurbSurface(const std::vector<std::string> &pt, NURBS::ParentSurface &ss) {
	int K0 = stoi(pt[1]);		//nBasis0 - 1
	int K1 = stoi(pt[2]);		//nBasis1 - 1
	int M0 = stoi(pt[3]);		//order0
	int M1 = stoi(pt[4]);		//order1
	bool closed0 = stoi(pt[5]) == 0;	//0-closed?
	bool closed1 = stoi(pt[6]) == 0;	//1-closed?
	bool peri0 = stoi(pt[7]) == 1;		//0-periodic
	bool peri1 = stoi(pt[8]) == 1;		//1-periodic

	int N0 = 1 + K0 - M0;		// n Interior Span 0;
	int N1 = 1 + K1 - M1;		// n Interior Span 1;
	int A = N0 + 2 * M0;		// n knotVector 0;
	int B = N1 + 2 * M1;		// n knotVector 1;
	int C = (1 + K0) * (1 + K1);// n Control Points (Total);
								//parse knotVectors
	std::vector<std::vector<double>> knotVector;
	knotVector.push_back({});
	for (int i = 10; i <= 10 + A; i++)
		knotVector[0].push_back(_stod(pt[i]));
	knotVector.push_back({});
	for (int i = 11 + A; i <= 11 + A + B; i++)
		knotVector[1].push_back(_stod(pt[i]));
	// parse weightVector
	std::vector<std::vector<double>> weightVector(K1 + 1);
	for (int i = 0; i < K1 + 1; i++)
		for (int j = 0; j < K0 + 1; j++)
			weightVector[i].push_back(_stod(pt[12 + A + B + i * (K0 + 1) + j]));
	weightVector = transpose(weightVector);
	//parse controlPoints
	std::vector<std::vector<NS::Point3D>> controlPoints(K1 + 1);
	int nPt = 0;
	for (int i = 0; i < K1 + 1; i++)
		for (int j = 0; j < K0 + 1; j++) {
			double x = _stod(pt[12 + A + B + C + nPt * 3]);
			double y = _stod(pt[12 + A + B + C + nPt * 3 + 1]);
			double z = _stod(pt[12 + A + B + C + nPt * 3 + 2]);
			nPt++;
			controlPoints[i].push_back(NS::Point3D(x, y, z));
		}
	controlPoints = transpose(controlPoints);
	//parse start-end info
	std::vector<double> u_begin; u_begin.push_back(_stod(pt[12 + A + B + 4 * C])); u_begin.push_back(_stod(pt[12 + A + B + 4 * C + 2]));
	std::vector<double> u_end; u_end.push_back(_stod(pt[12 + A + B + 4 * C + 1])); u_end.push_back(_stod(pt[12 + A + B + 4 * C + 3]));
	// get nurbsurface
	std::vector<NURBS::Basis> b(2);
	for (int i = 0; i < 2; i++) {
		b[i] = NURBS::Basis(knotVector[i]);
		b[i].setRange(u_begin[i], u_end[i]);
	}
	ss = NURBS::ParentSurface(b, controlPoints, weightVector);
}
void IGESReader::parse_matrix_tran(const std::vector<std::string> &dt, std::vector<std::vector<double>> &R, std::vector<double> &T) {
	R.clear();
	T.clear();
	std::string empty = "        ";
	int pid_begin = stoi(dt[1]);
	int pid_end = pid_begin + stoi(dt[12]);
	int matrix_tran = dt[6] == empty ? 0 : stoi(dt[6]);
	std::vector<std::string> pt;
	readPPart(pid_begin, pid_end, pt);
	R.resize(3);
	for (int i = 0; i < 3; i++) R[i].resize(3);
	T.resize(3);
	R[0][0] = _stod(pt[1]); R[0][1] = _stod(pt[2]); R[0][2] = _stod(pt[3]); T[0] = _stod(pt[4]);
	R[1][0] = _stod(pt[5]); R[1][1] = _stod(pt[6]); R[1][2] = _stod(pt[7]); T[1] = _stod(pt[8]);
	R[2][0] = _stod(pt[9]); R[2][1] = _stod(pt[10]); R[2][2] = _stod(pt[11]); T[2] = _stod(pt[12]);
}
//void IGESReader::set_shape_attributes(const Attributes &attributes,const Attributes &parentAttributes, const StatusNumber &status, const StatusNumber &parentStatus, Shape &shape) {
//	if (parentStatus.hierarchy == 0)
//		shape.color_id = parentAttributes.color_idx;
//	else if (parentStatus.hierarchy == 1)
//		shape.color_id = attributes.color_idx;
//	else
//		ASSERT(0);
//}
void IGESReader::set_attribute(const std::vector<std::string> &dt, Shape &s) {
	StatusNumber status(dt[8]);
	if (status.hierarchy == 0)	//use parent's attribute				
		s.setpropHierarchy(true);
	std::string empty = "        ";
	int color = dt[11] == empty ? 0 : stoi(dt[11]);
	if (color == 0)
		return;
	ASSERT(color < 0);
	int d_idx = (-color - 1) / 2;
	std::vector<std::string> dt_color = d_text[d_idx];
	ASSERT(stoi(dt_color[13]) == 0);
	int pid_begin = stoi(dt_color[1]);
	int pid_end = pid_begin + stoi(dt_color[12]);
	std::vector<std::string> pt;
	readPPart(pid_begin, pid_end, pt);
	s.set_color(_stod(pt[1]), _stod(pt[2]), _stod(pt[3]));

	int level = dt[4] == empty ? 0 : stoi(dt[4]);
	s.set_layer(layer_name[ stoi(dt[4]) ]);
	s.set_level(level);
}
void IGESReader::set_attribute(const std::vector<std::string> &dt, ShapeManager2DSurfs &s, const int surf_id) {
	std::string empty = "        ";
	int level = dt[4] == empty ? 0 : stoi(dt[4]);
	s.set_layer(layer_name[stoi(dt[4])],surf_id);
	s.set_level(level,surf_id);
}
void IGESReader::parseCircularArc(const std::vector<std::string> &dt, CircularArc &cc) {
	std::string empty = "        ";
	int pid_begin = stoi(dt[1]);
	int pid_end = pid_begin + stoi(dt[12]);
	int matrix_tran = dt[6] == empty ? 0 : stoi(dt[6]);
	std::vector<std::string> pt;
	readPPart(pid_begin, pid_end, pt);
	double z = _stod(pt[1]), x = _stod(pt[2]), y = _stod(pt[3]);
	double x2 = _stod(pt[4]), y2 = _stod(pt[5]), x3 = _stod(pt[6]), y3 = _stod(pt[7]);
	NS::Point3D pt1(x, y, z), pt2(x2, y2, z), pt3(x3, y3, z);
	cc = CircularArc(pt1, pt2, pt3);
	if (matrix_tran > 0) {
		std::vector<std::vector<double>> R; std::vector<double> T;
		parse_matrix_tran(d_text[(matrix_tran - 1) / 2], R, T);
		cc.set_tran_matrix(R, T);
	}
	set_attribute(dt, cc);
}
void IGESReader::parseNurbCurve(const std::vector<std::string> &dt, IGES_NURBSCurve &cc) {
	std::string empty = "        ";
	int pid_begin = stoi(dt[1]);
	int matrix_tran = dt[6] == empty ? 0 : stoi(dt[6]);		
	int pid_end = pid_begin + stoi(dt[12]);
	int form = stoi(dt[13]);
	std::vector<std::string> pt;
	readPPart(pid_begin, pid_end, pt);

	int  K = stoi(pt[1]);	//nBasis - 1
	int  M = stoi(pt[2]);	//order
	bool planner = stoi(pt[3]) == 1;	//planar?
	bool closed = stoi(pt[4]) == 1;	//closed?
	bool rational = stoi(pt[5]) == 0;	//rational?
	bool periodic = stoi(pt[6]) == 0;	//periodic
	int N = 1 + K - M;					//n interior span
	int A = N + 2 * M;					//n knotVector
										//parse knotVector
	std::vector<double> knotVector;
	for (int i = 7; i <= 7 + A; i++)
		knotVector.push_back(_stod(pt[i]));
	//parse weightVector
	std::vector<double> weightVector;
	if (rational)
		for (int i = 8 + A; i <= 8 + A + K; i++)
			weightVector.push_back(_stod(pt[i]));
	// parse control point
	std::vector<NS::Point3D> controlPoints;
	for (int i = 9 + A + K; i <= 9 + A + 4 * K; i += 3)
		controlPoints.push_back(NS::Point3D(_stod(pt[i]), _stod(pt[i + 1]), _stod(pt[i + 2])));
	// start and end position
	double u_begin = _stod(pt[12 + A + 4 * K]);
	double u_end = _stod(pt[13 + A + 4 * K]);
	//
	NURBS::Basis b = NURBS::Basis(knotVector);
	b.setRange(u_begin, u_end);
	if (rational)
		cc = IGES_NURBSCurve(b, controlPoints, weightVector);
	else
		cc = IGES_NURBSCurve(b, controlPoints);
	set_attribute(dt, cc);
}
void IGESReader::parseLine(const std::vector<std::string> &dt, IGES::Line &line) {
	std::string empty = "        ";
	int pid_begin	= stoi(dt[1]);
	int matrix_tran = dt[ 6] == empty ? 0 : stoi(dt[ 6]);
	int pid_end		= pid_begin + stoi(dt[12]);
	int form		= stoi(dt[13]);
	if (form != 0)
		std::cout << "form 0 of LINE(110) is not supported" << std::endl;
	std::vector<std::string> pt;
	readPPart(pid_begin, pid_end, pt);
	double x0 = _stod(pt[1]), y0 = _stod(pt[2]), z0 = _stod(pt[3]);
	double x1 = _stod(pt[4]), y1 = _stod(pt[5]), z1 = _stod(pt[6]);
	NS::Point3D pt1(x0, y0, z0), pt2(x1, y1, z1);
	line = IGES::Line(pt1, pt2);
	if (matrix_tran > 0) {
		std::vector<std::vector<double>> R; std::vector<double> T;
		parse_matrix_tran(d_text[(matrix_tran - 1) / 2], R, T);
		line.set_tran_matrix(R, T);
	}
	set_attribute(dt, line);
}
void IGESReader::parseSubFigure408(const std::vector<std::string> &dt, ShapeGroup &sg) {
	std::string empty = "        ";
	int pid_begin = stoi(dt[1]);
	int pid_end = pid_begin + stoi(dt[12]);
	int matrix_tran = dt[6] == empty ? 0 : stoi(dt[6]);
	std::vector<std::string> pt;
	readPPart(pid_begin, pid_end, pt);
	int DE_pointer = (stoi(pt[1]) - 1)/2;
	double x = _stod(pt[2]), y = _stod(pt[3]), z = _stod(pt[3]), S = _stod(pt[4]);
	std::vector<std::vector<double>> R = { { 1,0,0 },{ 0,1,0 },{ 0,0,1 } };
	std::vector<double> T = { 0,0,0 };
	if (matrix_tran > 0)
		parse_matrix_tran(d_text[(matrix_tran - 1) / 2], R, T);


	std::vector<std::string> dt_real = d_text[DE_pointer];				// defination of the group of the geometry idx
	ASSERT(stoi(dt_real[0]) == 308);
	pid_begin = stoi(dt_real[1]);
	pid_end = pid_begin + stoi(dt_real[12]);
	pt.clear();
	readPPart(pid_begin, pid_end, pt);
	int n_member = stoi(pt[3]);
	for (int i = 0; i < n_member; i++) {
		int member_idx = (stoi(pt[i + 4]) - 1) / 2;
		std::vector<std::string> dt_member = d_text[member_idx];
		int entityTyp = stoi(dt_member[0]);
		switch (entityTyp) {
		case 100: {										// Circular Arc
			CircularArc arc;
			parseCircularArc(dt_member, arc);
			arc.set_tran_matrix(R, T);
			sg.add_arc(arc);
			//return;
			break;
		}
		case 110: {										// Line
			Line line;
			parseLine(dt_member, line);
			line.set_tran_matrix(R, T);
			sg.add_line(line);
			break;
		}
		case 126: {
			IGES_NURBSCurve nurbs_curve;
			parseNurbCurve(dt_member,  nurbs_curve);
			nurbs_curve.set_tran_matrix(R, T);
			sg.add_nurbs_curves(nurbs_curve);
			break;
		}
		case 408: {										// Singular Subfigure Instance
			ShapeGroup _shapes;
			parseSubFigure408(dt_member,  _shapes);
			break;
		}
		default:
			std::cout << "entity typ number " << entityTyp << " not supported in parseSubFigure408" << std::endl;
		}
		//Shape s = parseEntityByDT(dt_member);
	}
	set_attribute(dt, sg);
}
void IGESReader::parseSurface2D143(const std::vector<std::string> &dt, ShapeManager2DSurfs &surface2D) {
	std::string empty = "        ";
	int pid_begin = stoi(dt[1]);
	int pid_end = pid_begin + stoi(dt[12]);
	int matrix_tran = dt[6] == empty ? 0 : stoi(dt[6]);
	std::vector<std::string> pt;
	readPPart(pid_begin, pid_end, pt);
	assert(stoi(pt[3]) == 1);
	int DE_pointer = (stoi(pt[4]) - 1) / 2;
	double x = _stod(pt[2]), y = _stod(pt[3]), z = _stod(pt[3]), S = _stod(pt[4]);
	std::vector<std::vector<double>> R = { {1,0,0},{0,1,0},{0,0,1} };
	std::vector<double> T = {0,0,0};
	if (matrix_tran > 0) 
		parse_matrix_tran(d_text[(matrix_tran - 1) / 2], R, T);
	std::vector<std::string> dt_real = d_text[DE_pointer];				// defination of the group of the geometry idx
	assert(stoi(dt_real[0]) == 141);
	pid_begin = stoi(dt_real[1]);
	pid_end = pid_begin + stoi(dt_real[12]);
	pt.clear();
	readPPart(pid_begin, pid_end, pt);
	int n_member = stoi(pt[4]);
	int surf_id = surface2D.getSurfSize();
	surface2D.addSurf();
	for (int i = 0; i < n_member; i++) {
		int start_p_idx = 5 + 3 * i;
		int member_idx = (stoi(pt[start_p_idx]) - 1) / 2;
		int dir = stoi(pt[start_p_idx + 1])==1 ? 1 : -1;
		assert(stoi(pt[start_p_idx + 2]) == 0);
		std::vector<std::string> dt_member = d_text[member_idx];
		int entityTyp = stoi(dt_member[0]);
		switch (entityTyp) {
		case 100: {										// Circular Arc
			CircularArc arc;
			parseCircularArc(dt_member, arc);
			arc.set_tran_matrix(R, T);
			surface2D.addArc(surf_id,arc,dir);
			//surface2D.add_arc(arc);
			break;
		}
		case 110: {										// Line
			Line line;
			parseLine(dt_member, line);
			line.set_tran_matrix(R, T);
			surface2D.addLine(surf_id, line, dir);
			//line.set_dir(dir);
			//surface2D.add_line(line);
			break;
		}
		case 126: {
			IGES_NURBSCurve nurbs_curve;
			parseNurbCurve(dt_member, nurbs_curve);
			nurbs_curve.set_tran_matrix(R, T);
			//nurbs_curve.set_dir(dir);
			surface2D.addNURBSCurve(surf_id,nurbs_curve,dir);
			break;
		}
		default:
			std::cout << "entity typ number " << entityTyp << " not supported in parseSurface2D-143" << std::endl;
		}
		//Shape s = parseEntityByDT(dt_member);
	}
	set_attribute(dt, surface2D, surf_id);
}
void IGESReader::parseProp406_layer(const std::vector<std::string> &dt, std::vector<std::string> &layer_name) {
	if (stoi(dt[13]) != 3)
		return;
	int level = stoi(dt[4]);
	int pid_begin = stoi(dt[1]);
	int pid_end = pid_begin + stoi(dt[12]);
	std::vector<std::string> pt;
	readPPart(pid_begin, pid_end, pt);
	ASSERT(stoi(pt[1]) == 2);
	layer_name[level] = pt[3];
}
	

IGESReader::IGESReader(std::string &path) {
	//
	std::cout << "start parsing iges input" << std::endl;
	// init variable
	std::string line;
	myfile = std::ifstream(path);
	if (myfile.is_open()) {
		while (getline(myfile, line))
		{
			reOrganizePart(line);
		}
		myfile.close();
	}
	else {
		std::cout << "unable to open file"<<std::endl;
		std::cout << "file path: " + path << std::endl;
		abort();
	}
	parseAllShapes();
	if (_surfs2D.getSurfSize() > 0) { // surface mode
		std::cout << "Input is surfaces" << std::endl;
		_surfs2D.removeConincideLines();
		_surfs2D.removeConincideArcs();
		_surfs2D.removeConcincideNURBSCrrves();
	}
	std::cout << "Reading all geometries successfully" << std::endl;
}
void IGESReader::parseAllShapes() {
	//110(line), 110(arc) and 126(BSpline) are extracted here
	layer_name.resize(100);
	for (size_t i = 0; i < d_text.size(); i++) {
		std::vector<std::string> dt = d_text[i];
		if (stoi(dt[0]) == 406)
			parseProp406_layer(dt, layer_name);
	}
	for (size_t i = 0; i < d_text.size(); i++) {
		std::vector<std::string> dt = d_text[i];
		if (dt[8].substr(2, 2) != "00")					// not independent
			continue;
		int entityTyp = stoi(dt[0]);
		switch (entityTyp) {
		case 143: {
			//Surface2D surf;
			parseSurface2D143(dt, _surfs2D);
			//shape_collection.addSurface(surf);
			break;
		}
		case 100: {
			CircularArc cc;
			parseCircularArc(dt,cc);
			shape_collection.addArc(cc); 
			break;
		}
		case 110: {
			Line line;
			parseLine(dt, line);
			shape_collection.addLine(line);
			break;
		}
		case 126: {
			IGES_NURBSCurve curve;
			parseNurbCurve(dt, curve);
			shape_collection.addNURBSCurve(curve);
			break;
		}			
		case 408: {
			ShapeGroup sg;
			parseSubFigure408(dt, sg);
			shape_collection.addShapeGroup(sg);
			break;
		}
		case 406:
		case 314:
		case 124:
			break;
		default:
			std::cout << "entity typ number " << entityTyp << " not supported" << std::endl;
		}
	}
}
//void IGESReader::parseEntities() {
//	Attributes empty_attributes;
//	StatusNumber init_status("00010001");
//	for (size_t i = 0; i < d_text.size(); i++) {
//		std::vector<std::string> dt = d_text[i];
//		if (dt[8].substr(2, 2) != "00")					// not independent
//			continue;
//		int entityTyp = stoi(dt[0]);
//		switch (entityTyp) {
//		case 100: {										// Circular Arc
//			CircularArc arc;
//			parseCircularArc(dt, init_status, empty_attributes, arc);
//			//arc.name = "arc";
//			shape_collection.addArc(arc);
//			break;
//		}
//		case 110: {										// Line
//			Line line;
//			parseLine(dt, init_status, empty_attributes, line);
//			//line.name = "line";
//			shape_collection.addLine(line);
//			break;
//		}
//		case 126: {
//			IGES_NURBSCurve nurbs_curve;
//			parseNurbCurve(dt, init_status, empty_attributes, nurbs_curve);
//			shape_collection.addNURBSCurve(nurbs_curve);
//			break;
//		}
//		case 408: {										// Singular Subfigure Instance
//			std::vector<Shape> shapes;
//			parseSubFigure408(dt, init_status, empty_attributes);
//			break;
//		}
//		case 124:										// matrix transformation, continue
//		case 308:										// Singular Subfigure defineation, use instance instead, continue
//		case 406:										// Property, continue
//			break;
//		default:
//			std::cout << "entity typ number " << entityTyp << " not supported" << std::endl;
//		}
//	}
//}

void IGESReader::parseEntities(std::vector<NURBS::ParentSurface> &s3D, std::vector<NURBS::Curve> &c2D, std::vector<IGES::Shape> &shp) {
	return;
	//for (size_t i = 0; i < p_text.size(); i++) {
	//	std::vector<std::string> pt = p_text[i];
	//	std::vector<std::string> dt = d_text[i];
	//	if (pt.size() == 0)
	//		continue;
	//	int entityTyp = stoi(pt[0]);
	//	switch (entityTyp) {
	//	case 128: {									// NURB Surface
	//		NURBS::ParentSurface s = NURBS::ParentSurface();
	//		parseNurbSurface(pt, s);
	//		s.id = i * 2 + 1;
	//		s3D.push_back(s);
	//		break;
	//	}
	//	case 126: {
	//		NURBS::Curve c;
	//		parseNurbCurve(pt, c);
	//		c.id = i * 2 + 1;
	//		c2D.push_back(c);
	//		break;
	//	}
	//	default:
	//		std::cout << "entity typ number " << entityTyp << " not supported" << std::endl;
	//	}
	//}
}
NS_END