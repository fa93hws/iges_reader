#include "IGESReader3D.h"


NSI_BEG
double IGESReader3D::_stod(const std::string &str) {
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
/**************************************** D,P,part **************************************************/
void IGESReader3D::readPPart(const std::string &line, std::vector<std::string> &p_text) const{
	p_text.push_back(line);
}
void IGESReader3D::readDPart(const std::string &line, std::vector<std::vector<std::string>> &d_text) const {
	//Directory Entry Section
	int idx = d_text.size();
	int line_idx = stoi(line.substr(73, 7));
	if (line_idx % 2 == 1)
		d_text.push_back({});
	else
		idx--;
	for (int i = 0; i < 9; i++) {
		d_text[idx].push_back(line.substr(i * 8, 8));
	}
}
void IGESReader3D::getPPart(const int pid_begin, const int pid_end, const std::vector<std::string> &p_text,std::vector<std::string> &pt) const {
	std::string str;
	for (int i = pid_begin; i < pid_end; i++) {
		str += p_text[i - 1].substr(0, 64);
	}
	std::string temp = "";
	for (size_t i = 0; i < str.size(); i++) {
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
/**************************************** Parse Shapes **************************************************/
void IGESReader3D::parseAllShapes(std::vector<std::vector<std::string>> &d_text, const std::vector<std::string> &pt, SolidsGroup &sg) {
	for (size_t i = 0; i < d_text.size(); i++) {
		std::vector<std::string> dt = d_text[i];
		if (dt[8].substr(2, 2) != "00")					// not independent
			continue;
		int entityTyp = stoi(dt[0]);
		switch (entityTyp) {
		case 186: {
			Solid186 solid;
			parseSolid186(d_text, i, pt, solid);
			_solidsGroup.addComp(solid,i);
		}
		case 314: // color
		case 124: // transformation matrix
		case 406: // properties
			break;
		default:
			std::cout << "entity typ number " << entityTyp << " not supported in parseAllShapes" << std::endl;
		}
	}
}
void IGESReader3D::parseSolid186(_PARSEINPUT, Solid186 &solid) {
	// p.229
	//int matrix_tran = dt[6] == empty ? 0 : stoi(dt[6]);
	std::vector<std::string> dt = dt_all[dt_idx];
	int pid_begin = stoi(dt[1]);
	int pid_end = pid_begin + stoi(dt[12]);
	int matrixTran = dt[6] == _empty ? 0 : stoi(dt[6]);
	std::vector<std::string> mypt;
	getPPart(pid_begin, pid_end ,pt,mypt);

	assert(stoi(mypt[3]) == 0);//no void solid
	int shell_pointer = (stoi(mypt[1]) - 1) / 2;
	Shell514 shell;
	parseShell514(dt_all, shell_pointer, pt, shell);
	solid.setShell(shell);
	solid.setDir(stoi(mypt[2]));
	parseMatrixTran(dt_all, matrixTran, pt, std::make_shared<Shape3D>(solid));
	// void shell
}
void IGESReader3D::parseShell514(_PARSEINPUT, Shell514 &shell) {
	// p.551
	std::vector<std::string> dt = dt_all[dt_idx];
	int pid_begin = stoi(dt[1]);
	int pid_end = pid_begin + stoi(dt[12]);	
	int formid = stoi(dt[13]);
	assert(formid == 1);//close shell
	std::vector<std::string> mypt;
	getPPart(pid_begin, pid_end, pt, mypt);
	
	int num_shell = stoi(mypt[1]);
	for (int i = 0; i < num_shell; i++) {
		Face510 face;
		int faceid	= stoi(mypt[i * 2 + 2]);
		faceid = (faceid - 1) / 2;
		int dir		= stoi(mypt[i * 2 + 3]);
		parseFace510(dt_all, faceid,pt,face);
		shell.addFace(face, dir);
	}
}
void IGESReader3D::parseFace510(_PARSEINPUT,Face510 &face) {
	// p.550
	std::vector<std::string> dt = dt_all[dt_idx];
	int pid_begin = stoi(dt[1]);
	int pid_end = pid_begin + stoi(dt[12]);
	std::vector<std::string> mypt;
	getPPart(pid_begin, pid_end, pt, mypt);

	int surfid = stoi(mypt[1]);
	surfid = (surfid - 1) / 2;
	face.setSurf(parseSurf(dt_all, surfid, pt, _surfaceGroup));
	// to be continue
	int numloops = stoi(mypt[2]);
	int outerloop = stoi(mypt[3]);
	face.setDir(outerloop);
	for (int i = 0; i < numloops; i++) {
		Loop508 loop;
		int loopid = stoi(mypt[4 + i]);
		loopid = (loopid - 1) / 2;
		parseLoop508(dt_all, loopid, pt, loop);
		face.addLoop(loop);
	}
}
int IGESReader3D::parseSurf(_PARSEINPUT, SurfacesGroup &sg) {
	std::vector<std::string> dt = dt_all[dt_idx];
	int dirIdx = sg.findDirIdx(dt_idx);
	if (dirIdx > -1) return dirIdx;
	int typ = stoi(dt[0]);
	int matrixTran = dt[6] == _empty ? 0 : stoi(dt[6]);
	switch (typ) {
	case 128: {
		Surf128* pss = new Surf128();
		parseNURBSSurf128(dt_all, dt_idx, pt, *pss);
		std::shared_ptr<Surface3D> sharedSS = std::shared_ptr<Surface3D>(pss);
		sg.addComp(sharedSS, dt_idx);
		std::shared_ptr<Shape3D> SharedS3 = std::shared_ptr<Shape3D>(sharedSS);
		parseMatrixTran(dt_all, matrixTran, pt, SharedS3);
		break;
	}
	default: 
		std::cout << "type #" << typ << " is not supported in parseSurf" << std::endl;
		assert(0);
	}
	dirIdx = sg.count();
	assert(dirIdx > -1);
	return dirIdx;
}
void IGESReader3D::parseNURBSSurf128(_PARSEINPUT, Surf128 &ss) {
	std::vector<std::string> dt = dt_all[dt_idx];
	int pid_begin = stoi(dt[1]);
	int pid_end = pid_begin + stoi(dt[12]);
	std::vector<std::string> txt;
	getPPart(pid_begin, pid_end, pt, txt);

	int K0 = stoi(txt[1]);		//nBasis0 - 1
	int K1 = stoi(txt[2]);		//nBasis1 - 1
	int M0 = stoi(txt[3]);		//order0
	int M1 = stoi(txt[4]);		//order1
	bool closed0 = stoi(txt[5]) == 0;	//0-closed?
	bool closed1 = stoi(txt[6]) == 0;	//1-closed?
	bool peri0 = stoi(txt[7]) == 1;		//0-periodic
	bool peri1 = stoi(txt[8]) == 1;		//1-periodic

	int N0 = 1 + K0 - M0;		// n Interior Span 0;
	int N1 = 1 + K1 - M1;		// n Interior Span 1;
	int A = N0 + 2 * M0;		// n knotVector 0;
	int B = N1 + 2 * M1;		// n knotVector 1;
	int C = (1 + K0) * (1 + K1);// n Control Points (Total);
								//parse knotVectors
	std::vector<std::vector<double>> knotVector;
	knotVector.push_back({});
	for (int i = 10; i <= 10 + A; i++)
		knotVector[0].push_back(_stod(txt[i]));
	knotVector.push_back({});
	for (int i = 11 + A; i <= 11 + A + B; i++)
		knotVector[1].push_back(_stod(txt[i]));
	// parse weightVector
	std::vector<std::vector<double>> weightVector(K0 + 1);
	for (int i = 0; i < K1 + 1; i++)
		for (int j = 0; j < K0 + 1; j++)
			weightVector[j].push_back(_stod(txt[12 + A + B + i * (K0 + 1) + j]));
	//weightVector = transpose(weightVector);
	//parse controlPoints
	std::vector<std::vector<NS::Point3D>> controlPoints(K0 + 1);
	int nPt = 0;
	for (int i = 0; i < K1 + 1; i++)
		for (int j = 0; j < K0 + 1; j++) {
			double x = _stod(txt[12 + A + B + C + nPt * 3]);
			double y = _stod(txt[12 + A + B + C + nPt * 3 + 1]);
			double z = _stod(txt[12 + A + B + C + nPt * 3 + 2]);
			nPt++;
			controlPoints[j].push_back(NS::Point3D(x, y, z));
		}
	//controlPoints = transpose(controlPoints);
	//parse start-end info
	std::vector<double> u_begin; u_begin.push_back(_stod(txt[12 + A + B + 4 * C])); u_begin.push_back(_stod(txt[12 + A + B + 4 * C + 2]));
	std::vector<double> u_end; u_end.push_back(_stod(txt[12 + A + B + 4 * C + 1])); u_end.push_back(_stod(txt[12 + A + B + 4 * C + 3]));
	// get nurbsurface
	std::vector<NURBS::Basis> b(2);
	for (int i = 0; i < 2; i++) {
		b[i] = NURBS::Basis(knotVector[i]);
		b[i].setRange(u_begin[i], u_end[i]);
	}
	NURBS::ParentSurface surf(b, controlPoints, weightVector);
	ss.setSurf(surf);
}
void IGESReader3D::parseLoop508(_PARSEINPUT, Loop508 &loop) {
	std::vector<std::string> dt = dt_all[dt_idx];
	int pid_begin = stoi(dt[1]);
	int pid_end = pid_begin + stoi(dt[12]);
	int formid = stoi(dt[13]);
	assert(formid == 1); // Bound of surface
	std::vector<std::string> mypt;
	getPPart(pid_begin, pid_end, pt, mypt);
	int numLoop = stoi(mypt[1]);
	auto ptr = mypt.cbegin()+1;
	for (int i = 0; i < numLoop; i++) {
		int typ = stoi(*(++ptr));
		assert(typ == 0); // only edge is considered
		int edgePtr = stoi(*(++ptr));
		edgePtr = (edgePtr - 1 ) / 2;
		int edgeIdx = stoi(*(++ptr));
		int dir		= stoi(*(++ptr));
		assert(stoi(*(++ptr)) == 0);
		Edge504 edge;
		parseEdge504(dt_all, edgePtr, pt, edgeIdx - 1, edge);
		loop.addEdge(edge,dir);
	}
}
int IGESReader3D::parseCurve(_PARSEINPUT, CurvesGroup &cg) {
	std::vector<std::string> dt = dt_all[dt_idx];
	int dirIdx = cg.findDirIdx(dt_idx);
	if (dirIdx > -1)
		return dirIdx;

	int typ = stoi(dt[0]);
	int matrixTran = dt[6] == _empty ? 0 : stoi(dt[6]);
	switch (typ) {
	case 126: {
		Curve126* pcc = new Curve126();
		parseCurve126(dt_all, dt_idx, pt, *pcc);
		std::shared_ptr<Curve3D> sharedCC = std::shared_ptr<Curve3D>(pcc);
		cg.addComp(sharedCC, dt_idx);
		std::shared_ptr<Shape3D> sharedC3 = std::shared_ptr<Shape3D>(sharedCC);
		parseMatrixTran(dt_all, matrixTran, pt, sharedC3);
		break;
	}
	default:
		std::cout << "type #" << typ << "not supported in parseCurve" << std::endl;
		assert(0);
	}
	dirIdx = cg.count();
	assert(dirIdx > -1);
	return dirIdx;
}
void IGESReader3D::parseCurve126(_PARSEINPUT, Curve126 &cc) {
	std::vector<std::string> dt = dt_all[dt_idx];
	int pid_begin = stoi(dt[1]);
	int matrix_tran = dt[6] == _empty ? 0 : stoi(dt[6]);
	int pid_end = pid_begin + stoi(dt[12]);
	int form = stoi(dt[13]);
	std::vector<std::string> mypt;
	getPPart(pid_begin, pid_end, pt,mypt);

	int  K = stoi(mypt[1]);	//nBasis - 1
	int  M = stoi(mypt[2]);	//order
	bool planner = stoi(mypt[3]) == 1;	//planar?
	bool closed = stoi(mypt[4]) == 1;	//closed?
	bool rational = stoi(mypt[5]) == 0;	//rational?
	bool periodic = stoi(mypt[6]) == 0;	//periodic
	int N = 1 + K - M;					//n interior span
	int A = N + 2 * M;					//n knotVector
										//parse knotVector
	std::vector<double> knotVector;
	for (int i = 7; i <= 7 + A; i++)
		knotVector.push_back(_stod(mypt[i]));
	//parse weightVector
	std::vector<double> weightVector;
	if (rational)
		for (int i = 8 + A; i <= 8 + A + K; i++)
			weightVector.push_back(_stod(mypt[i]));
	// parse control point
	std::vector<NS::Point3D> controlPoints;
	for (int i = 9 + A + K; i <= 9 + A + 4 * K; i += 3)
		controlPoints.push_back(NS::Point3D(_stod(mypt[i]), _stod(mypt[i + 1]), _stod(mypt[i + 2])));
	// start and end position
	double u_begin = _stod(mypt[12 + A + 4 * K]);
	double u_end = _stod(mypt[13 + A + 4 * K]);
	//
	NURBS::Basis b = NURBS::Basis(knotVector);
	b.setRange(u_begin, u_end);
	NURBS::ParentCurve pc;
	if (rational)
		pc = NURBS::ParentCurve(b, controlPoints, weightVector);
	else
		pc = NURBS::ParentCurve(b, controlPoints);
	cc.setCurve(pc);
}
int IGESReader3D::parseVer502(_PARSEINPUT, const int &ptidx, Vertex502 &ver) {
	// ptidx start from 0 instead of 1 (different from that in standard)
	std::vector<std::string> dt = dt_all[dt_idx];
	std::pair<int, int> idx = { dt_idx,ptidx };
	int dirIdx = ver.findDirIdx(idx);
	if (dirIdx > -1) return dirIdx;

	int pid_begin = stoi(dt[1]);
	int pid_end = pid_begin + stoi(dt[12]);
	std::vector<std::string> mypt;
	getPPart(pid_begin, pid_end, pt, mypt);
	auto ptr = mypt.begin() + 3 * ptidx + 2;
	double x = _stod(*(ptr++));
	double y = _stod(*(ptr++));
	double z = _stod(*(ptr++));
	ver.addComp(NS::Point3D(x, y, z), idx);
	dirIdx = ver.count();
	assert(dirIdx > -1);
	return dirIdx;
}
void IGESReader3D::parseEdge504(_PARSEINPUT, const int edgeIdx, Edge504 &edge) {
	// p.546
	// edgeIdx start from 0 instead of 1 (different from that in standard)
	std::vector<std::string> dt = dt_all[dt_idx];
	int pid_begin = stoi(dt[1]);
	int pid_end = pid_begin + stoi(dt[12]);
	std::vector<std::string> mypt;
	getPPart(pid_begin, pid_end, pt, mypt);

	assert(edgeIdx <= stoi(mypt[1]));
	int curvePtr = stoi(mypt[(edgeIdx) * 5 + 2]);
	curvePtr = (curvePtr - 1) / 2;
	int curveId = parseCurve(dt_all, curvePtr, pt,_curvesGroup);
	edge.setCurve(curveId);
	std::pair<int,int> ptPtr(stoi(mypt[edgeIdx * 5 + 3]) ,stoi(mypt[edgeIdx * 5 + 4])) ;
	ptPtr.first = (ptPtr.first - 1) / 2;
	int pt0 = parseVer502(dt_all, ptPtr.first, pt, ptPtr.second - 1,_vertex);
	ptPtr = std::make_pair(stoi(mypt[edgeIdx * 5 + 5]) ,stoi(mypt[edgeIdx * 5 + 6]));
	ptPtr.first = (ptPtr.first - 1) / 2;
	int pt1 = parseVer502(dt_all, ptPtr.first, pt, ptPtr.second - 1, _vertex);
	edge.setPts(pt0, pt1);
}
void IGESReader3D::parseMatrixTran(_PARSEINPUT, std::shared_ptr<Shape3D> &s3d) {
	arma::mat R(3, 3, arma::fill::eye);
	arma::vec T(3, 1, arma::fill::zeros);
	std::vector<std::string> dt = dt_all[dt_idx];
	int pid_begin = stoi(dt[1]);
	int pid_end = pid_begin + stoi(dt[12]);
	int matrix_tran = dt[6] == _empty ? 0 : stoi(dt[6]);
	std::vector<std::string> mypt;
	getPPart(pid_begin, pid_end, pt, mypt);
	R[0,0] = _stod(pt[1]); R[0,1] = _stod(pt[2]);  R[0,2] = _stod(pt[3]);  T[0] = _stod(pt[4]);
	R[1,0] = _stod(pt[5]); R[1,1] = _stod(pt[6]);  R[1,2] = _stod(pt[7]);  T[1] = _stod(pt[8]);
	R[2,0] = _stod(pt[9]); R[2,1] = _stod(pt[10]); R[2,2] = _stod(pt[11]); T[2] = _stod(pt[12]);
	s3d->setMatTran(R, T);
}
IGESReader3D::IGESReader3D(std::string &path) {
	//
	std::cout << "start parsing iges input" << std::endl;
	// init variable
	std::string line;
	std::ifstream myfile = std::ifstream(path);
	std::vector<std::vector<std::string>> d_text;
	std::vector<std::string>			  p_text;
	if (myfile.is_open()) {
		while (getline(myfile, line)) {
			switch (line[72]) {
			case 'D':
				readDPart(line,d_text); break;
			case 'P':
				readPPart(line,p_text); break;
			case 'S':
			case 'G':
			case 'T':
				break;
			default:
				std::cout << "f not right" << std::endl;
			}
		}			
		myfile.close();
	}
	else {
		std::cout << "unable to open file" << std::endl;
		std::cout << "file path: " + path << std::endl;
		assert(0);
	}
	parseAllShapes(d_text,p_text,_solidsGroup);
	//if (_surfs2D.getSurfSize() > 0) { // surface mode
	//	std::cout << "Input is surfaces" << std::endl;
	//	_surfs2D.removeConincideLines();
	//	_surfs2D.removeConincideArcs();
	//	_surfs2D.removeConcincideNURBSCrrves();
	//}
	//std::cout << "Reading all geometries successfully" << std::endl;
}
NS_END