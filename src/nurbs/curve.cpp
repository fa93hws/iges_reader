#include "curve.h"
#include "../shared/utility.h"


NSN_BEG
void Curve::curveDerivsAlg1(const std::vector<NS::Point3D> &P, const double u, const int d, std::vector<NS::Vector3D> &SKL) const {
	// not for rational
	int p = basis.order;
	int du = d > p ? p : d;
	SKL.resize(d + 1);
	int uspan = basis.findSpan(u);
	std::vector<std::vector<double>> Nu;
	basis.dersBasisFuns(u, uspan, du, Nu);
	for (int k = 0; k <= du; k++) {
		SKL[k].set(0.0, 0.0, 0.0);
		for (int r = 0; r <= p; r++)
			SKL[k] = SKL[k] + P[uspan - p + r]* Nu[k][r];
	}
}
int Curve::calCombiation(const int n, const int k)const {
	if (n == k || k == 0) return 1;
	if (k == 1) return n;
	assert(n > k);
	int denu = k;
	if (n - k < k) denu = n - k;
	int out = 1;
	for (int i = 0; i < denu; i++) out *= (n - i);
	for (int i = 2; i <= denu; i++) out /= i;
	return out;
}
void Curve::getComp(const int i, std::vector<double> &comp) const {
	if (rational) comp.push_back(weightVector[i]);
	else comp.push_back(1.0);
	comp.push_back(controlPoints[i].x() * comp[0]);
	comp.push_back(controlPoints[i].y() * comp[0]);
	comp.push_back(controlPoints[i].z() * comp[0]);
}
void Curve::calculateMRep(const int v) {
	int order = getOrder();
	if (order != controlPoints.size() - 1) return;
	assert(order == controlPoints.size() - 1);
	// cal S
	arma::mat  smat(order + v + 1, 4 * v + 4);
	smat.zeros();
	// finish cal comb,cal S
	for (int i = 0; i <= order; i++) {
		std::vector<double> comp;
		getComp(i, comp);
		for (int j = 0; j <= v; j++) {
			double comb = calCombiation(v, j);
			comb *= calCombiation(order, i);
			comb /= calCombiation(order + v, i + j);
			for (size_t k = 0; k < comp.size() ; k++)
				smat(i + j, k * (1 + v) + j) = comb * comp[k];
		}
	}
	// debug
	//arma::mat m(12,6);
	//m << 0 << 0 << 0 << 0 << 0 << 0 << arma::endr
	//	<< -1 << 1 << -1 << 1.5 << 0.5 << -1 << arma::endr
	//	<< 0 << -1 << 0 << -1 << -1 << 0 << arma::endr //
	//	<< 2 << -2 << 2 << -3 << -1 << 2 << arma::endr
	//	<< 0 << 0 << 0 << 0 << 0 << 1 << arma::endr
	//	<< 0 << 0 << 0 << 0 << 1 << 0 << arma::endr //
	//	<< 2 << -1 << 2 << -2 << 0 << 0  << arma::endr
	//	<< 1 << 0 << 0 << 0 << 0 << 0 << arma::endr
	//	<< 0 << 1 << 0 << 0 << 0 << 0 <<  arma::endr//
	//	<< 0 << 0 << 2 << -1 << 0 << 0 << arma::endr
	//	<< 0 << 0 << 1 << 0 << 0 << 0 <<  arma::endr
	//	<< 0 << 0 << 0 << 1 << 0 << 0 <<  arma::endr;
	//m.print();
	//arma::mat mult = smat * m;
	//mult.print();
	//// end debug
	// finish cal S, cal M
	//arma::mat U, V;
	//arma::vec S;
	//svd(U, S, V, smat);
	//U.print("U:");
	//S.print("S:");
	//V.print("V:");
	arma::mat m = null(smat,1e-6);
	_mRep.resize(4);
	for (int i = 0; i < 4; i++) {
		_mRep[i] = m.rows((1+v) * i, (1+v) * i + v);
		_mRep[i].debug_set_size();
	}
}
void Curve::ratCurveDerives(const int d, const double u, std::vector<NS::Vector3D> &SKL) const {
	//getPw();
	std::vector<NS::Point3D> Aw;
	for (size_t i = 0; i < controlPoints.size(); i++) 
		Aw.push_back(controlPoints[i] * weightVector[i]);
	std::vector<NS::Vector3D> aders;
	curveDerivsAlg1(Aw, u, d, aders);

	std::vector<NS::Vector3D> wders;
	std::vector<NS::Point3D> wp;
	for (size_t i = 0; i < weightVector.size(); i++)
		wp.push_back(NS::Point3D(weightVector[i], 0.0, 0.0));
	curveDerivsAlg1(wp, u, d, wders);
	//cal Bin
	PolynomialsFunctor pFunctor;
	std::vector<std::vector<int>> Bin(d + 1);
	for (int i = 0; i < d + 1; i++)
		for (int j = 0; j < d + 1; j++)
			Bin[i].push_back(pFunctor.binoCoe(i, j));
	// init skl
	if (SKL.size() == 0) {
		SKL.resize(d + 1);
		for (int i = 0; i < d; i++) 
			SKL[i].set(0.0, 0.0, 0.0);
	}

	NS::Point3D P3dTemp;
	getCoord(u, P3dTemp);
	SKL[0].set(P3dTemp.x(), P3dTemp.y(), P3dTemp.z());
	for (int k = 0; k <= d; k++) {
		NS::Vector3D v = aders[k];
		for (int i = 1; i <= k; i++)
			v = v - Bin[k][i] * wders[i].x() * SKL[k - i];
		SKL[k] = v / wders[0].x();
	}
}
void Curve::calMrepVal(const NS::Point3D &pt, arma::mat &Mp) const {
	Mp = _mRep[0] + _mRep[1] * pt.x() + _mRep[2] * pt.y() + _mRep[3] * pt.z();
}
void Curve::calMrepVal(const NS::Vector3D &pt, arma::mat &Mp) const {
	Mp = _mRep[0] + _mRep[1] * pt.x() + _mRep[2] * pt.y() + _mRep[3] * pt.z();
}

void Curve::getPw(const std::vector<NS::Point3D> &cp, std::vector<Point4f> &pw) const{
	pw.clear();
	if (!rational) {
		for (size_t i = 0; i < cp.size(); i++)
			pw.push_back(Point4f(cp[i], 1.0));
		return;
	}
	for (size_t i = 0; i < cp.size(); i++)
		pw.push_back(Point4f(cp[i], weightVector[i]));
}

void Curve::getCpFromPw(const std::vector<Point4f> &pw, std::vector<NS::Point3D> &cp, std::vector<double> &w) const {
	cp.clear();
	w.clear();
	for (Point4f c : pw) {
		w.push_back(c.w);
		c /= c.w;
		cp.push_back(NS::Point3D(c.x, c.y, c.z));
	}
}
void Curve::getCoord(const double knot, NS::Point3D &pt) const {
	pt.set(0.0, 0.0, 0.0);
	//auto basisFunsInfo = basis.basisFuns(kont);
	int i0, i1;
	std::vector<double> N;
	basis.basisFuns(knot, N, i0, i1);

	double denu = 0.0;
	for (int i = i0; i <= i1; i++) {
		NS::Point3D temp(0.0, 0.0, 0.0);
		if (rational) {
			temp = controlPoints[i] * N[i-i0] * weightVector[i];
			denu += N[i-i0] * weightVector[i];
		}
		else
			temp = controlPoints[i] * N[i-i0];
		pt += temp;
	}
	if (rational)
		pt /= denu;
}
void Curve::_getScatters(const int numPoints, CurveInfo &c) {
	//init knot vector;
	std::vector<double> k(numPoints);
	for (int i = 0; i < numPoints; i++)
		k[i] = -1.0 + i * 2.0 / (numPoints - 1);

	std::vector<std::vector<NS::Point3D>> pts = {};
	pts.push_back({});
	for (int i = 0; i < numPoints; i++) {
		NS::Point3D temp;
		getCoord(k[i],temp);
		pts[0].push_back(temp);
	}
	pts.push_back(controlPoints);
	c.pnts = pts;
	c.knots = k;
}
Curve::Curve(Basis &b, std::vector<NS::Point3D> &_cP) {
	assert(_cP.size() == b.nBasis);
	basis = b; controlPoints = _cP;
	weightVector = std::vector<double>(controlPoints.size(), 1.0);
	int ord = getOrder();
	calculateMRep(ord > 1 ? ord - 1 : 1);
	//calculateMRep(1);
}
Curve::Curve(Basis &b, std::vector<NS::Point3D> &_cP, std::vector<double> &w) {
	assert(_cP.size() == b.nBasis);
	assert(_cP.size() == w.size());
	basis = b; controlPoints = _cP; weightVector = w;
	rational = true;
	int ord = getOrder();
	calculateMRep(ord > 1 ? ord - 1 : 1);
}
Curve::Curve(const Curve &c) { 
	basis = c.basis; 
	controlPoints = c.controlPoints;
	weightVector = c.weightVector; 
	rational = c.rational; 
	c.getMRep(_mRep);
}
double Curve::arc_length(const double k0, const double k1) const{
	assert(k1 > k0);
	//gauss quadrature
	double k[6] = { -0.932469514203152, - 0.661209386466264, - 0.238619186083197, 0.238619186083197, 0.661209386466264, 0.932469514203152 };
	double w[6] = { 0.171324492379171, 0.360761573048139, 0.467913934572691, 0.467913934572691, 0.360761573048139, 0.171324492379171 };

	double scale  = (k1 - k0) / 2;
	double offset = k0 + scale;
	double len=0.0;
	for (int i = 0; i < 6; i++) {
		std::vector<NS::Vector3D> temp;
		ratCurveDerives(1, k[i]*scale+offset, temp);
		len += temp[1].length() * w[i];
	}
	return len*scale;
}
void Curve::reverse() {
	std::reverse(weightVector.begin(), weightVector.end());
	std::reverse(controlPoints.begin(), controlPoints.end());
	basis.reverse();
}
//advance algorithm
void Curve::curveKnotIns(double u, int r) {
	//insert knot k r times
	std::vector<Point4f> pw;
	getPw(controlPoints,pw);
	size_t np = controlPoints.size(); // number of control points
	int p = basis.order; // order of the basis
	int mp = np + p + 1;
	int nq = np + r;	// number of control points after knot insertion;
	std::vector<double> U = basis.knotVector;
	// get multiplicy of k
	int s = 0;
	for (double temp : basis.knotVector) {
		if (abs(u - temp) < absTOL)
			s++;
		if (u + 0.1 < temp)
			break;
	}
	// insert knot value into knot vector
	int k = basis.insert(u, r);
	// save unchanged control points
	std::vector<Point4f> _pw(nq);
	std::vector<Point4f> rw(p - s + 1);
	for (int i = 0; i <= k - p; i++) _pw[i] = pw[i];
	for (size_t i = k - s; i < np; i++) _pw[i + r] = pw[i];
	for (int i = 0; i <= p - s; i++) rw[i] = pw[k - p + i];
	int L;
	for (int j = 1; j <= r; j++) {
		L = k - p + j;
		for (int i = 0; i <= p - j - s; i++) {
			double alpha = (u - U[L + i]) / (U[i + k + 1] - U[L + i]);
			rw[i] = rw[i + 1] * alpha + rw[i] * (1.0 - alpha);
		}
		_pw[L] = rw[0];
		_pw[k + r - j - s] = rw[p - j - s];
	}
	for (int i = L + 1; i < k - s; i++)
		_pw[i] = rw[i - L];
	std::vector<NS::Point3D> pts;
	std::vector<double> w;
	getCpFromPw(_pw,pts,w);
	controlPoints = pts;
	weightVector = w;
}
double Curve::calculatePerImage(const NS::Point3D &pt) const {
	arma::mat M;
	calMrepVal(pt, M);
	//double dist = arma::det(M * M.t());
	M = M.t();
	
	arma::mat U, _V;
	arma::vec S;
	arma::svd(U, S, _V, M);
	_V.print("V:");
	arma::vec V= null(M, 1e-6);
	double u = V[1];
	u /= (V.n_rows - 1)*V[0] + V[1];
	return u * 2 - 1; // to [-1,1]
}
void Curve::calculateLineIntecs(const NS::Point3D &ptOnLine, const NS::Vector3D &slope, 
	std::vector<double> &intecUs, std::vector<NS::Point3D> &intecPts) const {
	_mRep[0].print("m0:");
	_mRep[1].print("m1:");
	_mRep[2].print("m2:");
	_mRep[3].print("m3:");
	bool debug = true;
	double TOL = 1e-6;
	arma::mat A, B;
	calMrepVal(ptOnLine, A);
	calMrepVal(slope, B); B = B - _mRep[0]; B = -B;
	// test

	//arma::mat M0, M1, M2, M3;
	//M0 << 1 << 0 << 0 << 0 << 0 << 0 << arma::endr
	//	<< -3 << 1 << 0 << 0 << 0 << 0 << arma::endr
	//	<< 1 << -3 << 0 << 0 << 0 << 0 << arma::endr
	//	<< 0 << 1 << 0 << 0 << 0 << 0 << arma::endr;
	//M1 << 1 << 0 << 3 << 0 << 0 << 0 << arma::endr
	//	<< 0 << 1 << -1 << 3 << 0 << 0 << arma::endr
	//	<< 0 << 0 << 1 << -1 << 0 << 0 << arma::endr
	//	<< 0 << 0 << 0 << 1 << 0 << 0 << arma::endr;
	//M2 << 0 << 0 << -3 << 0 << 2 << 0 << arma::endr
	//	<< 0 << 0 << -3 << -3 << 0 << 2 << arma::endr
	//	<< 0 << 0 << 3 << -3 << 0 << 0 << arma::endr
	//	<< 0 << 0 << 0 << 3 << 0 << 0 << arma::endr;
	//M3 << 0 << 0 << 0 << 0 << 2 << 0 << arma::endr
	//	<< 0 << 0 << 0 << 0 << -2 << -2 << arma::endr
	//	<< 0 << 0 << 0 << 0 << 1 << -2 << arma::endr
	//	<< 0 << 0 << 0 << 0 << 0 << 1 << arma::endr;
	//M0 << 0 << 0 << 0 << 0 << 0 << 0 << -1 << arma::endr
	//	<< 0 << 0 << 0 << 1 << 0 << 0 << 0 << arma::endr
	//	<< 1 << 0 << 0 << 0 << 0 << 0 << 0 << arma::endr
	//	<< 0 << 0 << 0 << 0 << 1 << 0 << 0 << arma::endr
	//	<< 0 << 1 << 0 << 0 << 0 << 1 << -1 << arma::endr
	//	<< 0 << 0 << 1 << 0 << 0 << 0 << 0 << arma::endr;
	//M1 << 0 << 0 << 0 << -1 << 0 << 0 << 0 << arma::endr
	//	<< 0 << 0 << 0 << 0 << -1 << 0 << 0 << arma::endr
	//	<< -1 << 0 << 0 << 0 << 0 << -1 << 0 << arma::endr
	//	<< 0 << 0 << 0 << 0 << 0 << 0 << -1 << arma::endr
	//	<< 0 << -1 << 0 << 0 << 0 << 0 << 0 << arma::endr
	//	<< 0 << 0 << -1 << 0 << 0 << 0 << 0 << arma::endr;
	//M2 << 0 << 0 << 0 << 0 << 0 << 0 << 1 << arma::endr
	//	<< 0 << 0 << 0 << 0 << 0 << 0 << 0 << arma::endr
	//	<< -1 << 0 << 0 << -1 << 0 << 0 << 0 << arma::endr
	//	<< 0 << 0 << 0 << 0 << 0 << 0 << 0 << arma::endr
	//	<< 0 << -1 << 0 << 0 << -1 << 0 << 0 << arma::endr
	//	<< 0 << 0 << -1 << 0 << 0 << -1 << 0 << arma::endr;
	//M3 << 0 << 0 << 0 << 1 << 0 << 0 << 0 << arma::endr
	//	<< 1 << 0 << 0 << 0 << 1 << 0 << 0 << arma::endr
	//	<< 0 << 0 << 0 << 0 << 0 << 1 << 0 << arma::endr
	//	<< 0 << 1 << 0 << 0 << 0 << 0 << 0 << arma::endr
	//	<< 0 << 0 << 1 << 0 << 0 << 0 << 0 << arma::endr
	//	<< 0 << 0 << 0 << 0 << 0 << 0 << 0 << arma::endr;
	//int nrow = M0.n_rows;
	//int ncol = M0.n_cols;
	//A = arma::mat(nrow * 3 + ncol, nrow * 4); A.zeros();
	//B = arma::mat(nrow * 3 + ncol, nrow * 4); B.zeros();
	//A.submat(arma::span(0, nrow - 1), arma::span(nrow, 2 * nrow - 1)) = arma::eye(nrow, nrow);
	//A.submat(arma::span(nrow, 2 * nrow - 1), arma::span(2 * nrow, 3 * nrow - 1)) = arma::eye(nrow, nrow);
	//A.submat(arma::span(nrow * 3, nrow * 3 + ncol - 1), arma::span(0, nrow - 1)) = M0.t();
	//A.submat(arma::span(nrow * 3, nrow * 3 + ncol - 1), arma::span(nrow, 2*nrow - 1)) = M1.t();
	//A.submat(arma::span(nrow * 3, nrow * 3 + ncol - 1), arma::span(2*nrow, 3*nrow - 1)) = M2.t();
	//A = A.t();
	//B.submat(arma::span(0, nrow-1), arma::span(0, nrow - 1)) = arma::eye(nrow, nrow);
	//B.submat(arma::span(nrow, 2 * nrow - 1), arma::span(nrow, 2 * nrow - 1)) = arma::eye(nrow, nrow);
	//B.submat(arma::span(nrow * 3, nrow * 3 + ncol - 1), arma::span(2 * nrow, 3 * nrow - 1)) = -M3.t();
	//B = B.t();

	//A = arma::mat(3, 4); A.zeros(); A(2, 1) = 1; A(1, 2) = -1; A(0, 3) = 1;
	//B = arma::mat(3, 4); B.zeros(); B(2, 1) = 1; B(1, 2) = 1; B(0, 3) = 1;
	//A = arma::mat(4, 4); A.zeros(); A(0, 0) = -1; A(1, 0) = -1; A(2, 0) = -1; A(0, 1) = 1; A(2, 1) = 2; A(3, 2) = -1; A(2, 3) = -1;
	//B = arma::mat(4, 4); B.zeros(); B(0, 1) = 1; B(0, 3) = 1; B(1, 0) = 1; B(1, 1) = -2; B(1, 2) = 1; B(2, 1) = 1; B(3, 0) = 1; B(3, 1) = -2;
	//A = arma::mat(3, 4); A.zeros(); A(2, 0) = -1; A(0, 1) = 1; A(0, 2) = -1; A(1, 2) = -1;
	//B = arma::mat(3, 4); B.zeros(); B(0, 0) = 1; B(1, 0) = 1; B(2, 0) = -1; B(0, 1) = 2; B(2, 1) = -2; B(0, 2) = -1; B(2, 2) = 2; B(1, 3) = -1; B(2, 3) = 1;
	//B = arma::mat(3, 4);
	//B << -1 << 0 << 1 << 1 << arma::endr
	//	<< 0 << -1 << 1 << -1 << arma::endr
	//	<< 1 << 1 << 1 << 0 << arma::endr;
	//A = arma::mat(3, 4); A.zeros();  A(2, 1) = 1; A(1, 2) = -1; A(0, 3) = 1;
	//B = -B;

	// test
	if (debug) {
		A.print("A:");
		B.print("B:");
	}
	int ro,ra12p;
	arma::mat U, V, A1, A11,A12,A12p, A1p, B1, B1p;
	arma::vec S;
	int times = 0;
	do {
		times++;
		ro = arma::rank(B, 1e-6);
		if (ro == B.n_cols) {
			if (B.n_rows == ro) break;
			else {
				A = A.t();
				B = B.t();
			}
		}
		arma::svd(U, S, V, B);
		B1 = B*V;
		if (debug) B1.print("B1:");
		A1 = A*V;
		if (debug) A1.print("A1:");
		A12 = A1.cols(ro, A1.n_cols - 1);
		if (debug) A12.print("A12:");
		arma::svd(U, S, V, A12);
		A12p = U.t() * A12;
		if (debug) A12p.print("A12p:");
		ra12p = arma::rank(A12p,1e-6);
		A1p = U.t()*A1;
		B1p = U.t()*B1;
		if (debug) { A1p.print("A1p:"); B1p.print("B1p:"); }
		ro = arma::rank(B1p);
		A = A1p.submat(arma::span(ra12p, A1p.n_rows - 1), arma::span(0, ro - 1));
		B = B1p.submat(arma::span(ra12p, A1p.n_rows - 1), arma::span(0, ro - 1));
		if (debug) A.print("A:");
		if (debug) B.print("B:");
		if (times > 10) { std::cout << "error in matrixRep find line,curve intersection"; assert(false); }
	} while (true);
	if (debug) A.print("A:");
	if (debug) B.print("B:");
	arma::cx_mat eig_vec;
	arma::cx_vec eig_val;
	arma::mat AB = A*B.i();
	arma::eig_gen(eig_val,eig_vec, AB);
	eig_val.print("eigval:");
	arma::vec eig = arma::real(eig_val);
	NS::Point3D intec(ptOnLine + eig[0] * slope);
	arma::mat M;
	calMrepVal(intec, M);
	double dist = arma::det(M*M.t());
	//double u = calculatePerImage(intec);
	auto s = 1;
	
}
//Export Curve2D
void Curve::getScattors(const int numPoints, CurveInfo &c) { _getScatters(numPoints + 1,c); }
void Curve::getScattors(CurveInfo &c) { _getScatters(101,c); }

/************************************ Parent Curve ***************************************************/
ParentCurve::ParentCurve(const NURBS::Curve _c) {
	basis = _c.basis; controlPoints = _c.controlPoints; weightVector = _c.weightVector;
	rational = _c.rational;
}
ParentCurve::ParentCurve(const NURBS::Basis &b, const std::vector<NS::Point3D> &_cP) {
	assert(controlPoints.size() != basis.nBasis);
	basis = b; controlPoints = _cP;
}
ParentCurve::ParentCurve(const NURBS::Basis &b, const std::vector<NS::Point3D> &_cP, const std::vector<double> &w) {
	assert(controlPoints.size() != basis.nBasis);
	assert(controlPoints.size() != basis.nBasis);
	basis = b; controlPoints = _cP; weightVector = w;
	rational = true;
}
void ParentCurve::subDivision() {
	children.clear();
	int npt = basis.order + 1;
	int nchild = basis.nBasis - npt;
	if (nchild == 0 && rational) // if there is no interior knot
		children.push_back(Curve(basis, controlPoints, weightVector));
	else if (nchild == 0 && !rational)
		children.push_back(Curve(basis, controlPoints));
	// find interior knot and its multiplicy
	if (nchild == 0) return;
	std::vector<double> iKnot(nchild);
	std::vector<int> s(nchild, 1);
	iKnot[0] = basis.knotVector[npt];
	int idx = 0;
	for (size_t i = npt + 1; i < basis.knotVector.size() - npt; i++) {
		if (abs(iKnot[idx] - basis.knotVector[i]) < absTOL)
			s[idx] ++;
		else {
			idx++;
			iKnot[idx] = basis.knotVector[i];
		}
	}
	Curve c2D(*this);
	// insert knots
	for (int i = 0; i < idx + 1; i++)
		if (basis.order - s[i] < 1)
			continue;
		else
			c2D.curveKnotIns(iKnot[i], basis.order - s[i]);
	// construct basis for subcurves
	std::vector<double> V;
	for (int i = 0; i <= basis.order; i++) V.push_back(-1.0);
	for (int i = 0; i <= basis.order; i++) V.push_back(1.0);
	Basis b(V);
	// Extract subcurve
	for (int i = 0; i < idx + 2; i++) {
		int begin_idx = i * (npt - 1);
		int end_idx = (i + 1) * (npt - 1);
		std::vector<NS::Point3D> cp(c2D.controlPoints.begin() + begin_idx, c2D.controlPoints.begin() + end_idx + 1);
		std::vector<double> wei(c2D.weightVector.begin() + begin_idx, c2D.weightVector.begin() + end_idx + 1);
		if (rational)
			children.push_back(Curve(b, cp, wei));
		else
			children.push_back(Curve(b, cp));
	}
}

NS_END
