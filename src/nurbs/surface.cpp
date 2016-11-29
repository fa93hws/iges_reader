#include "surface.h"


NSN_BEG
int Surface::calCombiation(const int n, const int k)const {
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
void Surface::getComp(const int i, const int j, std::vector<double> &comp) const {
	if (rational) comp.push_back(weightVector[i][j]);
	else comp.push_back(1.0);
	comp.push_back(controlPoints[i][j].x() * comp[0]);
	comp.push_back(controlPoints[i][j].y() * comp[0]);
	comp.push_back(controlPoints[i][j].z() * comp[0]);
}
void Surface::calMrepVal(const NS::Point3D &pt, arma::mat &Mp) const {
	Mp = _mRep[0] + _mRep[1] * pt.x() + _mRep[2] * pt.y() + _mRep[3] * pt.z();
}
void Surface::calMrepVal(const NS::Vector3D &pt, arma::mat &Mp) const {
	Mp = _mRep[0] + _mRep[1] * pt.x() + _mRep[2] * pt.y() + _mRep[3] * pt.z();
}
void Surface::calculateMRep() {
	int v1, v2;
	// cal S
	int d1 = basis[0].order;
	int d2 = basis[1].order;
	if (d1 > d2) {
		v1 = 2 * d2 - 1;
		v2 = d1 - 1;
	}
	else {
		v1 = d2 - 1;
		v2 = 2 * d1 - 1;
	}
	arma::mat  smat((v1 + d1 + 1)*(v2 + d2 + 1), 4 * (v1 + 1)*(v2 + 1));
	smat.zeros();
	// finish cal comb,cal S
	for (int i = 0; i <= d1; i++) {
		for (int j = 0; j <= d2; j++) {
			std::vector<double> comp;
			getComp(i,j, comp);
			for (int k = 0; k <= v1; k++) {
				for (int l = 0; l <= v2; l++) {
					double comb = calCombiation(v1, k);
					comb *= calCombiation(v2, l);
					comb *= calCombiation(d1, i);
					comb *= calCombiation(d2, j);
					//comb /= calCombiation(v1 + v2, i + k);
					//comb /= calCombiation(d1 + d2, j + l);
					comb /= calCombiation(v1 + d1, i + k);
					comb /= calCombiation(v2 + d2, j + l);
					//int row = (j+l)*(v1 + d1 + 1) + k + i;
					int row = (k + i)*(v2 + d2 + 1) + j + l;
					for (int n = 0; n < 4; n++) {		
						int col = n * (1 + v1)*(v2 + 1) + k * (v2 + 1) + l;		
						//int col = n * (1 + v1)*(v2 + 1) + l * (v1 + 1) + k;
						assert(smat(row, col) == 0);
						smat(row, col) = comb * comp[n];
					}
				}
			}
		}
	}
	// //verify
	//double u = -0.5;
	//double v = 0.25;
	//arma::mat L(1, (v1 + d1 + 1)*(v2 + d2 + 1));
	//arma::mat R(1, 4*(v1 + 1)*(v2 + 1));
	//std::vector<double> k = { -1,-1,-1,1,1,1 };
	//std::vector<Basis> b(2);
	//b[0] = Basis(k);
	//k = { -1,-1,-1,-1,1,1,1,1 };
	//b[1] = Basis(k);
	//for (int i = 0; i <= v1 + d1; i++) {
	//	for (int j = 0; j <= v2 + d2; j++) {
	//		std::vector<double> N;
	//		int pos,pos1;
	//		b[0].basisFuns(u, N, pos, pos1);
	//		double B0 = N[i];
	//		b[1].basisFuns(v, N, pos, pos1);
	//		double B1 = N[j];
	//		L(0, i*(v2 + d2 + 1) + j) = B0*B1;
	//	}
	//}
	//k = { -1,-1,1,1 };
	//b[0] = Basis(k);
	//b[1] = Basis(k);
	//for (int i = 0; i <= v1; i++) {
	//	for (int j = 0; j <= v2; j++) {
	//		std::vector<double> val(4,0);
	//		for (int k = 0; k <= basis[0].order; k++) {
	//			std::vector<double> N1, N2;
	//			int pos;
	//			basis[0].basisFuns(u, N1, pos, pos);
	//			for (int l = 0; l <= basis[1].order; l++) {
	//				basis[1].basisFuns(v, N2, pos, pos);
	//				std::vector<double> comp;
	//				getComp(k, l, comp);
	//				for (int m = 0; m < 4; m++)
	//					val[m] += N1[k] * N2[l] * comp[m];
	//			}
	//		}
	//		std::vector<double> N;
	//		int pos, pos1;
	//		b[0].basisFuns(u, N, pos, pos1);
	//		double B0 = N[i];
	//		b[1].basisFuns(v, N, pos, pos1);
	//		double B1 = N[j];
	//		for (int m = 0; m < 4; m++)
	//			R(0, m * (v1 + 1)*(v2 + 1) + i*(v2 + 1) + j) =B0*B1* val[m];
	//	}
	//}
	//L.print("L:");
	//(L*smat).print("L*S:");
	//R.print("R:");
	//(L*smat - R).print("LS-R");

	// verify
	//arma::mat _m(16,4);
	//_m << -1 << 1 << 0 << 0 << arma::endr
	//	<< -1 << 0 << 0 << 0 << arma::endr
	//	<< -1 << 2 << 0 << -1 << arma::endr
	//	<< 0 << 0 << -1 << 0 << arma::endr//
	//	<< 1 << -1 << 0 << 0 << arma::endr
	//	<< 1 << -2 << 0 << 0 << arma::endr
	//	<< 1 << -2 << 0 << 0 << arma::endr
	//	<< 0 << -2 << 0 << 0 << arma::endr//
	//	<< 0 << 0 << 0 << 1 << arma::endr
	//	<< 0 << 0 << 1 << 0 << arma::endr
	//	<< 0 << 0 << 0 << 1 << arma::endr
	//	<< 0 << 0 << 1 << 0 << arma::endr//
	//	<< 0 << 1 << 0 << 0 << arma::endr
	//	<< 1 << 0 << 0 << 0 << arma::endr
	//	<< 0 << 1 << 0 << 0 << arma::endr
	//	<< 1 << 0 << 0 << 0 << arma::endr;
	//(smat*_m).print("0?");
	//smat.print("S:");
	//_m.print("_m:");
	//(R*_m).print("R*_m");
	arma::mat m = null(smat, 1e-6);
	//m.print("M:");
	//(R*m).print("R*m");
	_mRep.resize(4);
	for (int i = 0; i < 4; i++) {
		_mRep[i] = m.rows((1 + v1)*(1+v2) * i, (1 + v1)*(1+v2) * (i+1)-1 );
		_mRep[i].debug_set_size();
	}
}
void Surface::getPw() {
	if (!rational) {
		// to be added for non-rational cases
		return;
	}
	pw.clear();
	for (size_t i = 0; i < controlPoints.size(); i++) {
		std::vector<Point4f> row_temp;
		for (size_t j = 0; j < controlPoints[i].size(); j++)
			row_temp.push_back(Point4f(controlPoints[i][j], weightVector[i][j]));
		pw.push_back(row_temp);
	}
}
void Surface::getCpFromPw() {
	if (!rational || pw.empty())
		return;
	weightVector.clear();
	controlPoints.clear();
	for (size_t i = 0; i < pw.size(); i++) {
		std::vector<double> w_row;
		std::vector<NS::Point3D> cp_row;
		for (int j = 0; j < pw[i].size(); j++) {
			w_row.push_back(pw[i][j].w);
			pw[i][j] /= pw[i][j].w;
			cp_row.push_back(NS::Point3D(pw[i][j].x, pw[i][j].y, pw[i][j].z));
		}
		weightVector.push_back(w_row);
		controlPoints.push_back(cp_row);
	}
}
void Surface::getCoord(const double knot0, const double knot1, NS::Point3D &pt) const {
	double knot[] = { knot0,knot1 };
	std::vector<std::vector<double>> N(2);
	std::vector<int> i0(2);
	std::vector<int> i1(2);
	for (int i = 0; i < 2; i++)
		basis[i].basisFuns(knot[i], N[i], i0[i], i1[i]);
	pt.set(0.0, 0.0, 0.0);

	double denu = 0.0; // only for rational
	for (int i = i0[0]; i <= i1[0]; i++) {
		for (int j = i0[1]; j <= i1[1]; j++) {
			if (rational) {
				pt += controlPoints[i][j] *
					(N[0][i - i0[0]] * N[1][j - i0[1]] * weightVector[i][j]);
				denu += N[0][i - i0[0]] * N[1][j - i0[1]] * weightVector[i][j];
			}
			else
				pt += controlPoints[i][j] * (N[0][i - i0[0]] * N[1][j - i0[1]]);
		}
	}
	if (rational)
		pt /= denu;
}
void Surface::_getGridPoints(const int numPoints, SurfaceInfo &s) {
	// Init Knot Vector
	std::vector<double> temp;
	for (int i = 0; i < numPoints; i++)
		temp.push_back(-1.0 + i * 2.0 / (numPoints - 1));
	std::vector<std::vector<double>> k = {};
	for (int i = 0; i < 2; i++)
		k.push_back(temp);
	// Cal points for each pair of knot value
	std::vector<std::vector<NS::Point3D>> gridPoints(numPoints);
	for (int i = 0; i < numPoints; i++) {
		gridPoints[i].resize(numPoints);
		for (int j = 0; j < numPoints; j++) {
			getCoord(k[0][i], k[1][j], gridPoints[i][j]);
		}
	}
	s.pnts = gridPoints;
	s.knots = k;
	//return SurfaceInfo(gridPoints, k);
}
void Surface::surfaceDerivsAlg1(const std::vector<std::vector<NS::Point3D>> &P, const double u,
	const double v, const int d, std::vector<std::vector<NS::Vector3D>> &SKL) {
	// not for rational
	int p = basis[0].order;
	int q = basis[1].order;
	int du = d > p ? p : d;
	int dv = d > q ? q : d;
	SKL.resize(d + 1);
	for (size_t i = 0; i < SKL.size(); i++) {
		SKL[i].resize(d + 1);
	}

	int uspan = basis[0].findSpan(u);
	std::vector<std::vector<double>> Nu;
	basis[0].dersBasisFuns(u, uspan, du, Nu);
	int vspan = basis[1].findSpan(v);
	std::vector<std::vector<double>> Nv;
	basis[1].dersBasisFuns(v, vspan, dv, Nv);

	std::vector<NS::Vector3D> temp(q + 1);
	for (int k = 0; k <= du; k++) {
		for (int s = 0; s <= q; s++) {
			temp[s].set(0.0, 0.0, 0.0);
			for (int r = 0; r <= p; r++)
				temp[s] = temp[s] + P[uspan - p + r][vspan - q + s] * Nu[k][r];
		}
		int dd = d - k > dv ? dv : d - k;
		for (int l = 0; l <= dd; l++) {
			SKL[k][l].set(0.0, 0.0, 0.0);
			for (int s = 0; s <= q; s++)
				SKL[k][l] = SKL[k][l] + temp[s] * Nv[l][s];
		}
	}
}

void Surface::ratSurfaceDerives(const int d, const double u, const double v,
	std::vector<std::vector<NS::Vector3D>> &SKL) {
	//getPw();
	std::vector<std::vector<NS::Point3D>> Aw;		
	for (size_t i = 0; i < controlPoints.size(); i++) {
		std::vector<NS::Point3D> Aw_row;
		for (size_t j = 0; j < controlPoints[i].size(); j++) {
			Aw_row.push_back( controlPoints[i][j] * weightVector[i][j] );
		}
		Aw.push_back(Aw_row);
	}
		
	std::vector<std::vector<NS::Vector3D>> aders;
	surfaceDerivsAlg1(Aw, u, v, d, aders);
	std::vector<std::vector<NS::Vector3D>> wders;
	std::vector<std::vector<NS::Point3D>> wp;
	for (size_t i = 0; i < weightVector.size(); i++)
	{
		std::vector<NS::Point3D> temp;
		for (size_t j = 0; j < weightVector[i].size(); j++)
			temp.push_back(NS::Point3D(weightVector[i][j], 0.0, 0.0));
		wp.push_back(temp);
	}
	surfaceDerivsAlg1(wp, u, v, d, wders);
	//cal Bin
	PolynomialsFunctor pFunctor;
	std::vector<std::vector<int>> Bin(d + 1);
	for (int i = 0; i < d + 1; i++)
		for (int j = 0; j < d + 1; j++)
			Bin[i].push_back(pFunctor.binoCoe(i, j));
	// init skl
	if (SKL.size() == 0) {
		SKL.resize(d + 1);
		for (int i = 0; i < d; i++) {
			SKL[i].resize(d + 1);
			for (int k = 0; k < d; k++)
					SKL[i][k].set(0.0,0.0,0.0);
		}
	}
	NS::Point3D P3dTemp;
	getCoord(u, v, P3dTemp);
	SKL[0][0].set(P3dTemp.x(), P3dTemp.y(), P3dTemp.z());
	//cal
	for (int k = 0; k <= d; k++)
	{
		for (int l = 0; l <= d - k; l++)
		{
			NS::Vector3D v = aders[k][l];
			for (int j = 1; j <= l; j++)
					v = v - Bin[l][j] * wders[0][j].x() * SKL[k][l - j];
			for (int i = 1; i <= k; i++){
				v = v - Bin[k][i] * wders[i][0].x() * SKL[k - i][l];
				NS::Vector3D v2(0.0,0.0,0.0);
				for (int j = 1; j <= l; j++)
						v2 = v2 + Bin[l][j] * wders[i][j].x() * SKL[k - i][l - j];
					v= v - Bin[k][i] * v2;
			}
			SKL[k][l]= v / wders[0][0].x();
		}
	}
}


Surface::Surface(std::vector<Basis> b, std::vector<std::vector<NS::Point3D>> cP) {
	assert(b.size() != 2);
	assert(b[1].nBasis != cP.size());
	for (std::vector<NS::Point3D> c : cP)
		assert(b[0].nBasis != c.size());
	basis = b; controlPoints = cP;
	//calculateMRep();
	//getConvexHull();
}
Surface::Surface(std::vector<Basis> b, std::vector<std::vector<NS::Point3D>> cP,
	std::vector<std::vector<double>> w) {
	assert(b.size() == 2);
	assert(b[0].nBasis == cP.size());
	assert(cP.size() == w.size());
	for (std::vector<NS::Point3D> c : cP)
		assert(b[1].nBasis == c.size());
	for (std::vector<double> _w : w)
		assert(b[1].nBasis == _w.size());
	basis = b; controlPoints = cP; weightVector = w; rational = true;
	calculateMRep();
	//getConvexHull();
}
void Surface::projPt_initialGuessByLS(const NS::Point3D &pt, double &u0, double &v0) {
	if (referPoints.size() == 0) {
		std::vector<double> u(numReferPoints);
		for (int i = 0; i < numReferPoints; i++) {
			u[i] = -1.0 + i * 2.0 / (numReferPoints - 1);
		}
		referPoints.resize(numReferPoints);
		for (int i = 0; i < numReferPoints;i++)
			referPoints[i].resize(numReferPoints);
		for (int i = 0; i < numReferPoints; i++)
			for (int j = 0; j < numReferPoints; j++) 
				getCoord(u[i], u[j], referPoints[i][j]);
	}
	// find nearest referencePoints
	int idx0, idx1;
	double dist0=99999.0;
	for (int i = 0; i < numReferPoints ; i++)
		for (int j = 0; j < numReferPoints; j++) {
			double dist = pt.distanceTo(referPoints[i][j]);
			if (dist < dist0) {
				dist0 = dist; idx0 = i; idx1 = j;
			}
		}
	u0 = -1.0 + idx0  * 2.0 / (numReferPoints - 1);
	v0 = -1.0 + idx1  * 2.0 / (numReferPoints - 1);
	return;


	//size_t ncp_row = controlPoints.size();
	//size_t npc_col = controlPoints[0].size();
	//// less than four points (point share the same coordinate)
	//std::vector<NS::Point3D> kpts;
	//kpts.push_back(controlPoints[0][0]);
	//kpts.push_back(controlPoints[ncp_row - 1][0]);
	//kpts.push_back(controlPoints[ncp_row - 1][npc_col - 1]);
	//kpts.push_back(controlPoints[0][npc_col - 1]);
	//for (int i = 0; i < 3; i++) {
	//	for (int j = i + 1; j < 4;j++)
	//		if (kpts[i].distanceTo(kpts[j]) < absTOL) {
	//			u0 = 0.0; v0 = 0.0; return;
	//		}
	//}
	//// four points
	//std::vector<std::vector<double>> _A = { { controlPoints[0][0].x(),controlPoints[0][0].y() ,controlPoints[0][0].z() ,1.0 }
	//	,{ controlPoints[ncp_row - 1][0].x(),controlPoints[ncp_row - 1][0].y() ,controlPoints[ncp_row - 1][0].z() ,1.0 }
	//	,{ controlPoints[ncp_row - 1][npc_col - 1].x(),controlPoints[ncp_row - 1][npc_col - 1].y() ,controlPoints[ncp_row - 1][npc_col - 1].z() ,1.0 }
	//,{ controlPoints[0][npc_col - 1].x(),controlPoints[0][npc_col - 1].y() ,controlPoints[0][npc_col - 1].z() ,1.0 } };
	//Matrix A(_A);
	//double det = A.calDet();
	//if (abs(det) < absTOL) {
	//	u0 = 0.0; v0 = 0.0; return;
	//}
	//A = A.transpose();

	//std::vector<std::vector<double>> _H1 = { { -1 },{ 1 },{ 1 },{ -1 } };//u
	//Matrix H1(_H1);
	//std::vector<std::vector<double>> _H2 = { { -1 },{ -1 },{ 1 },{ 1 } };//v
	//Matrix H2(_H2);

	//Matrix coeu, coev;
	//LinearAlgebra::leastSquare(A, H1, coeu);
	//LinearAlgebra::leastSquare(A, H2, coev);

	//u0 = pt.x() * coeu.values[0][0] + pt.y() * coeu.values[1][0] + pt.z() * coeu.values[2][0] + coeu.values[3][0];
	//v0 = pt.x() * coev.values[0][0] + pt.y() * coev.values[1][0] + pt.z() * coev.values[2][0] + coev.values[3][0];
	//if (u0 < -1.0) u0 = -1.0; if (u0 > 1.0) u0 = 1.0;
	//if (v0 < -1.0) v0 = -1.0; if (v0 > 1.0) v0 = 1.0;

	//Point3D initPt = getCoord(u0, v0);
}

void Surface::getScattors(const int numPoints, SurfaceInfo &s) { return _getGridPoints(numPoints + 1, s); }
void Surface::getScattors(SurfaceInfo &s) { return _getGridPoints(101, s); }

void Surface::surfaceKnotIns(const int dir, const double u, const int r) {
	getPw();
	//insert knot u into dir knotvector r times.
	std::vector<double> U = basis[dir].knotVector;
	//vector<double> V = basis[1].knotVector;
	int p = basis[dir].order;
	//int q = basis[1].order;
	int np = basis[dir].nBasis;
	int nq = basis[1 - dir].nBasis;
	int s = 0;
	for (double temp : U) {
		if (abs(u - temp) < absTOL)
			s++;
		if (u + 0.1 < temp)
			break;
	}
	int k = basis[dir].insert(u, r);
	//Init alpha
	double **alpha = new double*[p - s];
	for (int i = 0; i < p - s; i++)
		alpha[i] = new double[r + 1];
	std::vector<Point4f> rw(p - s + 1);
	//Calculate alpha for later use
	int L;
	for (int j = 1; j <= r; j++) {
		L = k - p + j;
		for (int i = 0; i <= p - j - s; i++)
			alpha[i][j] = (u - U[L + i]) / (U[i + k + 1] - U[L + i]);
	}
	if (dir == 0) {// insert a pt in each column
					//Init output pw
		std::vector<std::vector<Point4f>> _pw(np + r);
		for (std::vector<Point4f> &v : _pw) for (int i = 0; i < nq; i++) v.push_back(Point4f(0.0, 0.0, 0.0, 0.0));
		// calculate for each iteration
		for (int row = 0; row < nq; row++) {
			// front unchanged control points
			for (int i = 0; i <= k - p; i++)	_pw[i][row] = pw[i][row];
			// tail unchanged control points
			for (int i = k - s; i < np; i++)	_pw[i + r][row] = pw[i][row];
			// temp
			for (int i = 0; i <= p - s; i++)	rw[i] = pw[k - p + i][row];
			// insert r times
			for (int j = 1; j <= r; j++) {
				L = k - p + j;
				for (int i = 0; i <= p - j - s; i++)
					rw[i] = rw[i + 1] * alpha[i][j] + rw[i] * (1.0 - alpha[i][j]);
				_pw[L][row] = rw[0];
				_pw[k + r - j - s][row] = rw[p - j - s];
			}
			for (int i = L + 1; i < k - s; i++)	_pw[i][row] = rw[i - L];
		}
		pw = _pw;
	}

	else {// insert a pt in each row
			//Init output pw
		std::vector<std::vector<Point4f>> _pw(nq);
		for (std::vector<Point4f> &v : _pw) for (int i = 0; i < np + r; i++) v.push_back(Point4f(0.0, 0.0, 0.0, 0.0));
		// calculate for each iteration
		for (int col = 0; col < nq; col++) {
			// front unchanged control points
			for (int i = 0; i <= k - p; i++)	_pw[col][i] = pw[col][i];
			// tail unchanged control points
			for (int i = k - s; i < np; i++)	_pw[col][i + r] = pw[col][i];
			// temp
			for (int i = 0; i <= p - s; i++)	rw[i] = pw[col][k - p + i];
			// insert r times
			for (int j = 1; j <= r; j++) {
				L = k - p + j;
				for (int i = 0; i <= p - j - s; i++)
					rw[i] = rw[i + 1] * alpha[i][j] + rw[i] * (1.0 - alpha[i][j]);
				_pw[col][L] = rw[0];
				_pw[col][k + r - j - s] = rw[p - j - s];
			}
			for (int i = L + 1; i < k - s; i++)	_pw[col][i] = rw[i - L];
		}
		pw = _pw;
	}


	getCpFromPw();
}
void Surface::getConvexHull() {
	ConvexHull polyhedron(controlPoints);
	if (polyhedron.faces.size() == 0) {
		is_plane = true;
		NS::Point3D  pt0= controlPoints[0][0];
		NS::Point3D  pt1 = controlPoints[1][0];
		NS::Point3D  pt2 = controlPoints[1][1];
		NS::Vector3D v01(pt0, pt1);
		NS::Vector3D v02(pt0, pt2);
		NS::Vector3D n = v01.crossProduct(v02);
		n.normalize();

		//Matrix temp(3, 3);
		//temp.set(0, 0, pt0.x()); temp.set(0, 1, pt0.y()); temp.set(0, 2, pt0.z());
		//temp.set(1, 0, pt1.x()); temp.set(1, 1, pt1.y()); temp.set(1, 2, pt1.z());
		//temp.set(2, 0, pt2.x()); temp.set(2, 1, pt2.y()); temp.set(2, 2, pt2.z());
		//temp = temp.inverse();
		//Matrix rhs(3, 1);
		//rhs.set(0, 0, -1); rhs.set(1, 0, -1); rhs.set(2, 0, -1);
		//Matrix coe = temp*rhs;
		plane_coe.resize(4);
		plane_coe[0] = n.x(); plane_coe[1] = n.y(); plane_coe[2] = n.z();
		plane_coe[3] = -n.dotProduct(pt0);
	}
	convexHull = polyhedron.faces;
}
void Surface::surfacePointInv(const NS::Point3D &pt, NS::Point3D &projPt, double &u, double &v)
{
	if (this->is_plane) {
		NS::Point3D pt0 = controlPoints[0][0];
		NS::Point3D pt1 = controlPoints[1][0];
		NS::Point3D pt2 = controlPoints[1][1];
		projPt = NS::Geometry::ProjectPtOnPlane(pt, pt0, pt1, pt2);
		u = -2.0, v = -2.0;
	}
	else {
		double eps1 = 1e-10;
		double eps2 = 1e-10;
		// find initial u,v
		double u0, v0;
		projPt_initialGuessByLS(pt, u0, v0);

		std::vector <std::vector <NS::Vector3D> > SKL;
		SKL.resize(3);
		for (int i = 0; i < 3; i++) {
			SKL[i].resize(3);
			for (int j = 0; j < 3; j++)
				SKL[i][j].set(0.0, 0.0, 0.0);
		}
		double fu, fv, gu, gv, du, dv, fuv, guv;
		double u1, v1;
		NS::Vector3D r;
		for (int cnt = 0; cnt < 100; cnt++)
		{
			ratSurfaceDerives(2, u0, v0, SKL);
			//NurbsSurfaceDerivs(n, p, U, m, q, V, Pw, u0, v0, 2, SKL);
			r.set(SKL[0][0].x() - pt.x(), SKL[0][0].y() - pt.y(), SKL[0][0].z() - pt.z());
			fu = SKL[1][0].module() * SKL[1][0].module() + r * SKL[2][0];
			fv = SKL[1][0] * SKL[0][1] + r * SKL[1][1];
			gu = fv;
			gv = SKL[0][1].module() * SKL[0][1].module() + r * SKL[0][2];

			fuv = r * SKL[1][0];
			guv = r * SKL[0][1];

			double det = fu * gv - fv * gu;
			du = -(gv * fuv - fv * guv) / det;
			dv = -(-gu * fuv + fu * guv) / det;
			if (du != du) du = 0;

			if (dv != dv) dv = 0;
			u1 = u0 + du;
			v1 = v0 + dv;

			double er1 = abs(fuv) / SKL[1][0].module() / r.module();
			double er2 = abs(guv) / SKL[0][1].module() / r.module();
			NS::Point3D temp1(SKL[0][0].x(), SKL[0][0].y(), SKL[0][0].z());
			if ((temp1.distanceTo(pt)) <= eps1 || (er1 <= eps2 && er2 <= eps2))
			{
				break;
			}
			else
			{
				if (u1 < -1.0) u1 = -1.0;
				if (u1 > 1.0) u1 = 1.0;
				if (v1 < -1.0) v1 = -1.0;
				if (v1 > 1.0) v1 = 1.0;
				NS::Vector3D temp;
				temp = du * SKL[1][0] + dv * SKL[0][1];
				if (temp.module() < eps1)
				{
					u0 = u1;
					v0 = v1;
					break;
				}
			}
			u0 = u1;
			v0 = v1;
		}
		u = u0;
		v = v0;
		getCoord(u, v, projPt);
	}
}
void Surface::calculatePerImage(const NS::Point3D &pt, double &u, double &v) const {
	arma::mat M;
	calMrepVal(pt, M);
	double dist = arma::det(M * M.t());
	M = M.t();
}
void Surface::calculateLineIntecs(const NS::Point3D &ptOnLine, const NS::Vector3D &slope,
	std::vector<double> &intecUs, std::vector<NS::Point3D> &intecPts) const {
	_mRep[0].print("m0:");
	_mRep[1].print("m1:");
	_mRep[2].print("m2:");
	_mRep[3].print("m3:");
	bool debug = false;
	double TOL = 1e-6;
	arma::mat A, B;
	calMrepVal(ptOnLine, A); 
	calMrepVal(slope, B); B -= _mRep[0]; B = -B;
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
	int ro, ra12p;
	arma::mat U, V, A1, A11, A12, A12p, A1p, B1, B1p;
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
		ra12p = arma::rank(A12p, 1e-6);
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
	arma::eig_gen(eig_val, eig_vec, AB);
	eig_val.print("eigval:");
	arma::vec eig = arma::real(eig_val);
	NS::Point3D intec(ptOnLine + eig[0] * slope);
	arma::mat M;
	calMrepVal(intec, M);
	double dist = arma::det(M*M.t());
	//double u = calculatePerImage(intec);
	auto s = 1;

}
void Surface::calDist(const NS::Point3D &pt, NS::Point3D &projPt) {
	double u, v;
	surfacePointInv(pt, projPt, u, v);
}
bool Surface::lineIntersection(const NS::Point3D &nearPt, const NS::Point3D &seg0, const NS::Point3D &seg1, const double TOL, double &u0, double &v0, NS::Point3D &intPt) {
	// find initial u,v
	projPt_initialGuessByLS(nearPt, u0, v0);
	if (u0 == 1)
		u0 = 0.95;
	else if (u0 == -1)
		u0 = -0.95;
	if (v0 == 1)
		v0 = 0.95;
	else if (v0 == -1)
		v0 = -0.95;
	// find two orthogonal planes P1=(N1,d1),P2 = (N2,d2)
	NS::Vector3D lin(seg0, seg1);
	NS::Vector3D N0;
	abs(lin.y())<1e-8 ? N0.set(0, 1, 0) : N0.set(1, -lin.x() / lin.y(), 0);
	N0 /= N0.module();
	NS::Vector3D N1 = lin.crossProduct(N0);
	N1 /= N1.module();
	/*------------------------------------------------------------------------------------------------------------*/
	double d0 = -N0.dotProduct(seg0);
	double d1 = -N1.dotProduct(seg0);
	//
	std::vector <std::vector <NS::Vector3D> > SKL;
	SKL.resize(2);
	for (int i = 0; i < 2; i++) {
		SKL[i].resize(2);
		for (int j = 0; j < 2; j++)
			SKL[i][j].set(0.0, 0.0, 0.0);
	}
	ratSurfaceDerives(1, u0, v0, SKL);
	double R0 = N0.dotProduct(SKL[0][0]) + d0;
	double R1 = N1.dotProduct(SKL[0][0]) + d1;
	double du, dv;
	std::vector<double> R(2);
	std::vector<std::vector<double>> J(2);
	for (int i = 0; i < 2; i++)	J[i].resize(2);
	Matrix J_m,J_inv;
	for (int i = 0; i < 10; i++) {		
		if (R0 * R0 + R1 * R1 < TOL * TOL) {
			intPt = NS::Point3D(SKL[0][0].x(), SKL[0][0].y(), SKL[0][0].z());
			double k = (intPt - seg0) * (seg1 - seg0) / ((seg1 - seg0) * (seg1 - seg0));
			intPt = k * (seg1 - seg0) + seg0;
			return true;
		}

		R[0] = R0;	R[1] = R1;

		J[0][0] = N0.dotProduct(SKL[1][0]);
		J[0][1] = N0.dotProduct(SKL[0][1]);
		J[1][0] = N1.dotProduct(SKL[1][0]);
		J[1][1] = N1.dotProduct(SKL[0][1]);
		J_m = Matrix(J); 
		if (abs(J_m.calDet()) < 1e-10)
			return false;
		J_inv = J_m.inverse();
		du = -(J_inv.get(0, 0) * R[0] + J_inv.get(0, 1) * R[1]);
		dv = -(J_inv.get(1, 0) * R[0] + J_inv.get(1, 1) * R[1]);

		u0 += du; v0 += dv;
		if (u0 > 1)	u0 = 1;
		if (u0<-1)	u0 = -1;
		if (v0 > 1)	v0 = 1;
		if (v0<-1)	v0 = -1;
		ratSurfaceDerives(1, u0, v0, SKL);
		R0 = N0.dotProduct(SKL[0][0]) + d0;
		R1 = N1.dotProduct(SKL[0][0]) + d1;
		//if ( R0*R0 + R1*R1 > R[0] * R[0] + R[1] * R[1])
		//	return false;

	}
	return false;
}
bool Surface::lineIntersection(const NS::Point3D &nearPt, const NS::Point3D &seg0, const NS::Point3D &seg1, const double TOL, NS::Point3D &intPt) {
	double u0, u1;
	return lineIntersection(nearPt,seg0, seg1,TOL, u0, u1, intPt);
}
void Surface::discrete(const int n,std::vector<std::vector<int>> &faces, std::vector<NS::Point3D> &pts) const{
	pts.clear();
	pts.resize(n*n);
	faces.clear();
	// Init Knot Vector
	std::vector<double> k = {};
	for (int i = 0; i < n; i++)
		k.push_back(-1.0 + i * 2.0 / (n - 1));
	// Cal points for each pair of knot value
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++) 
			getCoord(k[i], k[j], pts[i*n+j]);
	for (int i = 0; i < n - 1; i++)
		for (int j = 0; j < n - 1; j++) {
			faces.push_back({
				i*n + j,
				i*n + j + 1,
				(i + 1)*n + j
			});
			faces.push_back({
				i*n + j + 1,
				(i + 1)*n + j + 1,
				(i + 1)*n + j
			});
		}
			

}
void Surface::offset(double dx, double dy, double dz) {
	for (size_t i = 0; i < controlPoints.size(); i++) {
		for (size_t j = 0; j < controlPoints[i].size(); j++) {
			controlPoints[i][j] += NS::Point3D(dx, dy, dz);
		}

	}
}

// parent surface
ParentSurface::ParentSurface(std::vector<NURBS::Basis> b, std::vector<std::vector<NS::Point3D>> cP) {
	assert(b.size() != 2);
	assert(b[1].nBasis != cP.size());
	for (std::vector<NS::Point3D> c : cP)
		assert(b[0].nBasis != c.size());
	basis = b; controlPoints = cP;
}

ParentSurface::ParentSurface(std::vector<NURBS::Basis> b, std::vector<std::vector<NS::Point3D>> cP,
	std::vector<std::vector<double>> w) {
	assert(b.size() == 2);
	assert(b[0].nBasis == cP.size());
	assert(cP.size() == w.size());
	for (std::vector<NS::Point3D> c : cP)
		assert(b[1].nBasis == c.size());
	for (std::vector<double> _w : w)
		assert(b[1].nBasis == _w.size());
	basis = b; controlPoints = cP; weightVector = w; rational = true;
}

void ParentSurface::subDivision() {
	children.clear();
	std::vector<int> npt(2);
	std::vector<int> nchild(2);
	for (int i = 0; i < 2; i++) {
		npt[i] = basis[i].order + 1;				// num of control points of the children
		nchild[i] = basis[i].nBasis - npt[i];		// num of the children (excluding multiple knot)
	}
	if (nchild[0] == 0 && nchild[1] == 0 && rational) // if there is no interior knot
		children.push_back(Surface(basis, controlPoints, weightVector));
	else if (nchild[0] == 0 && nchild[1] == 0 && !rational)
		children.push_back(Surface(basis, controlPoints));
	std::vector<std::vector<double>> iKnot(2);	// knot that need to be inserted
	std::vector<std::vector<int>> s(2);
	std::vector<int> idx = { 0,0 };				// idx of the knot to find the multi
	for (int j = 0; j < 2; j++) {
		if (nchild[j] == 0)
			continue;
		iKnot[j].push_back(basis[j].knotVector[npt[j]]);
		s[j].push_back(1);						// multiplicy
		for (size_t i = npt[j] + 1; i < basis[j].knotVector.size() - npt[j]; i++) {
			if (abs(iKnot[j].back() - basis[j].knotVector[i]) < absTOL)
				s[j][idx[j]] ++;
			else {
				idx[j]++;
				iKnot[j].push_back(basis[j].knotVector[i]);
				s[j].push_back(1);
			}
		}
		idx[j] ++; // number of the different knots to be inserted
	}
	Surface s3D(*this);
	// insert knots
	for (int j = 0; j < 2; j++) {
		if (nchild[j] == 0)
			continue;
		for (int i = 0; i < idx[j]; i++)
			if (basis[j].order - s[j][i] < 1)
				continue;
			else
				s3D.surfaceKnotIns(j, iKnot[j][i], basis[j].order - s[j][i]);
	}
	// construct basis for subsurfaces
	std::vector<double> V;
	std::vector<Basis> b(2);
	for (int j = 0; j < 2; j++) {
		for (int i = 0; i <= basis[j].order; i++) V.push_back(-1.0);
		for (int i = 0; i <= basis[j].order; i++) V.push_back(1.0);
		b[j] = Basis(V);
		V.clear();
	}
	// Extract subsurfaces
	for (int i = 0; i < idx[0] + 1; i++) {
		int b_idx1 = i * (npt[0] - 1);
		int e_idx1 = (i + 1) * (npt[0] - 1);

		for (int j = 0; j < idx[1] + 1; j++) {
			int b_idx2 = j*(npt[1] - 1);
			int e_idx2 = (j + 1) * (npt[1] - 1);

			std::vector<std::vector<NS::Point3D>> cp;
			std::vector<std::vector<double>> weight;

			for (int ii = b_idx1; ii <= e_idx1; ii++) {
				std::vector<NS::Point3D> cp_row;
				std::vector<double> weight_row;
				for (int jj = b_idx2; jj <= e_idx2; jj++) {
					cp_row.push_back(s3D.controlPoints[ii][jj]);
					weight_row.push_back(s3D.weightVector[ii][jj]);
				}
				cp.push_back(cp_row);
				weight.push_back(weight_row);
			}
			Surface temp(b, cp, weight);
			temp.getConvexHull();
			this->children.push_back(temp);
		}
	}
}
int ParentSurface::containsPoint(const NS::Point3D &pt, std::vector<int> &idx,double TOL_ratio) const {
	size_t nFace = convexHull.size();
	for (size_t i = 0; i < children.size(); i++) {
		std::vector<std::vector<NS::Point3D>> cp = children[i].controlPoints;
		std::vector<ConvexHull::face> cv = children[i].convexHull;			
		if (cv.size() == 0 && children[i].is_plane) { //plane
			std::vector<double> p_c = children[i].plane_coe;
			double dist = abs(p_c[0] * pt.x() + p_c[1] * pt.y() + p_c[2] * pt.z() + p_c[3]);
			if (dist<1e-4)
				idx.push_back(i);
			continue;
		}
		if (cv.size() == 0 && !children[i].is_plane)
			continue;
		for (size_t j = 0; j < cv.size(); j++) {
			double dist = cv[j].norm.X[0] * pt.x() + cv[j].norm.X[1] * pt.y() + cv[j].norm.X[2] * pt.z();
			if ( (cv[j].disc>0 && dist > cv[j].disc*(1+TOL_ratio)) || (cv[j].disc<=0 && dist > cv[j].disc/(1 + TOL_ratio))) { // outside
				break;
			}
			if (j == cv.size() - 1)
				idx.push_back(i);
		}
	}
	return idx.size();
}

void ParentSurface::offset(double dx, double dy, double dz) {
	for (size_t i = 0; i < controlPoints.size(); i++) {
		for (size_t j = 0; j < controlPoints[i].size(); j++) {
			controlPoints[i][j] += NS::Point3D(dx, dy, dz);
		}
	}
	for (Surface s : children)
		s.offset(dx, dy, dz);
}
NS_END