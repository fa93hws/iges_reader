#include "MatrixRep.h"


NSN_BEG
int MatrixRep::calCombiation(const int n, const int k)const { 
	if (n == k || k==0) return 1;
	if (k == 1) return n;
	assert(n > k);
	int denu = k;
	if (n - k < k) denu = n - k;
	int out = 1;
	for (int i = 0; i < denu; i++) out *= (n - i);
	for (int i = 2; i <= denu; i++) out /= i;
	return out;
}
void MatrixRep::getCompForF_curve(const Curve &c ,const int i, std::vector<double> &comp) const {
	if (c.rational) comp.push_back(c.weightVector[i]);
	else comp.push_back(1.0);
	comp.push_back(c.controlPoints[i].x());
	comp.push_back(c.controlPoints[i].y());
	comp.push_back(c.controlPoints[i].z());
}

MatrixRep::MatrixRep(const Curve &c,const int v) {
	int order = c.getOrder();
	assert(order == c.controlPoints.size() - 1);
	std::vector<arma::mat> mmat(4);
	// cal S
	arma::mat  smat(order + v + 1, 4 * v + 4);
	smat.zeros();
	// finish cal comb,cal S
	for (int i = 0; i <= order; i++) {
		std::vector<double> comp;
		getCompForF_curve(c, i,comp);
		for (int j = 0; j <= v; j++) {
			double comb = calCombiation(v, j);
			comb *= calCombiation(order, i);
			comb /= calCombiation(order + v, i + j);
			for (size_t k = 0; k < comp.size(); k++)
				smat(i + j, k * (1 + v) + j) = comb * comp[k];
		}
	}
	// finish cal S, cal M
	arma::mat m = null(smat);
	for (int i = 0; i < 4; i++) {
		mmat[i] = m.rows(2 * i, 2 * i + 1);
		mmat[i].debug_set_size();
	}
}
NS_END