#include "basis.h";

NSN_BEG


// calculate the span for later calculating the basis function.
int Basis::findSpan(double u) const
{
	//n	Num of Basis - 1
	//p	Order of B-Spline Basis
	//u	The Knot
	//U	Knot Vector
	int low, high, mid;
	if (abs(u - knotVector[nBasis + 1]) < absTOL)
		return nBasis - 1;
	low = order;
	high = nBasis + 1;
	mid = (low + high) / 2;
	while (u < knotVector[mid] || u >= knotVector[mid + 1])
	{
		if (u < knotVector[mid])
			high = mid;
		else
			low = mid;
		mid = (low + high) / 2;
	}
	return(mid);
}
//constructor
Basis::Basis(std::vector<double> U) {
	knotVector = U;
	int length = U.size();
	order = 0;
	for (int i = 1; i < length; i++)
		if (abs(U[i] - U[0]) < absTOL)
			order++;
		else
			break;
	nBasis = length - order - 1;
}
void Basis::_basisFuns(const double u, const int i,std::vector<double> &N) const {
	assert(u >= -1 && u <= 1);
	N.resize(order + 1);
	N[0] = 1.0;
	std::vector<double> left(order + 1), right(order + 1);
	for (int j = 1; j <= order; j++) {
		left[j] = u - knotVector[i + 1 - j];
		right[j] = knotVector[i + j] - u;
		double saved = 0.0;
		for (int r = 0; r < j; r++) {
			double temp = N[r] / (right[r + 1] + left[j - r]);
			N[r] = saved + right[r + 1] * temp;
			saved = left[j - r] * temp;
		}
		N[j] = saved;
	}
}

void Basis::dersBasisFuns( const double u,const int i, const int dk,
					std::vector <std::vector <double> > &ders) const{
	ders.resize(dk + 1);
	for (std::vector<double> &t : ders)	t.resize(order + 1);
	int p = order;
	int n = dk;
	std::vector<double> U = knotVector;
	std::vector <std::vector <double> > ndu;
	ndu.resize(p + 1);
	for (int i = 0; i < p + 1; i++) ndu[i].resize(p + 1);
	ndu[0][0] = 1.0;
	std::vector<double> left(p + 1);
	std::vector<double> right(p + 1);

	for (int j = 1; j <= p; j++)
	{
		left[j] = u - U[i + 1 - j];
		right[j] = U[i + j] - u;
		double saved = 0.0;

		for (int r = 0; r < j; r++)
		{
			ndu[j][r] = right[r + 1] + left[j - r];
			double temp = ndu[r][j - 1] / ndu[j][r];
			ndu[r][j] = saved + right[r + 1] * temp;
			saved = left[j - r] * temp;
		}
		ndu[j][j] = saved;
	}

	for (int j = 0; j <= p; j++)
	{
		ders[0][j] = ndu[j][p];
	}
	std::vector <std::vector <double> >a;
	a.resize(2);
	a[0].resize(p + 1); a[1].resize(p + 1);

	int s1, s2;
	for (int r = 0; r <= p; r++)
	{
		s1 = 0, s2 = 1;
		a[0][0] = 1.0;
		int j1, j2;
		int j;
		for (int k = 1; k <= n; k++)
		{
			double d = 0.0;
			int rk = r - k;
			int pk = p - k;
			if (r >= k)
			{
				a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
				d = a[s2][0] * ndu[rk][pk];
			}

			if (rk >= -1) j1 = 1;
			else j1 = -rk;
			if (r - 1 <= pk) j2 = k - 1;
			else j2 = p - r;
			for (j = j1; j <= j2; j++)
			{
				a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
				d += a[s2][j] * ndu[rk + j][pk];
			}
			if (r <= pk)
			{
				a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
				d += a[s2][k] * ndu[r][pk];
			}
			ders[k][r] = d;
			j = s1; s1 = s2; s2 = j;
		}
	}

	int r = p;
	for (int k = 1; k <= n; k++)
	{
		for (int j = 0; j <= p; j++) ders[k][j] = ders[k][j] * r;
		r *= (p - k);
	}
}

// normlize the knot vector to -1 to 1
void Basis::normlized(){
	double u0 = knotVector[0];
	double u1 = knotVector[knotVector.size() - 1];
	double k = 2.0 / (u1 - u0);
	double b = 1 - u1 * k;
	for (int i = 0; i < knotVector.size(); i++) {
		if (abs(knotVector[i] - u0) < absTOL)
			knotVector[i] = -1.0;
		else if (abs(knotVector[i] - u1) < absTOL)
			knotVector[i] = 1.0;
		else
			knotVector[i] = knotVector[i] * k + b;
	}
}

void Basis::setRange(double ubegin, double uend) {
	if (abs(knotVector[0] + 1) > absTOL || abs(knotVector[knotVector.size() - 1] - 1) > absTOL) {
		double u0 = knotVector[0];
		double u1 = knotVector[knotVector.size() - 1];
		double k = 2.0 / (u1 - u0);
		double b = 1 - u1 * k;
		ubegin = ubegin * k + b;
		uend = uend * k + b;
		if (abs(ubegin + 1) < absTOL)
			ubegin = -1.0;
		if (abs(uend - 1) < absTOL)
			uend = 1.0;
		this->normlized();
	}
	u_begin = ubegin;
	u_end = uend;
}

void Basis::reverse() {
	std::reverse(knotVector.begin(), knotVector.end());
	normlized();
}
// insert knot u r times,return findspan(u)
int Basis::insert(const double u,const int r) {
	int k = findSpan(u);
	knotVector.insert(knotVector.begin() + k + 1, r, u);
	nBasis += r;
	return k;
}
// return shapefun N which starts from i0-th to i1-th shape functions (including i0 and i1)
void Basis::basisFuns(const double u, std::vector<double> &N,int &i0, int &i1) const {
	int iSpan = findSpan(u);
	_basisFuns(u, iSpan, N);
	i0 = iSpan - order;
	i1 = iSpan;
}
//Export Basis plot
void Basis::exportBasisScatter(std::vector<std::vector<NS::Point3D>> &scatters) {
	int N = 100;
	double *k = new double[N];
	scatters.resize(nBasis);
	for (int i = 0; i <= N; i++) {
		double lclKnot = 2.0 / N * i - 1;
			
		std::vector<double> N;
		int i0, i1;
		basisFuns(lclKnot, N, i0, i1);

		for (int j = 0; j < nBasis; j++) {
			if (j >= i0 && j <= i1)
				scatters[j].push_back(NS::Point3D(lclKnot, N[j - i0], 0.0));
			else
				scatters[j].push_back(NS::Point3D(lclKnot, 0.0, 0.0));
		}
	}
}
NS_END