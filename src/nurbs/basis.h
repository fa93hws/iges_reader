#pragma once
#include "..\shared\config.h"
#include <vector>
#include <assert.h>
#include "../geo/point3d.h"

NSN_BEG
class Basis {
private:
	double absTOL = 1e-6;
	void _basisFuns(const double u, const int i,std::vector<double> &N) const;
		
	void normlized();

public:
	//variables
	std::vector<double> knotVector;
	//vector<double> weightVector;
	int order;
	int nBasis;
	double u_begin = -1;
	double u_end = 1;

	//Constructor
	Basis() {};
	Basis(std::vector<double> U);

	//simple algorithm
	int findSpan(double u) const;
	void setRange(double ubegin, double uend);
	int insert(const double u, const int r);
	void dersBasisFuns(const double u, const int i,const int dk,
		std::vector <std::vector <double> > &ders) const;
	void reverse();
	// return shapefun N which starts from i0-th to i1-th shape functions
	void basisFuns(const double u, std::vector<double> &N,int &i0, int &i1) const;

	//export plotting
	void exportBasisScatter(std::vector<std::vector<NS::Point3D>> &scatters);

};
NS_END