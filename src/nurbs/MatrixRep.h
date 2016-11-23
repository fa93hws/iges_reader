#include "..\shared\config.h"
#include "curve.h"
#include "surface.h"
#include "basis.h"
#include <vector>
#include "..\algebra\matrix.h"
#include <armadillo>

NSN_BEG
class MatrixRep {
private:
	int calCombiation(const int a, const int b) const;
	void getCompForF_curve(const Curve &c, const int i,std::vector<double> &comp) const;
public:
	MatrixRep(const Curve &c,const int v);

};
NS_END