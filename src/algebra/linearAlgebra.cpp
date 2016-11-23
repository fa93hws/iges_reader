#include "linearAlgebra.h"
namespace NURBS {
	void LinearAlgebra::leastSquare(Matrix &A, Matrix &B, Matrix &coe) {
		// A * coe = B;
		assert(A.row >= A.col);
		assert(A.row == B.row);
		if (A.row == A.col) {
			Matrix temp = A.inverse();
			coe = temp * B;
		}
		else {
			Matrix temp = A;
			temp = A.transpose();
			temp = temp * A;
			double det = temp.calDet();
			temp = temp.inverse();
			temp = temp * A.transpose();
			coe = temp * B;
		}
	}
}