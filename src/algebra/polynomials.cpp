#include "polynomials.h"
namespace NURBS {
	int PolynomialsFunctor::binoCoe(static int n, static int k) {
		if (n == k || k == 0) {
			return 1;
		}
		else if (n == 0)
			return 0;
		else {
			return binoCoe(n - 1, k) + binoCoe(n - 1, k - 1);
		}
	}
}