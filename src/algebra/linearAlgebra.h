#pragma once
#include "matrix.h"
#include <stdexcept>
#include <assert.h>

namespace NURBS {
	class LinearAlgebra {
	public:
		LinearAlgebra() {};

		static void leastSquare(Matrix &A, Matrix &B, Matrix &coe);
	};
}