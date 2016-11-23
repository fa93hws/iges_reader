#pragma once
#include <vector>
#include <stdexcept>
#include <assert.h>
#include "..\shared\config.h"


NSN_BEG

class Matrix {
private:
	double determ(const std::vector<std::vector<double>> &a, const int &n);
	const double absTOL = 1e-6;
public:
	size_t row;
	size_t col;
	std::vector<std::vector<double>> values;
	

	Matrix(std::vector<std::vector<double>> &val);
	Matrix(const int _row, const int _col);
	Matrix() {}
	double get(const int r, const int c) { return values[r][c]; };
	void set(const int r, const int c, const double v) { values[r][c] = v; }

	Matrix inverse();
	Matrix transpose();	
	double calDet();
	double calAdj(const int rId, const int cId);

	Matrix operator* (const Matrix &m);
	Matrix operator= (const Matrix &m);
};
NS_END