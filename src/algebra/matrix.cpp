#include "matrix.h"

NSN_BEG



double Matrix::calDet() {
	//return calAdj(-1, -1);
	assert(col == row);
	return determ(values, col);
	//if (col == 2)
	//	return values[0][0] * values[1][1] - values[1][0] * values[0][1];
	//if (col == 3) {
	//	double det = 0;
	//	for (size_t i = 0; i < row; i++) {
	//		double temp = 1.0;
	//		for (size_t j = 0; j < col; j++) {
	//			temp *= i + j < row ? values[i + j][j] : values[i + j - row][j];
	//		}
	//		det += temp;
	//	}
	//	for (size_t i = 0; i < row; i++) {
	//		double temp = 1.0;
	//		for (size_t j = 0; j < col; j++) {
	//			//bool justf = (i >= j);
	//			//int idx = i - j >= 0 ? i - j : i - j + col;
	//			temp *= i >= j ? values[i - j][j] : values[col + i - j][j];
	//		}
	//		det -= temp;
	//	}
	//	return det;
	//}
}

double Matrix::determ(const std::vector<std::vector<double>> &a, const int &n) {
	std::vector<std::vector<double>> temp(row);
	for (size_t i = 0; i < row; i++)	temp[i].resize(row);
	double det = 0.0;
	int  p, h, k, i, j;
	if (n == 1) {
		return values[0][0];
	}
	else if (n == 2) {
		det = (a[0][0] * a[1][1] - a[0][1] * a[1][0]);
		return det;
	}
	else {
		for (p = 0; p<n; p++) {
			h = 0;
			k = 0;
			for (i = 1; i<n; i++) {
				for (j = 0; j<n; j++) {
					if (j == p) {
						continue;
					}
					temp[h][k] = a[i][j];
					k++;
					if (k == n - 1) {
						h++;
						k = 0;
					}
				}
			}
			det = det + a[0][p] * pow(-1, p)*determ(temp, n - 1);
		}
		return det;
	}
}

double Matrix::calAdj(const int rIdx, const int cIdx) {
	assert(col == row);
	if (row == 2 && (rIdx + cIdx) % 2 == 0)
		return values[1 - cIdx][1 - rIdx];
	else if (row == 2 && (rIdx + cIdx) % 2 != 0)
		return -values[1 - cIdx][1 - rIdx];

	std::vector<std::vector<double>> _val;
	_val = values;
	_val.erase(_val.begin() + rIdx);
	for (std::vector<double> &r : _val)
		r.erase(r.begin() + cIdx);
	Matrix temp(_val);
	double res = temp.calDet();
	if ((rIdx + cIdx) % 2 == 0) {
		return res;
	}
	else
		return -res;
	
}

Matrix Matrix::inverse() {
	double det = this->calDet();
	assert( abs(det) > absTOL);
	std::vector<std::vector<double>> val;
	for (size_t i = 0; i < this->row; i++) {
		std::vector<double> val_row;
		for (size_t j = 0; j < this->col; j++) {
			double adj = calAdj(i, j);
			val_row.push_back(adj/ det);
		}
		val.push_back(val_row);
	}
	return Matrix(val);
}

Matrix Matrix::transpose() {
	std::vector<std::vector<double>> val;
	for (size_t i = 0; i < col; i++) {
		std::vector<double> val_row;
		for (size_t j = 0; j < row; j++) {
			val_row.push_back(this->values[j][i]);
		}
		val.push_back(val_row);
	}
	return Matrix(val);
}

Matrix Matrix::operator* (const Matrix &m) {
	assert(this->col = m.row);
	std::vector<std::vector<double>> val;
	for (size_t i = 0; i < row; i++) {
		std::vector<double> val_row;
		for (size_t j = 0; j < m.col; j++) {
			double temp = 0;
			for (size_t k = 0; k < col; k++) {
				temp += values[i][k] * m.values[k][j];
			}
			val_row.push_back(temp);
		}
		val.push_back(val_row);
	}
	Matrix tmp(val);
	return tmp;
}

Matrix Matrix::operator= (const Matrix &m) {
	values = m.values;
	row = m.row;
	col = m.col;
	return *this;
}

Matrix::Matrix(const int _row, const int _col) {
	std::vector<std::vector<double>> val(_row, std::vector<double>(_col, 0));
	values = val;
	row = _row;
	col = _col;
}
Matrix::Matrix(std::vector<std::vector<double>> &val) {
	values = val;
	row = val.size();
	col = val[0].size();
	for (size_t i = 0; i < row; i++)
		assert(col == val[i].size());
}
NS_END