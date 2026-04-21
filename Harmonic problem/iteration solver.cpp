//iteration solver.cpp

#include "iteration solver.h"

std::vector<BlockVector> LOSsolver::precondition(std::vector<Block> &M, std::vector<BlockVector> &vector) const {
	int n = (int)M.size();
	std::vector<BlockVector> result(n);

	for (int i = 0; i < n; i++) {
		result[i].p_ = vector[i].p_ * (M[i].p_ + M[i].c_ * M[i].c_ / M[i].p_) - vector[i].c_ * M[i].c_;
		result[i].c_ = vector[i].p_ * M[i].p_;
	}

	return result;
}

//double LOSsolver::dotProduct(std::vector<Block> &v1, std::vector<Block> &v2) const {
//	double result = 0;
//	int n = (int)v1.size();
//
//	for (int i = 0; i < n; result += v1[i].p_ * v2[i].p_ + v1[i].c_ * v2[i].c_);
//
//	return result;
//}
//
//std::vector<Block> LOSsolver::VplusV(std::vector<Block> &v1, std::vector<Block> &v2) const {
//	int n = (int)v1.size();
//	std::vector<Block> result(n);
//
//	for (int i = 0; i < n; result[i] = v1[i] + v2[i], i++);
//
//	return result;
//}
//
//std::vector<Block> LOSsolver::CmultV(double &c, std::vector<Block> &v) const {
//	int n = (int)v.size();
//	std::vector<Block> result(n);
//
//	for (int i = 0; i < n; result[i] = v[i] * c, i++);
//
//	return result;
//}

void LOSsolver::solve(BlockCSRMatrix &matrix, std::vector<double> &b, std::vector<double> &x) {
	int n = (int)b.size();
	std::vector<Block> r(n), z(n), p(n), M(n);
	double a = 0, beta = 0, discrepancy = 1;

	for (int i = 0; i < n; i++) {
		M[i].p_ = sqrt(matrix(i, i).p_);
		M[i].c_ = sqrt(matrix(i, i).c_);
	}
	
	//for (int i = 0; i < n; i++) {
	//	r[i].
	//}
}