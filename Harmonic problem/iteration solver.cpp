//iteration solver.cpp

#include <iostream>
#include "iteration solver.h"

std::vector<BlockVector> LOSsolver::precondition(std::vector<Block> &M, std::vector<BlockVector> &vector) const {
	int n = (int)M.size();
	std::vector<BlockVector> result(n);

	for (int i = 0; i < n; i++) {
		double det = M[i].p_ * M[i].p_ + M[i].c_ * M[i].c_;
		if (abs(det) < 1e-12) throw std::runtime_error("Zero diagonal block determinant");
		result[i].p_ = (M[i].p_ * vector[i].p_ + M[i].c_ * vector[i].c_) / det;
		result[i].c_ = (-M[i].c_ * vector[i].p_ + M[i].p_ * vector[i].c_) / det;
	}

	return result;
}

double LOSsolver::dotProduct(std::vector<BlockVector> &v1, std::vector<BlockVector> &v2) const {
	double result = 0;
	int n = (int)v1.size();

	for (int i = 0; i < n; result += v1[i].p_ * v2[i].p_ + v1[i].c_ * v2[i].c_, i++);

	return result;
}

std::vector<BlockVector> LOSsolver::VplusV(std::vector<BlockVector> &v1, std::vector<BlockVector> &v2) const {
	int n = (int)v1.size();
	std::vector<BlockVector> result(n);

	for (int i = 0; i < n; result[i] = v1[i] + v2[i], i++);

	return result;
}

std::vector<BlockVector> LOSsolver::CmultV(double &c, std::vector<BlockVector> &v) const {
	int n = (int)v.size();
	std::vector<BlockVector> result(n);

	for (int i = 0; i < n; result[i] = v[i] * c, i++);

	return result;
}

void LOSsolver::solve(BlockCSRMatrix &matrix, std::vector<BlockVector> &b, std::vector<BlockVector> &x) {
	int n = (int)b.size(),i = 1;
	std::vector<BlockVector> r(n), z(n), p(n);
	std::vector<Block> M(n);
	double a = 0, beta = 0, discrepancy = 1;
	std::fill(x.begin(), x.end(), BlockVector{ 0.0, 0.0 });

	for (int i = 0; i < n; i++) {
		M[i].p_ = matrix(i, i).p_;
		M[i].c_ = matrix(i, i).c_;
	}

	r = precondition(M, b);
	z = precondition(M, r);
	std::vector<BlockVector> y(n);
	matrix.multiply(z, y);
	p = precondition(M, y);;
	
	do
	{
		a = dotProduct(p, r) / dotProduct(p, p);
		std::vector<BlockVector> az = CmultV(a, z);
		x = VplusV(x, az);

		discrepancy = sqrt(dotProduct(r, r)) / sqrt(dotProduct(b, b));
		std::cout << "i: " << i << " " << discrepancy << std::endl;
		double ma = -a;
		std::vector<BlockVector> map = CmultV(ma, p);
		r = VplusV(r, map);
		std::vector<BlockVector> Mr = precondition(M, r);
		std::vector<BlockVector> t(n);
		matrix.multiply(Mr, t);
		Mr = precondition(M, t);

		beta = dotProduct(p, Mr) / dotProduct(p, p);
		std::vector<BlockVector> pMr = precondition(M, r), betaz = CmultV(beta, z), betap = CmultV(beta, p);
		z = VplusV(pMr, betaz);
		p = VplusV(Mr, betap);

		i++;
	} while (abs(discrepancy) > eps_ && i < maxiter_);
}