//iteration solver.h

#pragma once

#include <vector>
#include <stdexcept>
#include "sparse matrix.h"

class LOSsolver {
private:
	double eps_;
	int maxiter_;
	std::vector<BlockVector> precondition(std::vector<Block> &M, std::vector<BlockVector> &vector) const;
	double dotProduct(std::vector<BlockVector> &v1, std::vector<BlockVector> &v2) const;
	std::vector<BlockVector> VplusV(std::vector<BlockVector> &v1, std::vector<BlockVector> &v2) const;
	std::vector<BlockVector> CmultV(double &c, std::vector<BlockVector> &v) const;

public:
	LOSsolver(double eps, int maxiter) : eps_(eps), maxiter_(maxiter) {}
	void solve(BlockCSRMatrix &matrix, std::vector<BlockVector> &b, std::vector<BlockVector> &x);
};