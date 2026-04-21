//iteration solver.h

#pragma once

#include <vector>
#include <stdexcept>
#include "sparse matrix.h"

//class MyVector {
//private:
//	std::vector<double> data_;
//	int n_ = 0;
//
//public:
//	MyVector(int size) : n_(size) {};
//	int getSize();
//
//	const double& operator[](int i) const{
//		if (i >= n_) throw std::out_of_range("MyVector: index out of range");
//		return data_[i];
//	}
//
//	double &operator[](int i) {
//		if (i >= n_) throw std::out_of_range("MyVector: index out of range");
//		return data_[i];
//	}
//
//	MyVector operator+=(MyVector v) {
//		if (v.getSize() != n_) throw std::invalid_argument("MyVector+: different size");
//		for (int i = 0; i < n_; data_[i] += v[i], i++);
//	}
//};

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