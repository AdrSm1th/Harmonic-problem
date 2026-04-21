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
	std::vector<BlockVector> precondition(std::vector<Block> &M, std::vector<BlockVector> &vector) const;
	//double dotProduct(std::vector<Block> &v1, std::vector<Block> &v2) const;
	//std::vector<Block> VplusV(std::vector<Block> &v1, std::vector<Block> &v2) const;
	//std::vector<Block> CmultV(double &c, std::vector<Block> &v) const;

public:
	LOSsolver(double eps) : eps_(eps) {}
	void solve(BlockCSRMatrix &matrix, std::vector<BlockVector> &b, std::vector<BlockVector> &x);
};