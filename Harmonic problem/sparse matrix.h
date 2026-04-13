//sparse matrix.h

#pragma once

#include <vector>
#include "mesh.h"

class Block {
public:
	double p_;
	double c_;
	Block() : p_(0), c_(0){};
	Block(double p, double c) : p_(p), c_(c) {}
	
	Block operator-(const Block &block) const {
		return Block(p_ - block.p_, c_ - block.c_);
	}

	Block operator-=(const Block &block) {
		p_ -= block.p_;
		c_ -= block.c_;
		return *this;
	}

	Block operator*(const Block &block) const {
		return Block(p_ * block.p_, c_ * block.c_);
	}

	Block operator*=(const Block &block) {
		p_ *= block.p_;
		c_ *= block.c_;
		return *this;
	}

	Block operator/(const Block &block) const {
		return Block((p_ * block.p_ + c_ * block.c_) / (block.p_ * block.p_ + block.c_ * block.c_),
			(p_ * block.c_ - block.p_ * c_) / (block.p_ * block.p_ + block.c_ * block.c_));
	}

	Block operator/=(const Block &block) {
		p_ = (p_ * block.p_ + c_ * block.c_) / (block.p_ * block.p_ + block.c_ * block.c_);
		c_ = (p_ * block.c_ - block.p_ * c_) / (block.p_ * block.p_ + block.c_ * block.c_);
		return *this;
	}

	bool operator<(const double &number) const {
		double det = p_ * p_ + c_ * c_;
		return det < number ? true : false;
	}
};

class BlockCSRMatrix {
private:
	std::vector<Block> al_;
	std::vector<Block> au_;
	std::vector<Block> di_;
	std::vector<int> ia_;
	std::vector<int> ja_;
	int index(int i, int j) const;

public:
	BlockCSRMatrix(Mesh3D &mesh);
	BlockCSRMatrix() : ia_(std::vector<int>(0)) { };
	void addBlock(int row, int col, Block &block);
	void multiply(std::vector<double> &x, std::vector<double> &y);
	void convertToProfile();
	void getLU();

	Block& operator()(int i, int j) {
		if (i < j) {
			int idx = index(j, i);
			if (idx == -1) {
				throw std::invalid_argument("Element not found");
			}
			return au_[idx];
		}
		else if (j > i) {
			int idx = index(i, j);
			if (idx == -1) {
				throw std::invalid_argument("Element not found");
			}
			return al_[idx];
		}
		else {
			return di_[i];
		}
	}

	const Block operator()(int i, int j) const {
		if (i < j) {
			int idx = index(j, i);
			if (idx == -1) {
				throw std::invalid_argument("Element not found");
			}
			return au_[idx];
		}
		else if (j > i) {
			int idx = index(i, j);
			if (idx == -1) {
				throw std::invalid_argument("Element not found");
			}
			return al_[idx];
		}
		else {
			return di_[i];
		}
	}
};