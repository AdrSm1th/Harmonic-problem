//sparse matrix.h

#pragma once

#include <vector>
#include <stdexcept>
#include "mesh.h"

class Block {
public:
	double p_;
	double c_;
	Block() : p_(0), c_(0) {};
	Block(double p, double c) : p_(p), c_(c) {}

	Block operator-(const Block &block) const {
		return Block(p_ - block.p_, c_ - block.c_);
	}

	Block operator-=(const Block &block) {
		p_ -= block.p_;
		c_ -= block.c_;
		return *this;
	}

	Block operator+(const Block &block) {
		return Block(p_ + block.p_, c_ + block.c_);
	}

	Block operator*(const Block &block) const {
		return Block(p_ * block.p_ - c_ * block.c_,
			     p_ * block.c_ + c_ * block.p_);
	}

	Block operator*(const double &c) const {
		return Block(p_ * c, c_ * c);
	}

	Block operator*=(const Block &block) {
		double new_p = p_ * block.p_ - c_ * block.c_;
		double new_c = p_ * block.c_ + c_ * block.p_;
		p_ = new_p;
		c_ = new_c;
		return *this;
	}

	Block operator/(const Block &b) const {
		double den = b.p_ * b.p_ + b.c_ * b.c_;
		return Block(
			(p_ * b.p_ + c_ * b.c_) / den,
			(c_ * b.p_ - p_ * b.c_) / den
		);
	}

	Block operator/=(const Block &block) {
		double den = block.p_ * block.p_ + block.c_ * block.c_;
		p_ = (p_ * block.p_ + c_ * block.c_) / den;
		c_ = (c_ * block.p_ - p_ * block.c_) / den;
		return *this;
	}

	bool operator<(const double &number) const {
		double det = p_ * p_ + c_ * c_;
		return det < number ? true : false;
	}
};

class BlockVector {
public:
	double p_;
	double c_;
	BlockVector() : p_(0), c_(0) {};
	BlockVector(double p, double c) : p_(p), c_(c) {}

	BlockVector &operator=(const BlockVector &other) = default;

	BlockVector operator-(const Block &block) const {
		return BlockVector(p_ - block.p_, c_ - block.c_);
	}

	BlockVector operator-=(const BlockVector &block) {
		p_ -= block.p_;
		c_ -= block.c_;
		return *this;
	}

	BlockVector operator+(const BlockVector &block) const {
		return BlockVector(p_ + block.p_, c_ + block.c_);
	}

	BlockVector operator*(double c) const {
		return BlockVector(p_ * c, c_ * c);
	}
};

class BlockCSRMatrix {
private:
	std::vector<Block> al_;
	std::vector<Block> au_;
	std::vector<Block> di_;
	std::vector<int> ia_;
	std::vector<int> ja_;
	std::vector<int> firstCol_;   // первый столбец профиля для каждой строки
	int index(int i, int j) const;

public:
	BlockCSRMatrix(Mesh3D &mesh);
	BlockCSRMatrix() : ia_(std::vector<int>(0)) { };
	void addBlock(int row, int col, Block block);
	void changeBlock(int row, int col, Block block);
	void multiply(std::vector<BlockVector> &x, std::vector<BlockVector> &y);
	int getSize();
	void solve(std::vector<BlockVector> &b, std::vector<BlockVector> &x);

	Block& operator()(int i, int j) {
		if (i < j) {
			int idx = index(j, i);
			if (idx == -1) {
				Block b(0, 0);
				return(b);
			}
			return au_[idx];
		}
		else if (i > j) {
			int idx = index(i, j);
			if (idx == -1) {
				Block b(0, 0);
				return(b);
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
				return(Block(0, 0));
			}
			return au_[idx];
		}
		else if (j > i) {
			int idx = index(i, j);
			if (idx == -1) {
				return(Block(0, 0));
			}
			return al_[idx];
		}
		else {
			return di_[i];
		}
	}
};