//assembler.h

#pragma once

#include <array>
#include "mesh.h"
#include "quadrature.h"
#include "sparse matrix.h"

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

	BlockVector operator+(const BlockVector &block) const{
		return BlockVector(p_ + block.p_, c_ + block.c_);
	}

	BlockVector operator*(double c) const {
		return BlockVector(p_ * c, c_ * c);
	}
};

class HarmonicAssembler {
private:
	Mesh3D mesh_;
	Quadrature quadrature_;
	Coefficients coefficients_;
	BlockCSRMatrix matrix_;
	std::vector<BlockVector> global_b_;

	std::array<Block, 64> A_local_ = {};
	std::array<BlockVector, 8> b_local_ = {};

public:
	HarmonicAssembler(Mesh3D &mesh, BlockCSRMatrix &matrix, std::vector<BlockVector> &global_b);
	void computeLocalMatrices(int elemId);
	void addToGlobalMatrix();
};