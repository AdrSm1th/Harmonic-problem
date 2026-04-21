//assembler.h

#pragma once

#include <array>
#include "mesh.h"
#include "quadrature.h"
#include "sparse matrix.h"

class HarmonicAssembler {
private:
	Mesh3D mesh_;
	Quadrature quadrature_;
	Coefficients coefficients_;

	std::array<Block, 64> A_local_ = {};
	std::array<BlockVector, 8> b_local_ = {};

	void computeLocalMatrices(int elemId);

public:
	HarmonicAssembler(Mesh3D &mesh);

	void assembleSystem(BlockCSRMatrix &matrix, std::vector<BlockVector> &global_b);
};