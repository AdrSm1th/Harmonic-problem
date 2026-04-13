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
	BlockCSRMatrix matrix_;
	std::vector<Block> global_b_;

	std::array<Block, 64> A_local_ = {};
	std::array<Block, 8> b_local_ = {};

public:
	HarmonicAssembler(Mesh3D &mesh, BlockCSRMatrix &matrix, std::vector<Block> &global_b);
	void computeLocalMatrices(int elemId);
	void addToGlobalMatrix();
};