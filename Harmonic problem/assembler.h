//assembler.h

#pragma once

#include <array>;
#include "mesh.h"
#include "quadrature.h"

struct Elem {
	double p;
	double c;
};

class HarmonicAssembler {
private:
	Mesh3D mesh_;
	Quadrature quadrature_;
	Coefficients coefficients_;
	std::array<Elem, 64> A_local_ = {};
	std::array<Elem, 8> b_local_ = {};

public:
	HarmonicAssembler(Mesh3D &mesh);
	void ComputeLocalMatrices(int elemId);
	void AddToGlobalMatrix();
};