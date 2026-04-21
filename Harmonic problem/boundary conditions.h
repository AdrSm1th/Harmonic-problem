//boundary conditions.h

#pragma once

#include <array>
#include <vector>
#include "mesh.h"
#include "sparse matrix.h"

class BoundaryFunctions {
public:
	static double u_s(double x, double y, double z);
	static double u_c(double x, double y, double z);
};

class BCManager {
private:
	Mesh3D *mesh_;
	BlockCSRMatrix matrix_;
	std::vector<Face> boundaryFaces_;
	std::vector<Block> b_;

public:
	BCManager(Mesh3D &mesh, const char *filename, BlockCSRMatrix &matrix, std::vector<Block> &b);
	void ApplyDirichle();
};