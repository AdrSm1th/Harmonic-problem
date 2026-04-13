//boundary conditions.h

#pragma once

#include <array>
#include <vector>
#include "mesh.h"
#include "sparse matrix.h"

class BCManager {
private:
	Mesh3D *mesh_;
	std::array<BoundaryCondition, 8> boundaryConditions_;
	BlockCSRMatrix matrix_;
	std::vector<Face> boundaryFaces_;
	void ApplyDirichle();
	void ApplyNeumann();
	void ApplyRobin();

public:
	BCManager(Mesh3D &mesh, const char *filename, BlockCSRMatrix &matrix);
	void ApplyBC();
};