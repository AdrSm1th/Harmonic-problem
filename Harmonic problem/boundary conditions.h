//boundary conditions.h

#pragma once

#include <array>
#include <vector>
#include "mesh.h"
#include "sparse matrix.h"

class BoundaryFunctions {
public:
	double u_s(double x, double y, double z);
	double u_c(double x, double y, double z);
	double theta_s(double x, double y, double z);
	double theta_c(double x, double y, double z);
	double ubeta_s(double x, double y, double z);
	double ubeta_c(double x, double y, double z);
};

class BCManager {
private:
	Mesh3D *mesh_;
	std::array<BoundaryCondition, 8> boundaryConditions_;
	BlockCSRMatrix matrix_;
	std::vector<Face> boundaryFaces_;
	void readBoundaryConditions(const char *filename);
	void ApplyDirichle();
	void ApplyNeumann();
	void ApplyRobin();

public:
	BCManager(Mesh3D &mesh, const char *filename, BlockCSRMatrix &matrix);
	void ApplyBC();
};