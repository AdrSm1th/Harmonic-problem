//boundary conditions.cpp

#include <fstream>
#include "boundary conditions.h"

BCManager::BCManager(Mesh3D &mesh, const char *filename, BlockCSRMatrix &matrix) {
	mesh_ = &mesh;
	readBoundaryConditions(filename);
	matrix_ = matrix;
	boundaryFaces_ = mesh_->getBoundaryFaces();
}

void BCManager::ApplyDirichle() {
	for (int i = 0; i < 8; i++) {
		if (boundaryConditions_[i].type != 1) continue;
		for (Face &face : boundaryFaces_) {
			if (face.localFaceId != i) continue;

		}
	}
}

void BCManager::readBoundaryConditions(const char *filename) {
	std::ifstream input(filename);
	for (int i = 0; i < 8; i++)
	{
		input >> boundaryConditions_[i].type;

		if (boundaryConditions_[i].type == 3) input >> boundaryConditions_[i].beta;
	}
}

double BoundaryFunctions::u_s(double x, double y, double z) {
	return 0;
}

double BoundaryFunctions::u_c(double x, double y, double z) {
	return 0;
}

double BoundaryFunctions::theta_s(double x, double y, double z) {
	return 0;
}

double BoundaryFunctions::theta_c(double x, double y, double z) {
	return 0;
}

double BoundaryFunctions::ubeta_s(double x, double y, double z) {
	return 0;
}

double BoundaryFunctions::ubeta_c(double x, double y, double z) {
	return 0;
}