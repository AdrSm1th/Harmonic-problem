//boundary conditions.cpp

#include "boundary conditions.h"

BCManager::BCManager(Mesh3D &mesh, const char *filename, BlockCSRMatrix &matrix) {
	mesh_ = &mesh;
	mesh_->readBoundaryConditions(filename);
	boundaryConditions_ = mesh_->getBoundaryConditions();
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