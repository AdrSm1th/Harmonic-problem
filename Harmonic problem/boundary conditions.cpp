//boundary conditions.cpp

#include <fstream>
#include "boundary conditions.h"
#include "basis.h"
#include "quadrature.h"

BCManager::BCManager(Mesh3D &mesh, const char *filename, BlockCSRMatrix &matrix, std::vector<Block> &b) {
	mesh_ = &mesh;
	matrix_ = matrix;
	boundaryFaces_ = mesh_->getBoundaryFaces();
	b_ = b;
}

void BCManager::ApplyDirichle() {
	for (int i = 0; i < 8; i++) {
		for (Face &face : boundaryFaces_) {
			if (face.localFaceId != i) continue;
			for (int j = 0; j < 4; j++) {
				int node = face.nodes[j];
				matrix_.addBlock(node, node, Block(1, 0));
				b_[node].p_ = BoundaryFunctions::u_s(mesh_->getNodeCoord(node, 0), mesh_->getNodeCoord(node, 1), mesh_->getNodeCoord(node, 2));
				b_[node].c_ = BoundaryFunctions::u_c(mesh_->getNodeCoord(node, 0), mesh_->getNodeCoord(node, 1), mesh_->getNodeCoord(node, 2));
				std::vector<int> element = mesh_->getElementNodes(face.elementId);
				for (int k = 0; k < element.size(); k++) {
					matrix_.addBlock(node, k, Block(0, 0));
				}
			}
			
		}
	}
}

double BoundaryFunctions::u_s(double x, double y, double z) {
	return 0;
}

double BoundaryFunctions::u_c(double x, double y, double z) {
	return 0;
}