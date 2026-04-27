//boundary conditions.cpp

#include <fstream>
#include <set>
#include "boundary conditions.h"
#include "basis.h"
#include "quadrature.h"

BCManager::BCManager(Mesh3D &mesh) {
	mesh_ = mesh;
	boundaryFaces_ = mesh_.getBoundaryFaces();
}

void BCManager::ApplyDirichle(BlockCSRMatrix &matrix, std::vector<BlockVector> &b) {
	for (int i = 0; i < 6; i++) {
		for (Face &face : boundaryFaces_) {
			if (face.localFaceId != i) continue;
			for (int j = 0; j < 4; j++) {
				int node = face.nodes[j];
				matrix.changeBlock(node, node, Block(1, 0));
				b[node].p_ = BoundaryFunctions::u_s(mesh_.getNodeCoord(node, 0), mesh_.getNodeCoord(node, 1), mesh_.getNodeCoord(node, 2));
				b[node].c_ = BoundaryFunctions::u_c(mesh_.getNodeCoord(node, 0), mesh_.getNodeCoord(node, 1), mesh_.getNodeCoord(node, 2));
				std::vector<int> element = mesh_.getElementNodes(face.elementId);
				for (int k = 0; k < element.size(); k++) {
					if (node != element[k]) matrix.changeBlock(node, element[k], Block(0, 0));
				}
			}
		}
	}
}

double BoundaryFunctions::u_s(double x, double y, double z) {
	//return sin(PI * x) * sin(PI * y) * sin(PI * z);
	return x + y + z;
	//return 0;
	//return x * x + y * y + z * z;
}

double BoundaryFunctions::u_c(double x, double y, double z) {
	return 0;
	//return x + y + z;
	//return x * x - y * y;
}