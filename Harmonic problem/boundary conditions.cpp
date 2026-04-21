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
	for (int i = 0; i < 8; i++) {
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

//void BCManager::ApplyDirichle(BlockCSRMatrix &matrix, std::vector<BlockVector> &b) {
//   int n = mesh_.getNumNodes();
//   std::vector<bool> isDirichlet(n, false);
//   std::vector<BlockVector> u_boundary(n);
//   std::vector<std::set<int>> dirNeighbors(n);
//
//   // 1. Определяем дирихле-узлы и их значения
//   for (const Face &face : boundaryFaces_) {
//      for (int j = 0; j < 4; ++j) {
//         int node = face.nodes[j];
//         if (!isDirichlet[node]) {
//            isDirichlet[node] = true;
//            double us = BoundaryFunctions::u_s(mesh_.getNodeCoord(node, 0),
//               mesh_.getNodeCoord(node, 1),
//               mesh_.getNodeCoord(node, 2));
//            double uc = BoundaryFunctions::u_c(mesh_.getNodeCoord(node, 0),
//               mesh_.getNodeCoord(node, 1),
//               mesh_.getNodeCoord(node, 2));
//            u_boundary[node] = BlockVector(us, uc);
//         }
//      }
//   }
//
//   // 2. Строим соседей для дирихле-узлов через элементы граней
//   for (const Face &face : boundaryFaces_) {
//      std::vector<int> elemNodes = mesh_.getElementNodes(face.elementId);
//      for (int j = 0; j < 4; ++j) {
//         int node = face.nodes[j];
//         for (int other : elemNodes) {
//            if (other != node) {
//               dirNeighbors[node].insert(other);
//            }
//         }
//      }
//   }
//
//   // 3. Применяем условия
//   for (int node = 0; node < n; ++node) {
//      if (!isDirichlet[node]) continue;
//      BlockVector u_val = u_boundary[node];
//
//      for (int neighbor : dirNeighbors[node]) {
//         // Обнуляем строку node
//         matrix.changeBlock(node, neighbor, Block(0, 0));
//
//         // Получаем блок из столбца node
//         Block A_neighbor_node;
//         try {
//            A_neighbor_node = matrix(neighbor, node);
//         }
//         catch (const std::invalid_argument &) {
//            continue; // связь отсутствует (не должно случаться)
//         }
//
//         // Вычитаем вклад из правой части соседа
//         b[neighbor].p_ -= A_neighbor_node.p_ * u_val.p_ - A_neighbor_node.c_ * u_val.c_;
//         b[neighbor].c_ -= A_neighbor_node.c_ * u_val.p_ + A_neighbor_node.p_ * u_val.c_;
//
//         // Обнуляем столбец node
//         matrix.changeBlock(neighbor, node, Block(0, 0));
//      }
//
//      // Диагональ и правая часть для самого узла
//      matrix.changeBlock(node, node, Block(1, 0));
//      b[node] = u_val;
//   }
//}

double BoundaryFunctions::u_s(double x, double y, double z) {
	//return sin(PI * x) * sin(PI * y) * sin(PI * z);
	//return x + y + z;
	return x * x + y * y + z * z;
}

double BoundaryFunctions::u_c(double x, double y, double z) {
	//return 0;
	return x * x - y * y;
}