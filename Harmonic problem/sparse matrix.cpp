//sparse matrix.cpp

#include <set>
#include <stdexcept>
#include "sparse matrix.h"

const int ELEM_SIZE = 8;

BlockCSRMatrix::BlockCSRMatrix(Mesh3D &mesh) {
	int n = mesh.getNumNodes();
	ia_.resize(n + 1);
	di_.resize(n);
	std::vector<std::set<int>> neighbours(mesh.getNumNodes());
	for (int k = 0; k < mesh.getNumElements(); k++) {
		std::vector<int> elem = mesh.getElementNodes(k);

		for (int i = 0; i < ELEM_SIZE; i++) {
			for (int j = 0; j < i; j++) {
				neighbours[elem[i]].insert(elem[j]);
			}
		}
	}

	for (int i = 1; i < neighbours.size(); i++) {
		ia_[i + 1] = ia_[i] + (int)neighbours[i].size();
		for (int neighbour : neighbours[i]) {
			ja_.push_back(neighbour);
		}
	}
	al_.resize(ia_[ia_.size() - 1]);
	au_.resize(ia_[ia_.size() - 1]);
}

void BlockCSRMatrix::addBlock(int row, int col, Block &block) {
	(*this)(row, col) = block;
}

void BlockCSRMatrix::multiply(std::vector<double> &x, std::vector<double> &y) {
	int n = di_.size();
	std::fill(y.begin(), y.end(), 0.0);

	for (int i = 0; i < n; ++i) {
		y[2 * i] += di_[i].p_ * x[2 * i] + di_[i].c_ * x[2 * i + 1];
		y[2 * i + 1] += di_[i].c_ * x[2 * i] + di_[i].p_ * x[2 * i + 1];

		for (int j_idx = ia_[i]; j_idx < ia_[i + 1]; ++j_idx) {
			int j = ja_[j_idx];
			if (i > j) {
				const Block &block = al_[j_idx];
				y[2 * i] += block.p_ * x[2 * j] + block.c_ * x[2 * j + 1];
				y[2 * i + 1] += block.c_ * x[2 * j] + block.p_ * x[2 * j + 1];
			}

			else if (i < j) {
				const Block &block = au_[j_idx];
				y[2 * i] += block.p_ * x[2 * j] + block.c_ * x[2 * j + 1];
				y[2 * i + 1] += block.c_ * x[2 * j] + block.p_ * x[2 * j + 1];
			}
		}
	}
}

void BlockCSRMatrix::convertToProfile() {
	std::vector<Block> newAl;
	std::vector<int> newIa;

	int jm = 0;
	int ia = 0;
	for (int i = 0; i < di_.size(); i++) {
		int end = jm + i;
		newIa.push_back(ia);
		for (int j = ia_[i]; jm < end; jm++, j++) {
			ia++;
			if (j == ja_[jm]) {
				newAl.push_back(al_[jm]);
			}
			else {
				Block nullBlock;
				nullBlock.p_ = 0;
				nullBlock.c_ = 0;
				newAl.push_back(nullBlock);
			}
		}
	}
	newIa.push_back(ia);

	al_ = newAl;
	ia_ = newIa;
}

int BlockCSRMatrix::index(int i, int j) const {
	for (int jm = ia_[i]; jm < ia_[i + 1]; jm++) {
		if (ja_[jm] == j) {
			return jm;
		}
	}
	return -1;
}

void BlockCSRMatrix::getLU() {
	int n = di_.size();
	for (int i = 1; i < n; i++) {
		for (int k = ia_[i]; k < ia_[i + 1]; k++) {
			//al
			al_[k] /= di_[ja_[k]];

			for (int j = ja_[k]; j < n; j++) {
				//di
				if (j == i) {
					for (int idx = ia_[j]; idx < ia_[j + 1]; idx++) {
						if (ja_[idx] == k) {
							di_[i] -= al_[k] * au_[idx];
						}
					}
				}

				//au
				for (int idx = ia_[j]; idx < ia_[j]; idx++) {
					if (ja_[idx] == i) {
						for (int idx_2 = ia_[j]; idx_2 < ia_[j + 1]; idx_2++) {
							if (ja_[idx_2] == k) {
								au_[idx] -= al_[k] * au_[idx_2];
							}
						}
					}
				}
			}
		}
	}
}