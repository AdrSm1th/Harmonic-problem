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

void BlockCSRMatrix::addBlock(int row, int col, Block block) {
	(*this)(row, col) = (*this)(row, col) + block;
}

void BlockCSRMatrix::changeBlock(int row, int col, Block block) {
	(*this)(row, col) = block;
}

void BlockCSRMatrix::multiply(std::vector<BlockVector> &x, std::vector<BlockVector> &y) {
	int n = (int)di_.size();
	std::fill(y.begin(), y.end(), BlockVector());

	for (int i = 0; i < n; ++i) {
		y[i].p_ += di_[i].p_ * x[i].p_ - di_[i].c_ * x[i].c_;
		y[i].c_ += di_[i].c_ * x[i].p_ + di_[i].p_ * x[i].c_;

		for (int j_idx = ia_[i]; j_idx < ia_[i + 1]; ++j_idx) {
			int j = ja_[j_idx];
			if (i > j) {
				const Block &block = al_[j_idx];
				y[i].p_ += block.p_ * x[j].p_ - block.c_ * x[j].c_;
				y[i].c_ += block.c_ * x[j].p_ + block.p_ * x[j].c_;

				const Block &blocku = au_[j_idx];
				y[j].p_ += blocku.p_ * x[i].p_ - blocku.c_ * x[i].c_;
				y[j].c_ += blocku.c_ * x[i].p_ + blocku.p_ * x[i].c_;
			}
		}
	}
}

int BlockCSRMatrix::index(int i, int j) const {
	for (int jm = ia_[i]; jm < ia_[i + 1]; jm++) {
		if (ja_[jm] == j) {
			return jm;
		}
	}
	return -1;
}

int BlockCSRMatrix::getSize() {
	return di_.size();
}