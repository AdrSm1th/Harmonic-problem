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

	for (int i = 1; i < (int)neighbours.size(); i++) {
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
	return (int)di_.size();
}

void BlockCSRMatrix::convertToProfile() {
	int n = (int)di_.size();

	firstCol_.resize(n);
	ia_profile_.resize(n + 1);
	ia_profile_[0] = 0;

	for (int i = 0; i < n; i++) {
		int min_col = i;

		for (int k = ia_[i]; k < ia_[i + 1]; k++) {
			int j = ja_[k];
			if (j < min_col) {
				min_col = j;
			}
		}

		firstCol_[i] = min_col;
	}

	for (int i = 0; i < n; i++) {
		int width = i - firstCol_[i];
		ia_profile_[i + 1] = ia_profile_[i] + width;
	}

	int totalSize = ia_profile_[n];
	al_profile_.assign(totalSize, Block());
	au_profile_.assign(totalSize, Block());

	for (int i = 0; i < n; i++) {
		int start = ia_profile_[i];
		int fc = firstCol_[i];

		for (int j = fc; j < i; j++) {
			int idx = start + (j - fc);

			int csr_idx = index(i, j);
			if (csr_idx != -1) {
				al_profile_[idx] = al_[csr_idx];
				au_profile_[idx] = au_[csr_idx];
			}
			else {
				al_profile_[idx] = Block(0, 0);
				au_profile_[idx] = Block(0, 0);
			}
		}
	}
}

void BlockCSRMatrix::LUdecomposeProfile() {
	int n = (int)di_.size();

	for (int i = 0; i < n; i++) {
		int fc_i = firstCol_[i];

		for (int j = fc_i; j < i; j++) {
			int idx_ij = ia_profile_[i] + (j - fc_i);
			Block sum = al_profile_[idx_ij];

			int fc_j = firstCol_[j];
			int k_start = std::max(fc_i, fc_j);

			for (int k = k_start; k < j; k++) {
				int idx_ik = ia_profile_[i] + (k - fc_i);
				int idx_kj = ia_profile_[j] + (k - fc_j);
				sum -= al_profile_[idx_ik] * au_profile_[idx_kj];
			}

			al_profile_[idx_ij] = sum / di_[j];
		}

		for (int j = fc_i; j < i; j++) {
			int idx_ij = ia_profile_[i] + (j - fc_i);
			Block sum = au_profile_[idx_ij];

			int fc_j = firstCol_[j];
			int k_start = std::max(fc_i, fc_j);

			for (int k = k_start; k < j; k++) {
				int idx_jk = ia_profile_[j] + (k - fc_j);
				int idx_ki = ia_profile_[i] + (k - fc_i);
				sum -= al_profile_[idx_jk] * au_profile_[idx_ki];
			}

			au_profile_[idx_ij] = sum;
		}

		Block sum = di_[i];

		for (int k = fc_i; k < i; k++) {
			int idx_ik = ia_profile_[i] + (k - fc_i);
			sum -= al_profile_[idx_ik] * au_profile_[idx_ik];
		}

		di_[i] = sum;
	}
}

void BlockCSRMatrix::solveProfileLU(std::vector<BlockVector> &b, std::vector<BlockVector> &x) {
	int n = (int)di_.size();

	std::vector<BlockVector> y(n);

	for (int i = 0; i < n; i++) {
		BlockVector sum = b[i];
		int fc_i = firstCol_[i];

		for (int j = fc_i; j < i; j++) {
			int idx = ia_profile_[i] + (j - fc_i);
			Block lij = al_profile_[idx];

			sum.p_ -= lij.p_ * y[j].p_ - lij.c_ * y[j].c_;
			sum.c_ -= lij.p_ * y[j].c_ + lij.c_ * y[j].p_;
		}

		y[i] = sum;
	}

	x.resize(n);

	for (int i = n - 1; i >= 0; i--) {
		BlockVector sum = y[i];

		for (int j = i + 1; j < n; j++) {
			int fc_j = firstCol_[j];

			if (i >= fc_j) {
				int idx = ia_profile_[j] + (i - fc_j);
				Block uij = au_profile_[idx];

				sum.p_ -= uij.p_ * x[j].p_ - uij.c_ * x[j].c_;
				sum.c_ -= uij.p_ * x[j].c_ + uij.c_ * x[j].p_;
			}
		}

		Block diag = di_[i];
		double den = diag.p_ * diag.p_ + diag.c_ * diag.c_;

		if (den < 1e-14)
			throw std::runtime_error("Zero diagonal in back substitution");

		x[i].p_ = (sum.p_ * diag.p_ + sum.c_ * diag.c_) / den;
		x[i].c_ = (sum.c_ * diag.p_ - sum.p_ * diag.c_) / den;
	}
}