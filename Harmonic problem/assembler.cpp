//assembler.cpp

#include "assembler.h"
#include "basis.h"

HarmonicAssembler::HarmonicAssembler(Mesh3D &mesh) {
	mesh_ = mesh;
	coefficients_ = mesh_.getCoefficients();
}

void HarmonicAssembler::computeLocalMatrices(int elemId) {
	std::fill(A_local_.begin(), A_local_.end(), Block());
	std::fill(b_local_.begin(), b_local_.end(), BlockVector());

	std::vector<int> element = mesh_.getElementNodes(elemId);
	double x1 = mesh_.getNodeCoord(element[0], 0);
	double x2 = mesh_.getNodeCoord(element[1], 0);
	double y1 = mesh_.getNodeCoord(element[0], 1);
	double y2 = mesh_.getNodeCoord(element[2], 1);
	double z1 = mesh_.getNodeCoord(element[0], 2);
	double z2 = mesh_.getNodeCoord(element[4], 2);
	double detJ = (x2 - x1) * (y2 - y1) * (z2 - z1) / 8;


	std::vector<double> points = quadrature_.getPoints(2);

	for (int ixi = 0; ixi < 2; ixi++)
	{
		//double xi = 2 * (points[ixi] - x1) / (x2 - x1) / 2;
		double xi = points[ixi];
		for (int ieta = 0; ieta < 2; ieta++)
		{
			//double eta = 2 * (points[ieta] - y1) / (y2 - y1) / 2;
			double eta = points[ieta];
			for (int izeta = 0; izeta < 2; izeta++)
			{
				//double zeta = 2 * (points[izeta] - z1) / (z2 - z1) / 2;
				double zeta = points[izeta];
				Basis3D bs(xi, eta, zeta);
				std::array<double, 8> psi = bs.getBasisFunctions();
				auto [dpsi_dxi, dpsi_deta, dpsi_dzeta] = bs.getBasisDerivatives();

				double h = x2 - x1;
				for (double &e : dpsi_dxi) {
					e *= 2 / h;
				}

				h = y2 - y1;
				for (double &e : dpsi_deta) {
					e *= 2 / h;
				}

				h = z2 - z1;
				for (double &e : dpsi_dzeta) {
					e *= 2 / h;
				}

				for (int i = 0; i < 8; i++) {
					for (int j = 0; j < 8; j++) {
						double dot = dpsi_dxi[i] * dpsi_dxi[j] + dpsi_deta[i] * dpsi_deta[j] + dpsi_dzeta[i] * dpsi_dzeta[j];
						A_local_[i * 8 + j].p_ += (coefficients_.lambda * dot - coefficients_.omega * coefficients_.omega * coefficients_.chi * psi[i] * psi[j]) * detJ;
						A_local_[i * 8 + j].c_ += coefficients_.sigma * coefficients_.omega * psi[i] * psi[j] * detJ;
					}
					double x_global = 0, y_global = 0, z_global = 0;
					for (int k = 0; k < 8; ++k) {
						x_global += psi[k] * mesh_.getNodeCoord(element[k], 0);
						y_global += psi[k] * mesh_.getNodeCoord(element[k], 1);
						z_global += psi[k] * mesh_.getNodeCoord(element[k], 2);
					}
					b_local_[i].p_ += mesh_.f_s(x_global, y_global, z_global) * psi[i] * detJ;
					b_local_[i].c_ += mesh_.f_c(x_global, y_global, z_global) * psi[i] * detJ;
				}
			}
		}
	}
}

void HarmonicAssembler::assembleSystem(BlockCSRMatrix &matrix, std::vector<BlockVector> &global_b) {
	for (int k = 0; k < mesh_.getNumElements(); k++) {
		std::vector<int> element = mesh_.getElementNodes(k);
		computeLocalMatrices(k);
		for (int i = 0; i < 8; i++) {
			int global_i = element[i];
			global_b[global_i].p_ += b_local_[i].p_;
			global_b[global_i].c_ += b_local_[i].c_;

			for (int j = 0; j < 8; j++) {
				int global_j = element[j];
				matrix.addBlock(global_i, global_j, A_local_[i * 8 + j]);
			}
		}
	}
}