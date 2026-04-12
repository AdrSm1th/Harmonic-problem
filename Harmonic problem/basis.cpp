//Basis.cpp

#include <tuple>
#include "basis.h"

Basis3D::Basis3D(double xi, double eta, double dzeta) {
	basis_[0] = (1 - xi) * (1 - eta) * (1 - dzeta) / 8;
	basis_[1] = (1 + xi) * (1 - eta) * (1 - dzeta) / 8;
	basis_[2] = (1 + xi) * (1 + eta) * (1 - dzeta) / 8;
	basis_[3] = (1 - xi) * (1 + eta) * (1 - dzeta) / 8;
	basis_[4] = (1 - xi) * (1 - eta) * (1 + dzeta) / 8;
	basis_[5] = (1 + xi) * (1 - eta) * (1 + dzeta) / 8;
	basis_[6] = (1 + xi) * (1 + eta) * (1 + dzeta) / 8;
	basis_[7] = (1 - xi) * (1 + eta) * (1 + dzeta) / 8;

	dpsi_dxi_[0] = -(1 - eta) * (1 - dzeta) / 8;
	dpsi_dxi_[1] = (1 - eta) * (1 - dzeta) / 8;
	dpsi_dxi_[2] = (1 + eta) * (1 - dzeta) / 8;
	dpsi_dxi_[3] = -(1 + eta) * (1 - dzeta) / 8;
	dpsi_dxi_[4] = -(1 - eta) * (1 + dzeta) / 8;
	dpsi_dxi_[5] = (1 - eta) * (1 + dzeta) / 8;
	dpsi_dxi_[6] = (1 + eta) * (1 + dzeta) / 8;
	dpsi_dxi_[7] = -(1 + eta) * (1 + dzeta) / 8;

	dpsi_deta_[0] = -(1 - xi) * (1 - dzeta) / 8;
	dpsi_deta_[1] = -(1 + xi) * (1 - dzeta) / 8;
	dpsi_deta_[2] = (1 + xi) * (1 - dzeta) / 8;
	dpsi_deta_[3] = (1 - xi) * (1 - dzeta) / 8;
	dpsi_deta_[4] = -(1 - xi) * (1 + dzeta) / 8;
	dpsi_deta_[5] = -(1 + xi) * (1 + dzeta) / 8;
	dpsi_deta_[6] = (1 + xi) * (1 + dzeta) / 8;
	dpsi_deta_[7] = (1 - xi) * (1 + dzeta) / 8;

	dpsi_dzeta_[0] = -(1 - xi) * (1 - eta) / 8;
	dpsi_dzeta_[1] = -(1 + xi) * (1 - eta) / 8;
	dpsi_dzeta_[2] = -(1 + xi) * (1 + eta) / 8;
	dpsi_dzeta_[3] = -(1 - xi) * (1 + eta) / 8;
	dpsi_dzeta_[4] = (1 - xi) * (1 - eta) / 8;
	dpsi_dzeta_[5] = (1 + xi) * (1 - eta) / 8;
	dpsi_dzeta_[6] = (1 + xi) * (1 + eta) / 8;
	dpsi_dzeta_[7] = (1 - xi) * (1 + eta) / 8;
}

std::array<double, 8>  Basis3D::getBasisFunctions() {
	return basis_;
}

std::tuple<std::array<double, 8>, std::array<double, 8>, std::array<double, 8>> Basis3D::getBasisDerivatives() {
	return std::make_tuple(dpsi_dxi_, dpsi_deta_, dpsi_dzeta_);
}