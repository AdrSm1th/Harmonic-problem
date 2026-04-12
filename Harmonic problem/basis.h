//Basis.h

#pragma once

#include <array>

class Basis3D {
private:
	std::array<double, 8> basis_;
	std::array<double, 8> dpsi_dxi_;
	std::array<double, 8> dpsi_deta_;
	std::array<double, 8> dpsi_dzeta_;

public:
	Basis3D(double xi, double eta, double dzeta);
	std::array<double, 8> getBasisFunctions();
	std::tuple<std::array<double, 8>, std::array<double, 8>, std::array<double, 8>> getBasisDerivatives();
};