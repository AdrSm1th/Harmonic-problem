//quadrature.h

#pragma once

#include <vector>

class Quadrature {
public :
	static std::vector<double> getPoints(int order);
	static std::vector<double> getWeights(int order);
};