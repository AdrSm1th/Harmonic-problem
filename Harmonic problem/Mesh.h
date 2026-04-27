//Mesh.h

#pragma once

#include <vector>
#include <array>
#include <functional>

const double PI = 3.14159265358979323846;

struct Face {
	std::array<int, 4> nodes;
	int elementId;
	int localFaceId;
};

struct Coefficients {
	double lambda = 0;
	double sigma = 0;
	double chi = 0;
	double omega = 0;
};

class Mesh3D {
private:
	std::vector<double> x_;
	std::vector<double> y_;
	std::vector<double> z_;
	Coefficients coefficients_;
	std::vector<std::vector<int>> elements_;
	std::vector<Face> boundaryFaces_;
	void CalculateNonuniformDimension(std::vector<double> &dimension, double a, double b, double q, int n, const char *dimName);

public:
	Mesh3D();
	Mesh3D(const char *filename);
	void readCoefficients(const char *filename);
	Coefficients getCoefficients() const;
	int getNumNodes() const;
	int getNumElements() const;
	double getNodeCoord(int nodeId, int dimension) const;
	std::vector<int> getElementNodes(int elementId) const;
	std::vector<Face> getBoundaryFaces() const;

	double f_s(double x, double y, double z) {
		return (3 * coefficients_.lambda * PI * PI - coefficients_.omega * coefficients_.omega * coefficients_.chi) * (sin(PI * x) * sin(PI * y) * sin(PI * z)) -
			coefficients_.omega * coefficients_.sigma * (cos(PI * x) * cos(PI * y) * cos(PI * z));
		//return -(coefficients_.omega * coefficients_.sigma + coefficients_.omega * coefficients_.omega * coefficients_.chi) * (x + y + z);
		//return -6 * coefficients_.lambda - coefficients_.omega * coefficients_.sigma * (x * x - y * y) - coefficients_.omega * coefficients_.omega * coefficients_.chi * (x * x + y * y + z * z);
	}

	double f_c(double x, double y, double z) {
		return (3 * coefficients_.lambda * PI * PI - coefficients_.omega * coefficients_.omega * coefficients_.chi) * (cos(PI * x) * cos(PI * y) * cos(PI * z)) +
			coefficients_.omega * coefficients_.sigma * (sin(PI * x) * sin(PI * y) * sin(PI * z));
		//return (coefficients_.omega * coefficients_.sigma - coefficients_.omega * coefficients_.omega * coefficients_.chi) * (x + y + z);
		//return coefficients_.omega * coefficients_.sigma * (x * x + y * y + z * z) - coefficients_.omega * coefficients_.omega * coefficients_.chi * (x * x - y * y);
	}
};