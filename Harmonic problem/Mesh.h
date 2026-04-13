//Mesh.h

#pragma once

#include <vector>
#include <array>

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

struct BoundaryCondition {
	int type;
	double u_s = 0;
	double u_c = 0;
	double theta_s = 0;
	double theta_c = 0;
	double beta = 0;
	double ubeta_s = 0;
	double ubeta_c = 0;
};

class Mesh3D {
private:
	std::vector<double> x_;
	std::vector<double> y_;
	std::vector<double> z_;
	Coefficients coefficients_;
	std::array<BoundaryCondition, 8> boundaryConditions_;
	std::vector<std::vector<int>> elements_;
	std::vector<Face> boundaryFaces_;
	void CalculateNonuniformDimension(std::vector<double> &dimension, double a, double b, double q, int n, const char *dimName);

public:
	Mesh3D();
	Mesh3D(const char *filename);
	void readCoefficients(const char *filename);
	Coefficients getCoefficients() const;
	void readBoundaryConditions(const char *filename);
	int getNumNodes() const;
	int getNumElements() const;
	double getNodeCoord(int nodeId, int dimension) const;
	std::vector<int> getElementNodes(int elementId) const;
	double getElementVolume(int elementId) const;
	std::vector<Face> getBoundaryFaces() const;

	double f_s(double x, double y, double z) {
		return 0;
	}

	double f_c(double x, double y, double z) {
		return 0;
	}
};