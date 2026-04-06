//Mesh.h

#pragma once

#include <vector>
#include <array>

struct Face {
	std::array<int, 4> nodes;
	int elementId;
	int localFaceId;
};

class Mesh3D {
private:
	std::vector<double> x_;
	std::vector<double> y_;
	std::vector<double> z_;
	std::vector<std::vector<int>> elements_;
	std::vector<Face> boundaryFaces_;

	void CalculateNonuniformDimension(std::vector<double> &dimension, double a, double b, double q, int n, const char *dimName);
public:
	Mesh3D(const char *filename);
	int getNumNodes();
	int getNumElements();
	double getNodeCoord(int nodeId, int dimension);
	std::vector<int> getElementNodes(int elementId);
	double getElementVolume(int elementId);
	std::vector<Face> getBoundaryFaces();
};