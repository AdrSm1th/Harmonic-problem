//Mesh.cpp

#include <stdexcept>
#include <fstream>
#include "Mesh.h"

void Mesh3D::CalculateNonuniformDimension(std::vector<double> &dimension, double a, double b, double q, int n, char dimName) {
	if (q <= 0) throw std::invalid_argument("Invalid q" + dimName);
	double h1 = (b - a) * (1 - q) / (1 - std::pow(q, n - 1));
	dimension[0] = a;
	for (int i = 1; i < n; i++) {
		double h = h1 * std::pow(q, i - 1);
		dimension[i] = dimension[i - 1] + h;
	}
}

Mesh3D::Mesh3D(const char *filename)
{
	std::ifstream input(filename);
	double ax = 0, bx = 0, qx = 0;
	double ay = 0, by = 0, qy = 0;
	double az = 0, bz = 0, qz = 0;
	int nx = 0, ny = 0, nz = 0;
	bool uniformX, uniformY, uniformZ;
	input >> ax >> bx >> ay >> by >> az >> bz >> nx >> ny >> nz >> uniformX >> uniformY >> uniformZ >> qx >> qy >> qz;

	if (nx < 2) throw std::invalid_argument("Invalid nx");
	if (ny < 2) throw std::invalid_argument("Invalid ny");
	if (nz < 2) throw std::invalid_argument("Invalid nz");

	int nxe = nx - 1, nye = ny - 1, nze = nz - 1;
	x_.resize(nx);
	y_.resize(ny);
	z_.resize(nz);
	elements_.resize((nxe * nye * nze), std::vector<int>(8));
	boundaryFaces_.resize(2 * ((nxe * nye) + (nxe * nze) + (nye * nze)));

	double hx = (bx - ax) / nx;
	double hy = (by - ay) / ny;
	double hz = (bz - az) / nz;

	if (uniformX) {
		for (int i = 0; i < nx - 1; x_[i++] = ax + hx * i);
		x_[nx - 1] = bx;
	}
	else CalculateNonuniformDimension(x_, ax, bx, qx, nx, 'x');
	if (uniformY) {
		for (int i = 0; i < ny - 1; y_[i++] = ay + hy * i);
		y_[ny - 1] = by;
	}
	else CalculateNonuniformDimension(y_, ay, by, qy, ny, 'y');
	if (uniformZ) {
		for (int i = 0; i < nz - 1; z_[i++] = az + hz * i);
		z_[nz - 1] = bz;
	}
	else CalculateNonuniformDimension(z_, az, bz, qz, nz, 'z');

	int fidx = 0;
	for (int i = 0; i < nx - 1; i++) {
		for (int j = 0; j < ny - 1; j++) {
			for (int k = 0; k < nz - 1; k++) {
				int nxy = nx * ny;
				int idx = i + j * nxe + k * nxe * nye;

				elements_[idx][0] = idx;
				elements_[idx][1] = idx + 1;
				elements_[idx][2] = idx + nx;
				elements_[idx][3] = idx + nx + 1;
				elements_[idx][4] = idx + nxy;
				elements_[idx][5] = idx + nxy;
				elements_[idx][6] = idx + nxy + nx;
				elements_[idx][7] = idx + nxy + nx + 1;

				if (k == 0) {
					boundaryFaces_[fidx].nodes[0] = i;
					boundaryFaces_[fidx].nodes[1] = i + 1;
					boundaryFaces_[fidx].nodes[2] = i + nx;
					boundaryFaces_[fidx].nodes[3] = i + nx + 1;
					boundaryFaces_[fidx].elementId = idx;
					boundaryFaces_[fidx].localFaceId = 0;
					fidx++;
				}
				if (j == 0) {
					boundaryFaces_[fidx].nodes[0] = i;
					boundaryFaces_[fidx].nodes[1] = i + 1;
					boundaryFaces_[fidx].nodes[2] = i + nxy;
					boundaryFaces_[fidx].nodes[3] = i + nxy + 1;
					boundaryFaces_[fidx].elementId = idx;
					boundaryFaces_[fidx].localFaceId = 1;
					fidx++;
				}
				if (i == 0) {
					boundaryFaces_[fidx].nodes[0] = i;
					boundaryFaces_[fidx].nodes[1] = i + nx;
					boundaryFaces_[fidx].nodes[2] = i + nxy;
					boundaryFaces_[fidx].nodes[3] = i + nxy + 1;
					boundaryFaces_[fidx].elementId = idx;
					boundaryFaces_[fidx].localFaceId = 2;
					fidx++;
				}
				if (k + 1 == nze) {
					int s = (k + 1) * nx * ny;
					boundaryFaces_[fidx].nodes[0] = s;
					boundaryFaces_[fidx].nodes[1] = s + 1;
					int s2 = s + ny;
					boundaryFaces_[fidx].nodes[2] = s2;
					boundaryFaces_[fidx].nodes[3] = s2 + 1;
					boundaryFaces_[fidx].elementId = idx;
					boundaryFaces_[fidx].localFaceId = 3;
					fidx++;
				}
				if (j + 1 == nye) {
					int s = (j + 1) * nx;
					boundaryFaces_[fidx].nodes[0] = s;
					boundaryFaces_[fidx].nodes[1] = s + 1;
					int s2 = s + nx * ny;
					boundaryFaces_[fidx].nodes[2] = s2;
					boundaryFaces_[fidx].nodes[3] = s2 + 1;
					boundaryFaces_[fidx].elementId = idx;
					boundaryFaces_[fidx].localFaceId = 4;
					fidx++;
				}
				if (i + 1 == nxe) {
					int s = i + 1;
					boundaryFaces_[fidx].nodes[0] = s;
					boundaryFaces_[fidx].nodes[1] = s + nx;
					int s2 = s + nx * ny;
					boundaryFaces_[fidx].nodes[2] = s2;
					boundaryFaces_[fidx].nodes[3] = s2 + 1;
					boundaryFaces_[fidx].elementId = idx;
					boundaryFaces_[fidx].localFaceId = 5;
					fidx++;
				}
			}
		}
	}

	input.close();
}

int Mesh3D::getNumNodes() {
	return x_.size() * y_.size() * z_.size();
}

int Mesh3D::getNumElements() {
	return elements_.size();
}

double Mesh3D::getNodeCoord(int nodeId, int dimension) {
	switch (dimension)
	{
	case 0: return x_[nodeId % x_.size()];
	case 1: return y_[(nodeId / x_.size()) % y_.size()];
	case 2: return z_[nodeId / (x_.size() * y_.size())];
	default:
		throw std::invalid_argument("Dimension out of range");
	}
}

std::vector<int> Mesh3D::getElementNodes(int elementId) {
	return elements_[elementId];
}

double Mesh3D::getElementVolume(int elementId) {
	std::vector<int> element = elements_[elementId];
	double dx = x_[getNodeCoord(element[0], 0)] - x_[getNodeCoord(element[1], 0)];
	double dy = y_[getNodeCoord(element[0], 1)] - y_[getNodeCoord(element[2], 1)];
	double dz = z_[getNodeCoord(element[0], 2)] - z_[getNodeCoord(element[4], 2)];
	return dx * dy * dz;
}

std::vector<Face> Mesh3D::getBoundaryFaces()
{
	return boundaryFaces_;
}