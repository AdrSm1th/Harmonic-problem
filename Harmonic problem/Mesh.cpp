//Mesh.cpp

#include <stdexcept>
#include <fstream>
#include "mesh.h"

void Mesh3D::CalculateNonuniformDimension(std::vector<double> &dimension, double a, double b, double q, int n, const char *dimName) {
	if (q <= 0) throw std::invalid_argument(std::string("Invalid q") + dimName);
	double h1 = (b - a) * (1 - q) / (1 - std::pow(q, n - 1));
	dimension[0] = a;
	for (int i = 1; i < n; i++) {
		double h = h1 * std::pow(q, i - 1);
		dimension[i] = dimension[i - 1] + h;
	}
}

Mesh3D::Mesh3D() {

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

	double hx = (bx - ax) / (nx - 1);
	double hy = (by - ay) / (ny - 1);
	double hz = (bz - az) / (nz - 1);

	if (uniformX) {
		for (int i = 0; i < nx - 1; x_[i] = ax + hx * i++);
		x_[nx - 1] = bx;
	}
	else CalculateNonuniformDimension(x_, ax, bx, qx, nx, "x");
	if (uniformY) {
		for (int i = 0; i < ny - 1; y_[i] = ay + hy * i++);
		y_[ny - 1] = by;
	}
	else CalculateNonuniformDimension(y_, ay, by, qy, ny, "y");
	if (uniformZ) {
		for (int i = 0; i < nz - 1; z_[i] = az + hz * i++);
		z_[nz - 1] = bz;
	}
	else CalculateNonuniformDimension(z_, az, bz, qz, nz, "z");

	int fidx = 0;
	for (int k = 0; k < nz - 1; k++) {
		for (int j = 0; j < ny - 1; j++) {
			for (int i = 0; i < nx - 1; i++) {
				int nxy = nx * ny;
				int eidx = i + j * nxe + k * nxe * nye;
				int idx = i + j * nx + k * nx * ny;

				elements_[eidx][0] = idx;
				elements_[eidx][1] = idx + 1;
				elements_[eidx][2] = idx + nx;
				elements_[eidx][3] = idx + nx + 1;
				elements_[eidx][4] = idx + nxy;
				elements_[eidx][5] = idx + nxy + 1;
				elements_[eidx][6] = idx + nxy + nx;
				elements_[eidx][7] = idx + nxy + nx + 1;

				if (k == 0) {
					boundaryFaces_[fidx].nodes[0] = idx;
					boundaryFaces_[fidx].nodes[1] = idx + 1;
					boundaryFaces_[fidx].nodes[2] = idx + nx;
					boundaryFaces_[fidx].nodes[3] = idx + nx + 1;
					boundaryFaces_[fidx].elementId = eidx;
					boundaryFaces_[fidx].localFaceId = 0;
					fidx++;
				}
				if (j == 0) {
					boundaryFaces_[fidx].nodes[0] = idx;
					boundaryFaces_[fidx].nodes[1] = idx + nx;
					boundaryFaces_[fidx].nodes[2] = idx + nxy;
					boundaryFaces_[fidx].nodes[3] = idx + nxy + nx;
					boundaryFaces_[fidx].elementId = eidx;
					boundaryFaces_[fidx].localFaceId = 1;
					fidx++;
				}
				if (i == 0) {
					boundaryFaces_[fidx].nodes[0] = idx;
					boundaryFaces_[fidx].nodes[1] = idx + 1;
					boundaryFaces_[fidx].nodes[2] = idx + nxy;
					boundaryFaces_[fidx].nodes[3] = idx + nxy + 1;

					boundaryFaces_[fidx].elementId = eidx;
					boundaryFaces_[fidx].localFaceId = 2;
					fidx++;
				}
				if (k + 1 == nze) {
					boundaryFaces_[fidx].nodes[0] = idx + nxy;
					boundaryFaces_[fidx].nodes[1] = idx + nxy + 1;
					boundaryFaces_[fidx].nodes[2] = idx + nxy + nx;
					boundaryFaces_[fidx].nodes[3] = idx + nxy + nx + 1;
					boundaryFaces_[fidx].elementId = eidx;
					boundaryFaces_[fidx].localFaceId = 3;
					fidx++;
				}
				if (j + 1 == nye) {
					boundaryFaces_[fidx].nodes[0] = idx + nx;
					boundaryFaces_[fidx].nodes[1] = idx + nx + 1;
					boundaryFaces_[fidx].nodes[2] = idx + nxy + nx;
					boundaryFaces_[fidx].nodes[3] = idx + nxy + nx + 1;
					boundaryFaces_[fidx].elementId = eidx;
					boundaryFaces_[fidx].localFaceId = 4;
					fidx++;
				}
				if (i + 1 == nxe) {
					boundaryFaces_[fidx].nodes[0] = idx + 1;
					boundaryFaces_[fidx].nodes[1] = idx + nx + 1;
					boundaryFaces_[fidx].nodes[2] = idx + nxy + 1;
					boundaryFaces_[fidx].nodes[3] = idx + nxy + nx + 1;
					boundaryFaces_[fidx].elementId = eidx;
					boundaryFaces_[fidx].localFaceId = 5;
					fidx++;
				}
				//localFaceId:
				//0 - bottom
				//1 - backward
				//2 - left
				//3 - up
				//4 - forward
				//5 - rigth
			}
		}
	}

	input.close();
}

void Mesh3D::readCoefficients(const char *filename) {
	std::ifstream input(filename);
	input >> coefficients_.lambda >> coefficients_.omega >> coefficients_.chi >> coefficients_.sigma;
	input.close();
}

void Mesh3D::readBoundaryConditions(const char *filename) {
	std::ifstream input(filename);
	for (int i = 0; i < 8; i++)
	{
		input >> boundaryConditions_[i].type;

		switch (boundaryConditions_[i].type) {
		case 1:
			input >> boundaryConditions_[i].u_s >> boundaryConditions_[i].u_c;
			break;

		case 2:
			input >> boundaryConditions_[i].theta_s >> boundaryConditions_[i].theta_c;
			break;

		case 3:
			input >> boundaryConditions_[i].beta >> boundaryConditions_[i].ubeta_s >> boundaryConditions_[i].ubeta_c;
			break;
		}
	}
}

Coefficients Mesh3D::getCoefficients() const {
	return coefficients_;
}

int Mesh3D::getNumNodes() const {
	return (int)x_.size() * (int)y_.size() * (int)z_.size();
}

int Mesh3D::getNumElements() const {
	return (int)elements_.size();
}

double Mesh3D::getNodeCoord(int nodeId, int dimension) const {
	switch (dimension) {
	case 0: return x_[nodeId % x_.size()];
	case 1: return y_[(nodeId / x_.size()) % y_.size()];
	case 2: return z_[nodeId / (x_.size() * y_.size())];
	default:
		throw std::invalid_argument("Dimension out of range");
	}
}

std::vector<int> Mesh3D::getElementNodes(int elementId) const {
	return elements_[elementId];
}

double Mesh3D::getElementVolume(int elementId) const {
	std::vector<int> element = elements_[elementId];
	double dx = getNodeCoord(element[0], 0) - getNodeCoord(element[1], 0);
	double dy = getNodeCoord(element[0], 1) - getNodeCoord(element[2], 1);
	double dz = getNodeCoord(element[0], 2) - getNodeCoord(element[4], 2);
	return dx * dy * dz;
}

std::vector<Face> Mesh3D::getBoundaryFaces() const {
	return boundaryFaces_;
}