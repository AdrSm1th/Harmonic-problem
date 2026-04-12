//sparse matrix.h

#pragma once

#include <vector>

struct Block
{
	double p = 0;
	double c = 0;
};

class BlockCSRMatrix {
private:
	std::vector<Block> data;
	std::vector<int> ia;
	std::vector<int> ja;

public:
	BlockCSRMatrix();
	void addBlock(int row, int col, Block block);
	void multiply(std::vector<double> &x, std::vector<double> &y);
	void convertToProfile();
	std::vector<Block> getLU();
};