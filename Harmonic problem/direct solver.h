//direct solver.h

#pragma once

#include <vector>	
#include "sparse matrix.h"

class GaussSolver {
private:
public:
	void solve(BlockCSRMatrix &matrix, std::vector<BlockVector> &b, std::vector<BlockVector> &x);
};