//direct solver.h

#pragma once

#include <vector>	
#include <utility>
#include "sparse matrix.h"

class GaussSolver {
private:
public:
	void solve(BlockCSRMatrix &matrix, std::vector<BlockVector> &b, std::vector<BlockVector> &x);
};
