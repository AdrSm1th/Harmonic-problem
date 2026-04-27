//direct solver.cpp

#include <cmath>
#include <stdexcept>
#include <algorithm>
#include "direct solver.h"

inline double abs2(const Block &a) {
   return a.p_ * a.p_ + a.c_ * a.c_;
}

inline void sub_mul(BlockVector &target, const Block &coeff, const BlockVector &src) {
   target.p_ -= coeff.p_ * src.p_ - coeff.c_ * src.c_;
   target.c_ -= coeff.p_ * src.c_ + coeff.c_ * src.p_;
}

void GaussSolver::solve(BlockCSRMatrix &matrix, std::vector<BlockVector> &b, std::vector<BlockVector> &x) {
   int n = b.size();
   std::vector<Block> A(n * n);
   for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
         A[i * n + j] = matrix(i, j);
      }
   }

   for (int col = 0; col < n - 1; ++col) {
      Block diag = A[col * n + col];
      for (int i = col + 1; i < n; ++i) {
         Block factor = A[i * n + col] / diag;
         A[i * n + col] = factor;
         for (int j = col + 1; j < n; ++j) {
            Block prod = factor * A[col * n + j];
            A[i * n + j].p_ -= prod.p_;
            A[i * n + j].c_ -= prod.c_;
         }
         sub_mul(b[i], factor, b[col]);
      }
   }

   x.resize(n);
   for (int i = n - 1; i >= 0; --i) {
      BlockVector sum = b[i];
      for (int j = i + 1; j < n; ++j) {
         Block u = A[i * n + j];
         sub_mul(sum, u, x[j]);
      }
      Block diag = A[i * n + i];
      double den = diag.p_ * diag.p_ + diag.c_ * diag.c_;
      x[i].p_ = (sum.p_ * diag.p_ + sum.c_ * diag.c_) / den;
      x[i].c_ = (sum.c_ * diag.p_ - sum.p_ * diag.c_) / den;
   }
}