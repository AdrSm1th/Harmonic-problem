//iteration solver.cpp

#include <iostream>
#include <cmath>
#include <algorithm>
#include "iteration solver.h"

std::vector<BlockVector> LOSsolver::precondition(std::vector<Block> &M, std::vector<BlockVector> &vector) const {
	int n = (int)M.size();
	std::vector<BlockVector> result(n);

	for (int i = 0; i < n; i++) {
		double det = M[i].p_ * M[i].p_ + M[i].c_ * M[i].c_;
		if (abs(det) < 1e-12) throw std::runtime_error("Zero diagonal block determinant");
		result[i].p_ = (M[i].p_ * vector[i].p_ + M[i].c_ * vector[i].c_) / det;
		result[i].c_ = (-M[i].c_ * vector[i].p_ + M[i].p_ * vector[i].c_) / det;
	}

	return result;
}

double LOSsolver::dotProduct(std::vector<BlockVector> &v1, std::vector<BlockVector> &v2) const {
	double result = 0;
	int n = (int)v1.size();

	for (int i = 0; i < n; result += v1[i].p_ * v2[i].p_ + v1[i].c_ * v2[i].c_, i++);

	return result;
}

std::vector<BlockVector> LOSsolver::VplusV(std::vector<BlockVector> &v1, std::vector<BlockVector> &v2) const {
	int n = (int)v1.size();
	std::vector<BlockVector> result(n);

	for (int i = 0; i < n; result[i] = v1[i] + v2[i], i++);

	return result;
}

std::vector<BlockVector> LOSsolver::CmultV(double &c, std::vector<BlockVector> &v) const {
	int n = (int)v.size();
	std::vector<BlockVector> result(n);

	for (int i = 0; i < n; result[i] = v[i] * c, i++);

	return result;
}

void LOSsolver::solve(BlockCSRMatrix &matrix,
   std::vector<BlockVector> &b,
   std::vector<BlockVector> &x)
{
   int n = (int)b.size();

   std::vector<BlockVector> r(n), z(n), p(n);
   std::vector<Block> M(n);

   x.assign(n, BlockVector());

   for (int i = 0; i < n; i++)
   {
      M[i].p_ = matrix(i, i).p_;
      M[i].c_ = matrix(i, i).c_;
   }

   double bNorm = std::sqrt(dotProduct(b, b));
   if (bNorm < 1e-30)
      bNorm = 1.0;

   r = precondition(M, b);
   z = r;

   std::vector<BlockVector> t(n);
   matrix.multiply(z, t);
   p = precondition(M, t);

   double discrepancy = std::sqrt(dotProduct(r, r)) / bNorm;

   int iter = 0;

   while (discrepancy > eps_ && iter < maxiter_)
   {
      double pp = dotProduct(p, p);

      if (std::abs(pp) < 1e-30)
         break;

      double a = dotProduct(p, r) / pp;

      for (int i = 0; i < n; i++)
      {
         x[i].p_ += a * z[i].p_;
         x[i].c_ += a * z[i].c_;

         r[i].p_ -= a * p[i].p_;
         r[i].c_ -= a * p[i].c_;
      }

      discrepancy = std::sqrt(dotProduct(r, r)) / bNorm;

      if (discrepancy < eps_)
         break;

      std::vector<BlockVector> Mr = precondition(M, r);

      matrix.multiply(Mr, t);
      Mr = precondition(M, t);

      double betaDen = dotProduct(p, p);

      if (std::abs(betaDen) < 1e-30)
         break;

      double beta = -dotProduct(p, Mr) / betaDen;

      std::vector<BlockVector> MinvR = precondition(M, r);

      for (int i = 0; i < n; i++)
      {
         z[i].p_ = MinvR[i].p_ + beta * z[i].p_;
         z[i].c_ = MinvR[i].c_ + beta * z[i].c_;

         p[i].p_ = Mr[i].p_ + beta * p[i].p_;
         p[i].c_ = Mr[i].c_ + beta * p[i].c_;
      }

      iter++;
   }

   std::cout << "LOS iterations: " << iter
      << " residual: " << discrepancy << std::endl;
}

double BSGSTAB_LUsolver::dotProduct(const std::vector<BlockVector> &v1,
   const std::vector<BlockVector> &v2) const
{
   double result = 0.0;
   int n = (int)v1.size();

   for (int i = 0; i < n; i++)
   {
      result += v1[i].p_ * v2[i].p_ + v1[i].c_ * v2[i].c_;
   }

   return result;
}

double BSGSTAB_LUsolver::norm(const std::vector<BlockVector> &v) const
{
   return std::sqrt(dotProduct(v, v));
}

double BSGSTAB_LUsolver::relativeResidual(BlockCSRMatrix &matrix,
   const std::vector<BlockVector> &b,
   std::vector<BlockVector> &x) const
{
   int n = (int)b.size();

   std::vector<BlockVector> Ax(n);
   std::vector<BlockVector> r(n);

   matrix.multiply(x, Ax);

   for (int i = 0; i < n; i++)
   {
      r[i].p_ = b[i].p_ - Ax[i].p_;
      r[i].c_ = b[i].c_ - Ax[i].c_;
   }

   double bNorm = norm(b);
   if (bNorm < 1e-30)
      bNorm = 1.0;

   return norm(r) / bNorm;
}

void BSGSTAB_LUsolver::applyPreconditionedMatrix(BlockCSRMatrix &matrix,
   BlockCSRMatrix &lu,
   const std::vector<BlockVector> &v,
   std::vector<BlockVector> &result) const
{

   std::vector<BlockVector> tmp1;
   std::vector<BlockVector> tmp2;

   lu.solveProfileU(v, tmp1);
   matrix.multiply(tmp1, tmp2);
   lu.solveProfileL(tmp2, result);
}

void BSGSTAB_LUsolver::solve(BlockCSRMatrix &matrix,
   std::vector<BlockVector> &b,
   std::vector<BlockVector> &x)
{
   int n = (int)b.size();

   if ((int)x.size() != n)
      x.assign(n, BlockVector());
   else
      std::fill(x.begin(), x.end(), BlockVector());

   BlockCSRMatrix lu = matrix;
   lu.convertToProfile();
   lu.LUdecomposeProfile();

   std::vector<BlockVector> y;
   lu.multiplyProfileU(x, y);

   std::vector<BlockVector> Ax(n);
   std::vector<BlockVector> r_original(n);

   matrix.multiply(x, Ax);

   for (int i = 0; i < n; i++)
   {
      r_original[i].p_ = b[i].p_ - Ax[i].p_;
      r_original[i].c_ = b[i].c_ - Ax[i].c_;
   }

   std::vector<BlockVector> r;
   lu.solveProfileL(r_original, r);

   std::vector<BlockVector> r0 = r;
   std::vector<BlockVector> z = r;

   double discrepancy = relativeResidual(matrix, b, x);

   if (discrepancy < eps_)
   {
      std::cout << "BSGSTAB LU iterations: 0 residual: "
         << discrepancy << std::endl;
      return;
   }

   std::vector<BlockVector> Az(n);
   std::vector<BlockVector> p(n);
   std::vector<BlockVector> Ap(n);
   std::vector<BlockVector> r_new(n);
   std::vector<BlockVector> z_new(n);
   std::vector<BlockVector> x_current(n);

   int iter = 0;

   for (iter = 1; iter <= maxiter_; iter++)
   {
      applyPreconditionedMatrix(matrix, lu, z, Az);

      double rho_old = dotProduct(r, r0);
      double alpha_den = dotProduct(r0, Az);

      if (std::abs(alpha_den) < 1e-30)
         throw std::runtime_error("BSGSTAB LU breakdown: alpha denominator is zero");

      double alpha = rho_old / alpha_den;

      for (int i = 0; i < n; i++)
      {
         p[i].p_ = r[i].p_ - alpha * Az[i].p_;
         p[i].c_ = r[i].c_ - alpha * Az[i].c_;
      }

      std::vector<BlockVector> y_alpha = y;

      for (int i = 0; i < n; i++)
      {
         y_alpha[i].p_ += alpha * z[i].p_;
         y_alpha[i].c_ += alpha * z[i].c_;
      }

      lu.solveProfileU(y_alpha, x_current);

      discrepancy = relativeResidual(matrix, b, x_current);

      if (discrepancy < eps_)
      {
         x = x_current;

         std::cout << "BSGSTAB LU iterations: "
            << iter
            << " residual: "
            << discrepancy
            << std::endl;

         return;
      }

      applyPreconditionedMatrix(matrix, lu, p, Ap);

      double gamma_den = dotProduct(Ap, Ap);

      if (std::abs(gamma_den) < 1e-30)
         throw std::runtime_error("BSGSTAB LU breakdown: gamma denominator is zero");

      double gamma = dotProduct(p, Ap) / gamma_den;

      for (int i = 0; i < n; i++)
      {
         y[i].p_ += alpha * z[i].p_ + gamma * p[i].p_;
         y[i].c_ += alpha * z[i].c_ + gamma * p[i].c_;

         r_new[i].p_ = p[i].p_ - gamma * Ap[i].p_;
         r_new[i].c_ = p[i].c_ - gamma * Ap[i].c_;
      }

      lu.solveProfileU(y, x_current);

      discrepancy = relativeResidual(matrix, b, x_current);

      if (discrepancy < eps_)
      {
         x = x_current;

         std::cout << "BSGSTAB LU iterations: "
            << iter
            << " residual: "
            << discrepancy
            << std::endl;

         return;
      }

      double rho_new = dotProduct(r_new, r0);

      if (std::abs(gamma * rho_old) < 1e-30)
         throw std::runtime_error("BSGSTAB LU breakdown: beta denominator is zero");

      double beta = alpha * rho_new / (gamma * rho_old);

      for (int i = 0; i < n; i++)
      {
         z_new[i].p_ = r_new[i].p_
            + beta * z[i].p_
            - beta * gamma * Az[i].p_;

         z_new[i].c_ = r_new[i].c_
            + beta * z[i].c_
            - beta * gamma * Az[i].c_;
      }

      r = r_new;
      z = z_new;
   }

   lu.solveProfileU(y, x);

   discrepancy = relativeResidual(matrix, b, x);

   std::cout << "BSGSTAB LU iterations: "
      << maxiter_
      << " residual: "
      << discrepancy
      << std::endl;
}