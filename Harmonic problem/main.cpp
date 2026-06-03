#include <iostream>
#include <fstream>
#include <chrono>
#include "mesh.h"
#include "sparse matrix.h"
#include "assembler.h"
#include "boundary conditions.h"
#include "iteration solver.h"
#include "direct solver.h"

struct ErrorInfo {
   double max_s = 0.0;
   double max_c = 0.0;

   double l2_s = 0.0;
   double l2_c = 0.0;

   double rel_l2_s = 0.0;
   double rel_l2_c = 0.0;
};

ErrorInfo calculateError(Mesh3D &mesh,
   const std::vector<BlockVector> &x,
   const char *filename = "output.txt")
{
   int n = mesh.getNumNodes();

   if ((int)x.size() != n)
      throw std::runtime_error("Wrong solution vector size in calculateError");

   std::ofstream ofs(filename);

   ofs << "x y z "
      << "u_s_num u_c_num "
      << "u_s_exact u_c_exact "
      << "err_s err_c"
      << std::endl;

   ErrorInfo error;

   double sumErrS2 = 0.0;
   double sumErrC2 = 0.0;

   double sumExactS2 = 0.0;
   double sumExactC2 = 0.0;

   for (int i = 0; i < n; i++)
   {
      double xc = mesh.getNodeCoord(i, 0);
      double y = mesh.getNodeCoord(i, 1);
      double z = mesh.getNodeCoord(i, 2);

      double us = BoundaryFunctions::u_s(xc, y, z);
      double uc = BoundaryFunctions::u_c(xc, y, z);

      double err_s = x[i].p_ - us;
      double err_c = x[i].c_ - uc;

      ofs << xc << " "
         << y << " "
         << z << " "
         << x[i].p_ << " "
         << x[i].c_ << " "
         << us << " "
         << uc << " "
         << err_s << " "
         << err_c << std::endl;

      error.max_s = std::max(error.max_s, std::abs(err_s));
      error.max_c = std::max(error.max_c, std::abs(err_c));

      sumErrS2 += err_s * err_s;
      sumErrC2 += err_c * err_c;

      sumExactS2 += us * us;
      sumExactC2 += uc * uc;
   }

   error.l2_s = std::sqrt(sumErrS2 / n);
   error.l2_c = std::sqrt(sumErrC2 / n);

   if (sumExactS2 > 1e-30)
      error.rel_l2_s = std::sqrt(sumErrS2 / sumExactS2);

   if (sumExactC2 > 1e-30)
      error.rel_l2_c = std::sqrt(sumErrC2 / sumExactC2);

   return error;
}

int main()
{
   try {
      Mesh3D mesh("mesh.txt");
      mesh.readCoefficients("coefs.txt");
      BlockCSRMatrix matrix(mesh);
      std::vector<BlockVector> b(mesh.getNumNodes());
      HarmonicAssembler assembler(mesh);
      assembler.assembleSystem(matrix, b);
      BCManager bcmanager(mesh);
      bcmanager.ApplyDirichle(matrix, b);
      std::vector<BlockVector> x(mesh.getNumNodes());

      LOSsolver solverlos(1e-12, 3000);
      GaussSolver solvergauss;
      BSGSTAB_LUsolver solverbsgstab(1e-12, 3000);

      int choice = 0;
      std::cin >> choice;

      std::chrono::steady_clock::time_point start, end;
      std::chrono::duration<double> duration;
      start = std::chrono::high_resolution_clock::now();

      if (choice == 1) solverlos.solve(matrix, b, x);
      else if (choice == 2) solvergauss.solve(matrix, b, x);
      else if (choice == 3) {
         matrix.convertToProfile();
         matrix.LUdecomposeProfile();
         matrix.solveProfileLU(b, x);
      }
      else if (choice == 4) {
         solverbsgstab.solve(matrix, b, x);
      }

      end = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration<double>(end - start);
      std::cout << "Runtime: " << duration.count() << " s\n";

      ErrorInfo error = calculateError(mesh, x, "disc.txt");

      std::cout << "Max error u_s: " << error.max_s << std::endl;
      std::cout << "Max error u_c: " << error.max_c << std::endl;

      std::cout << "L2 error u_s: " << error.l2_s << std::endl;
      std::cout << "L2 error u_c: " << error.l2_c << std::endl;

      std::cout << "Relative L2 error u_s: " << error.rel_l2_s << std::endl;
      std::cout << "Relative L2 error u_c: " << error.rel_l2_c << std::endl;
   }

   catch (const std::invalid_argument &e){
      std::cerr << "Error: " << e.what() << std::endl;
   }

   return 0;
}