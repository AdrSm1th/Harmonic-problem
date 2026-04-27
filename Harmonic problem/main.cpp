#include <iostream>
#include <fstream>
#include <chrono>
#include "mesh.h"
#include "sparse matrix.h"
#include "assembler.h"
#include "boundary conditions.h"
#include "iteration solver.h"
#include "direct solver.h"

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

      std::chrono::steady_clock::time_point start, end;
      std::chrono::duration<double> duration;
      start = std::chrono::high_resolution_clock::now();
      //LOSsolver solver(1e-15, 1000);
      //solver.solve(matrix, b, x);
      GaussSolver solver;
      solver.solve(matrix, b, x);

      end = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration<double>(end - start);
      std::cout << "Runtime: " << duration.count() << " s\n";

      std::ofstream ofs("disc.txt");

      double maxdxs = 0, maxdxc = 0, normxs = 0, normxc = 0;
      for (int i = 0; i < mesh.getNumNodes(); i++) {
         double xc = mesh.getNodeCoord(i, 0), y = mesh.getNodeCoord(i, 1), z = mesh.getNodeCoord(i, 2);
         double us = BoundaryFunctions::u_s(xc, y, z), uc = BoundaryFunctions::u_c(xc, y, z);
         //ofs << x[i].p_ << " " << x[i].c_ << " " << us << " " << uc << " " << x[i].p_ - us << " " << x[i].c_ - uc << std::endl;
         ofs << x[i].p_ - us << " " << x[i].c_ - uc << std::endl;

         double dxs = abs(x[i].p_ - us);
         if (maxdxs < dxs) maxdxs = dxs;
         double dxc = abs(x[i].c_ - uc);
         if (maxdxc < dxc) maxdxc = dxc;

         normxs += dxs; normxc += dxc;
      }
      std::cout << maxdxs << " " << maxdxc << std::endl;
      std::cout << sqrt(normxs) / x.size() << " " << sqrt(normxc) / x.size() << std::endl;
   }

   catch (const std::invalid_argument &e){
      std::cerr << "Error: " << e.what() << std::endl;
   }

   return 0;
}