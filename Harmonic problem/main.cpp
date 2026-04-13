#include <iostream>
#include "mesh.h"
#include "sparse matrix.h"

int main()
{
   try {
      Mesh3D mesh("mesh.txt");
      BlockCSRMatrix matrix(mesh);
   }

   catch (const std::invalid_argument &e){
      std::cerr << "Error: " << e.what() << std::endl;
   }

   return 0;
}