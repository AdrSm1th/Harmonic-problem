#include <iostream>
#include "mesh.h"

int main()
{
   try {
      Mesh3D mesh("mesh.txt");
   }

   catch (const std::invalid_argument &e){
      std::cerr << "Error: " << e.what() << std::endl;
   }

   return 0;
}