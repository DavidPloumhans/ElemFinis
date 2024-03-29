#include <iostream>
#include <gmsh.h>

int main(int argc, char **argv)
{
  if(argc < 2) {
    std::cout << "Usage: " << argv[0] << " file.geo [options]" << std::endl;
    return 0;
  }

  gmsh::initialize();
  try{
    gmsh::open(argv[1]);
  }
  catch(...){
    exit(0);
  }
  gmsh::model::mesh::generate(3);
  gmsh::write("test.msh");
  gmsh::finalize();

  return 0;
}
