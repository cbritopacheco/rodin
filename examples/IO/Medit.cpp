#include <Rodin/Geometry.h>

using namespace Rodin;
using namespace Geometry;

int main(int, char** argv)
{
  Mesh mesh;
  mesh.load(argv[1], IO::FileFormat::MEDIT);

  std::cout << "Saved Medit file to mfem.mesh in Mfem format" << std::endl;

  mesh.save("mfem.mesh", IO::FileFormat::MFEM);
};
