#include <iostream>

#include <Rodin/Alert.h>
#include <Rodin/Variational.h>
#include <RodinExternal/MMG.h>

using namespace Rodin;
using namespace Rodin::External;
using namespace Rodin::Variational;

int main(int argc, char** argv)
{
  if (argc != 3)
  {
    auto ex = Rodin::Alert::Exception()
      << "Bad parameters.\nUsage:\n\t"
      << argv[0] << " <filename>.mesh" << " <filename>.sol";
    ex.raise();
  }

  // Load Rodin mesh and solution
  auto rodinMesh = Mesh::load(argv[1]);
  auto gf = GridFunction<H1>::load(argv[2]);

  // Convert the Rodin mesh to MMG
  auto mmgMesh = Cast(rodinMesh).to<MMG::Mesh2D>();

  // Perform the cast and set the mesh we just casted
  auto sol = Cast(gf).to<MMG::IncompleteScalarSolution2D>().setMesh(mmgMesh);

  // Save mesh and solution
  mmgMesh.save("multi-mat.mesh");
  sol.save("multi-mat.sol");

  return 0;
}

