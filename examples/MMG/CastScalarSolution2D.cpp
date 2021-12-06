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

  // Load MMG mesh and solution
  auto mmgMesh = MMG::Mesh2D::load(argv[1]);
  auto mmgSol  = MMG::ScalarSolution2D::load(argv[2]).setMesh(mmgMesh);

  // Convert the MMG mesh to Rodin
  auto rodinMesh = Cast<MMG::Mesh2D, Rodin::Mesh>().cast(mmgMesh);

  // Build finite element space
  H1 Vh(rodinMesh);

  // Perform the cast using the finite element space we just built
  auto gf = Cast<MMG::ScalarSolution2D, GridFunction<H1>>(Vh).cast(mmgSol);

  rodinMesh.save("multi-mat.mesh");
  gf.save("multi-mat.gf");

  return 0;
}
