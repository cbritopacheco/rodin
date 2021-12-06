#include <iostream>

#include <Rodin/Alert.h>
#include <Rodin/Variational.h>
#include <RodinExternal/MMG.h>

using namespace Rodin;
using namespace Rodin::External;

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
  MMG::Mesh2D mesh = MMG::Mesh2D::load(argv[1]);
  MMG::ScalarSolution2D sol =
    MMG::ScalarSolution2D::load(argv[2]).setMesh(mesh);

  // Convert the MMG mesh to Rodin
  auto rodinMesh = Cast<MMG::Mesh2D, Rodin::Mesh>().cast(mesh);

  // Build finite element space
  Variational::H1             Vh(rodinMesh);

  // Perform the cast using the finite element space we just built
  auto gf = Cast<
    MMG::ScalarSolution2D,
    Variational::GridFunction<Variational::H1>>(Vh).cast(sol);

  rodinMesh.save("multi-mat.mesh");
  gf.save("multi-mat.gf");

  return 0;
}
