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
    Rodin::Alert::Exception()
      << "Bad parameters.\nUsage:\n\t"
      << argv[0] << " <filename>.mesh" << " <filename>.gf"
      << Alert::Raise;
  }

  // Load Rodin mesh and solution
  auto rodinMesh = Mesh::load(argv[1]);
  auto gf = GridFunction<H1>::load(argv[2]);

  switch (rodinMesh.getDimension())
  {
    case 2:
    {
      switch (gf.getVectorDimension())
      {
        case 1:
        {
          // Convert the Rodin mesh to MMG
          auto mmgMesh = Cast(rodinMesh).to<MMG::Mesh2D>();

          // Perform the cast and set the mesh we just casted
          auto sol = Cast(gf).to<MMG::IncompleteScalarSolution2D>().setMesh(mmgMesh);

          // Save mesh and solution
          mmgMesh.save("mmg.mesh");
          sol.save("mmg.sol");
        }
      }
      break;
    }
    case 3:
    {
      switch (gf.getVectorDimension())
      {
        case 1:
        {
          // Convert the Rodin mesh to MMG
          auto mmgMesh = Cast(rodinMesh).to<MMG::Mesh3D>();

          // Perform the cast and set the mesh we just casted
          auto sol = Cast(gf).to<MMG::IncompleteScalarSolution3D>().setMesh(mmgMesh);

          // Save mesh and solution
          mmgMesh.save("mmg.mesh");
          sol.save("mmg.sol");
        }
      }
      break;
    }
    default:
      Alert::Exception("Bad mesh dimension").raise();
  }

  return 0;
}

