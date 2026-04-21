#include <iostream>

#include <petscsys.h>

#include "CoupledLV0DCoronary3D.h"

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);

  try
  {
    Examples::Heart::CoupledLV0DCoronary3D::Config cfg;
    cfg.meshPath = "CoronaryArtery.mesh";
    cfg.xdmfBasename = "CoronaryArtery";
    cfg.csvPath = "CoronaryArtery.csv";

    Examples::Heart::CoupledLV0DCoronary3D simulation(cfg);
    const int status = simulation.initialize().run();
    PetscFinalize();
    return status;
  }
  catch (const std::exception& e)
  {
    std::cerr << "Fatal error: " << e.what() << "\n";
    PetscFinalize();
    return 1;
  }
}
