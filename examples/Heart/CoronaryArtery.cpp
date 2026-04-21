/**
 * @file CoronaryArtery.cpp
 * @brief Driver for the coupled LV-0D / coronary 3D Newtonian flow example.
 *
 * This executable configures and runs `CoupledLV0DCoronary3D` with default
 * paths for the example mesh and output files:
 * - Mesh: `../resources/examples/Heart/CoronaryArtery_Fluid.medit.mesh`
 * - XDMF basename: `CoronaryArtery`
 * - CSV output: `CoronaryArtery.csv`
 */
#include <exception>
#include <iostream>

#include <petscsys.h>

#include "CoronaryArtery/CoupledLV0DCoronary3D.h"

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);

  try
  {
    Rodin::Examples::Heart::CoupledLV0DCoronary3D::Config cfg;
    cfg.meshPath = "../resources/examples/Heart/CoronaryArtery_Fluid.medit.mesh";
    cfg.xdmfBasename = "CoronaryArtery";
    cfg.csvPath = "CoronaryArtery.csv";

    Rodin::Examples::Heart::CoupledLV0DCoronary3D simulation(cfg);
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
