/*
 * Installation Test - Poisson Equation Solver
 * 
 * This test verifies that the installed Rodin library can be:
 * 1. Found by CMake using find_package()
 * 2. Linked against using CMake targets
 * 3. Used to solve a simple Poisson equation
 */

#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Assembly.h>
#include <Rodin/Variational.h>
#include <iostream>

using namespace Rodin;
using namespace Rodin::Solver;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main()
{
  std::cout << "Starting Rodin installation test..." << std::endl;
  
  try {
    // Create a 2D triangular mesh
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, {16, 16});
    mesh.getConnectivity().compute(1, 2);

    // Define P1 finite element space
    P1 Vh(mesh);
    
    // Define trial and test functions
    TrialFunction u(Vh);
    TestFunction v(Vh);
    
    // Right-hand side function f = 1
    RealFunction f = 1.0;
    
    // Assemble the Poisson problem:
    // -Δu = f in Ω
    // u = 0 on ∂Ω
    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, Zero());
    
    // Solve the system using CG solver
    CG(poisson).solve();
    
    // Print results
    std::cout << "========================================" << std::endl;
    std::cout << "✓ Poisson equation solved successfully!" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Mesh cells: " << mesh.getCellCount() << std::endl;
    std::cout << "DOFs: " << Vh.getSize() << std::endl;
    std::cout << "========================================" << std::endl;
    
    // Verify we got reasonable results
    if (mesh.getCellCount() == 0 || Vh.getSize() == 0) {
      std::cerr << "✗ ERROR: Invalid results" << std::endl;
      return 1;
    }
    
    std::cout << "✓ Installation test PASSED" << std::endl;
    return 0;
    
  } catch (const std::exception& e) {
    std::cerr << "✗ ERROR: " << e.what() << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "✗ ERROR: Unknown exception" << std::endl;
    return 1;
  }
}
