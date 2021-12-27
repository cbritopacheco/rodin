/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Mesh.h>
#include <Rodin/Solver.h>
#include <Rodin/Variational.h>
#include <RodinExternal/MMG.h>

using namespace Rodin;
using namespace Rodin::External;
using namespace Rodin::Variational;

template <class M, class L>
class Compliance
{
  public:
    Compliance(H1& vh,
        const ScalarCoefficient<M>& mu, const ScalarCoefficient<L>& lambda)
      : m_bf(vh), m_mu(mu), m_lambda(lambda)
    {};

    double operator()(GridFunction<H1>& v)
    {
      m_bf = ElasticityIntegrator(m_lambda, m_mu);
      return m_bf(v, v);
    }

  private:
    BilinearForm<H1>      m_bf;
    ScalarCoefficient<M>  m_mu;
    ScalarCoefficient<L>  m_lambda;
};

int main(int, char**)
{
  const char* meshFile = "../resources/mfem/meshes/holes.mesh";

  // Define boundary attributes
  int Gamma0 = 1, GammaD = 2, GammaN = 3;

  // Load mesh
  Mesh Omega = Mesh::load(meshFile);

  // Lamé coefficients
  auto mu     = ScalarCoefficient(0.3846),
       lambda = ScalarCoefficient(0.5769);

  auto ell = ScalarCoefficient(5);
  auto alpha = ScalarCoefficient(0.1);

  // Preconditioned Conjugate Gradient solver
  auto pcg = Solver::PCG().setMaxIterations(500)
                          .setRelativeTolerance(1e-12);

  size_t maxIt = 30;
  double eps = 1e-6;
  double coef = 0.1;
  double meshsize = 0.1;
  double oldObj, newObj;

  // Optimization loop
  for (size_t i = 0; i < maxIt; i++)
  {
    // Finite element spaces
    int d = 2;
    H1 Vh(Omega, d);
    H1 Sh(Omega);

    // Compliance
    Compliance compliance(Vh, mu, lambda);

    // Elasticity equation
    GridFunction u(Vh);
    Problem elasticity(u);
    elasticity = ElasticityIntegrator(lambda, mu)
               + DirichletBC(GammaD, VectorCoefficient{0, 0})
               + NeumannBC(GammaN, VectorCoefficient{0, -1});
    pcg.solve(elasticity);

    // Hilbert extension-regularization procedure
    GridFunction theta(Vh);
    auto e = ScalarCoefficient(0.5) * (Jacobian(u) + Jacobian(u).T());
    auto Ae = ScalarCoefficient(2.0) * mu * e + lambda * Trace(e) * IdentityMatrix(d);
    Problem hilbert(theta);
    hilbert = VectorDiffusionIntegrator(alpha)
            + VectorMassIntegrator()
            - VectorBoundaryFluxLFIntegrator({ Gamma0 }, Dot(Ae, e) - ell)
            + DirichletBC(GammaD, VectorCoefficient{0, 0})
            + DirichletBC(GammaN, VectorCoefficient{0, 0});
    pcg.solve(hilbert);

    // Update objective
    oldObj = newObj;
    newObj = compliance(u) + ell.value() * Omega.getVolume();

    std::cout << "[" << i << "] Objective: " << newObj << std::endl;

    // Test for convergence
    if (i > 0 && abs(oldObj - newObj) < eps)
      break;

    // Make the displacement
    GridFunction normGrad(Sh);
    normGrad = Magnitude(theta);
    double ngMax = normGrad.max();
    double step = coef * meshsize / (eps * eps + ngMax);
    theta *= step;
    double t = Omega.getMaximumDisplacement(theta);
    if (t < 1e-2)
    {
      // If the max displacement is too small reject the iteration
      coef /= 2;
      continue;
    }
    else
    {
      Omega.displace(theta);

      // Refine the mesh using MMG
      auto mmgMesh = Cast(Omega).to<MMG::Mesh2D>();
      MMG::MeshOptimizer2D().optimize(mmgMesh);
      Omega = Cast(mmgMesh).to<Mesh>();
    }
    coef = std::max(0.05, 1.02 * coef);

    Omega.save("Omega.mesh");
  }

  std::cout << "Saved final mesh to Omega.mesh" << std::endl;

  return 0;
}
