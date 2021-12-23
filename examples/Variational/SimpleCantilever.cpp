#include <Rodin/Mesh.h>
#include <Rodin/Solver.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Variational;

template <class M, class L>
class Compliance
{
  public:
    Compliance(Mesh& mesh,
        const ScalarCoefficient<M>& mu, const ScalarCoefficient<L>& lambda)
      : m_mesh(mesh), m_fes(m_mesh), m_d(m_mesh.getDimension()),
        m_mu(mu), m_lambda(lambda),
        m_one(m_fes)
    {
      m_one = ScalarCoefficient{1.0};
    };

    double operator()(GridFunction<H1>& v)
    {
      LinearForm lf(m_fes);
      auto e = ScalarCoefficient(0.5) * (Jacobian(v) + Jacobian(v).T());
      auto Ae = ScalarCoefficient(2.0) * m_mu * e + m_lambda * Trace(e) * IdentityMatrix(m_d);
      lf = DomainLFIntegrator(Dot(Ae, e));
      return lf(m_one);
    }

  private:
    Mesh&                 m_mesh;
    H1                    m_fes;
    int                   m_d;
    ScalarCoefficient<M>  m_mu;
    ScalarCoefficient<L>  m_lambda;
    GridFunction<H1>      m_one;
};

int main(int, char**)
{
  const char* meshFile = "../resources/mfem/meshes/holes.mesh";

  // Define boundary attributes
  int Gamma0 = 1, GammaD = 2, GammaN = 3;

  // Load mesh
  Mesh Omega = Mesh::load(meshFile);

  // Build finite element spaces
  int d = 2;
  H1 Vh(Omega, d);

  // Lamé coefficients
  auto mu     = ScalarCoefficient(0.3846),
       lambda = ScalarCoefficient(0.5769);

  // Compliance
  Compliance compliance(Omega, mu, lambda);

  // Elasticity equation
  GridFunction u(Vh);
  auto g = VectorCoefficient{0, -1};
  Problem elasticity(u);
  elasticity = ElasticityIntegrator(lambda, mu)
             + DirichletBC(GammaD, VectorCoefficient{0, 0})
             + NeumannBC(GammaN, g);

  // Hilbert extension-regularization
  GridFunction theta(Vh);
  auto alpha = ScalarCoefficient(0.001);
  auto e = ScalarCoefficient(0.5) * (Jacobian(u) + Jacobian(u).T());
  auto Ae = ScalarCoefficient(2.0) * mu * e + lambda * Trace(e) * IdentityMatrix(d);

  Problem hilbert(theta);
  hilbert = VectorDiffusionIntegrator(alpha)
          + VectorMassIntegrator()
          + VectorBoundaryFluxLFIntegrator(Dot(Ae, e))
          + DirichletBC(GammaD, VectorCoefficient{0, 0});

  // Solve both problems
  Solver::PCG().setMaxIterations(200)
               .setRelativeTolerance(1e-12)
               .solve(elasticity);

  Solver::PCG().setMaxIterations(200)
               .setRelativeTolerance(1e-12)
               .solve(hilbert);

  BilinearForm bf(Vh);
  bf = ElasticityIntegrator(lambda, mu);

  std::cout << compliance(u) << std::endl;
  std::cout << bf(u, u) << std::endl;

  Omega.save("Omega.mesh");
  u.save("u.gf");
  theta.save("theta.gf");

  return 0;
}
