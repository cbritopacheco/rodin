#include <gtest/gtest.h>
#include <random>
#include "Rodin/QF/OrthogonalBasis.h"
#include "Rodin/QF/NodeElimination.h"
using namespace Rodin; using namespace Rodin::QF; using namespace Rodin::Geometry;

TEST(OrthogonalBasisTest, IsOrthogonalOnTheElement)
{
  for (const auto g : {Polytope::Type::Triangle, Polytope::Type::Tetrahedron})
  {
    const size_t d = Polytope::Traits(g).getDimension();
    for (size_t p : {4u, 10u, 16u, 20u})
    {
      if (d == 3 && p > 10) continue;
      const auto idx = OrthogonalBasis::indices(d, p);
      const auto rule = NodeElimination::productSeed(g, 2 * p + 2);
      const Eigen::Index n = (Eigen::Index)idx.size(), nq = (Eigen::Index)rule.getSize();
      Math::Matrix<Real> V(n, nq);
      for (Eigen::Index i = 0; i < n; ++i)
        for (Eigen::Index q = 0; q < nq; ++q)
          V(i, q) = OrthogonalBasis::evaluate(idx[(size_t)i], rule.getPoint((size_t)q))
                    * std::sqrt(rule.getWeight((size_t)q));
      const Math::Matrix<Real> G = V * V.transpose();
      Real off = 0;
      for (Eigen::Index i = 0; i < n; ++i)
        for (Eigen::Index j = 0; j < n; ++j)
          if (i != j) off = std::max(off, std::abs(G(i,j)) / std::sqrt(G(i,i)*G(j,j)));
      EXPECT_LT(off, 1e-12) << (d == 2 ? "tri" : "tet") << " degree " << p;
    }
  }
}

/// @brief FD-consistency gate on the basis gradient, as the house pattern
/// requires for a residual and tangent pair.
TEST(OrthogonalBasisTest, GradientMatchesFiniteDifferences)
{
  std::mt19937 rng(11u);
  std::uniform_real_distribution<Real> uni(0.05, 0.9);
  for (const auto g : {Polytope::Type::Triangle, Polytope::Type::Tetrahedron})
  {
    const size_t d = Polytope::Traits(g).getDimension();
    const auto idx = OrthogonalBasis::indices(d, 6);
    Real worst = 0;
    for (int trial = 0; trial < 40; ++trial)
    {
      Math::SpatialVector<Real> x;
      x.resize((Eigen::Index)d);
      Real tot;
      do {
        tot = 0;
        for (size_t k = 0; k < d; ++k) { x[(Eigen::Index)k] = uni(rng); tot += x[(Eigen::Index)k]; }
      } while (tot > 0.9);
      for (const auto& id : idx)
      {
        const auto ana = OrthogonalBasis::gradient(id, x);
        for (size_t k = 0; k < d; ++k)
        {
          const Real h = 1e-6;
          auto xp = x, xm = x;
          xp[(Eigen::Index)k] += h; xm[(Eigen::Index)k] -= h;
          const Real fd = (OrthogonalBasis::evaluate(id, xp)
                         - OrthogonalBasis::evaluate(id, xm)) / (2 * h);
          const Real scale = std::max({std::abs(fd), std::abs(ana[(Eigen::Index)k]), Real(1)});
          worst = std::max(worst, std::abs(fd - ana[(Eigen::Index)k]) / scale);
        }
      }
    }
    EXPECT_LT(worst, 1e-6) << (d == 2 ? "triangle" : "tetrahedron");
    std::printf("%s worst FD mismatch %.2e\n", d == 2 ? "tri" : "tet", worst);
  }
}
