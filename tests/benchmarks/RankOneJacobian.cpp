/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#include <benchmark/benchmark.h>

#include "Rodin/Assembly.h"
#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"
#include "Rodin/Variational/H1.h"
#include "Rodin/Variational/P1.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Benchmarks
{
  namespace
  {
    class OtherIdentityMatrix final : public MatrixFunctionBase<Real, OtherIdentityMatrix>
    {
      public:
        explicit OtherIdentityMatrix(size_t dimension)
          : m_dimension(dimension)
        {}

        Math::SpatialMatrix<Real> getValue(const Geometry::Point&) const
        {
          return Math::SpatialMatrix<Real>::Identity(
            static_cast<std::uint8_t>(m_dimension),
            static_cast<std::uint8_t>(m_dimension));
        }

        size_t getRows() const
        {
          return m_dimension;
        }
        size_t getColumns() const
        {
          return m_dimension;
        }

        Optional<size_t> getOrder(const Geometry::Polytope&) const noexcept
        {
          return 0;
        }

        OtherIdentityMatrix* copy() const noexcept override
        {
          return new OtherIdentityMatrix(*this);
        }

      private:
        size_t m_dimension;
    };

    LocalMesh unitSquare(size_t n)
    {
      LocalMesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {n, n});
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 0);
      return mesh;
    }

    template <class FES>
    void specializedAssembly(benchmark::State& state, FES& fes)
    {
      TrialFunction u(fes);
      TestFunction v(fes);
      const IdentityMatrix A(2);
      BilinearForm form(u, v);
      form = Integral(Dot(Dot(A, Jacobian(u)), Dot(A, Jacobian(v))));

      for (auto _ : state)
      {
        form.assemble();
        benchmark::DoNotOptimize(form.getOperator().nonZeros());
      }
    }

    template <class FES>
    void genericAssembly(benchmark::State& state, FES& fes)
    {
      TrialFunction u(fes);
      TestFunction v(fes);
      const IdentityMatrix A(2);
      const OtherIdentityMatrix B(2);
      BilinearForm form(u, v);
        // The distinct coefficient types deliberately select the generic rule.
      form = Integral(Dot(Dot(A, Jacobian(u)), Dot(B, Jacobian(v))));

      for (auto _ : state)
      {
        form.assemble();
        benchmark::DoNotOptimize(form.getOperator().nonZeros());
      }
    }
  }

  static void P1RankOneSpecialized(benchmark::State& state)
  {
    auto mesh = unitSquare(32);
    P1 fes(mesh, 2);
    specializedAssembly(state, fes);
  }
  BENCHMARK(P1RankOneSpecialized);

  static void P1RankOneGeneric(benchmark::State& state)
  {
    auto mesh = unitSquare(32);
    P1 fes(mesh, 2);
    genericAssembly(state, fes);
  }
  BENCHMARK(P1RankOneGeneric);

  static void H1P2RankOneSpecialized(benchmark::State& state)
  {
    auto mesh = unitSquare(32);
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, 2);
    specializedAssembly(state, fes);
  }
  BENCHMARK(H1P2RankOneSpecialized);

  static void H1P2RankOneGeneric(benchmark::State& state)
  {
    auto mesh = unitSquare(32);
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, 2);
    genericAssembly(state, fes);
  }
  BENCHMARK(H1P2RankOneGeneric);
}
