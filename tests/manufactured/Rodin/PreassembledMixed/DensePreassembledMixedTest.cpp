/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/// @file @brief Affine manufactured fields for dense preassembled mixed blocks.
#include <array>
#include <type_traits>
#include <boost/bimap.hpp>
#include <gtest/gtest.h>
#include "Rodin/Assembly.h"
#include "Rodin/Configure.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  class DensePreassembledMixed : public testing::TestWithParam<Polytope::Type>
  {
    protected:
      template <template <class, class> class Assembler>
      void check()
      {
        auto mesh = LocalMesh::UniformGrid(GetParam(), {3, 3});
        mesh.scale(0.5);
        mesh.getConnectivity().compute(2, 1);
        mesh.getConnectivity().compute(1, 0);
        mesh.getConnectivity().compute(1, 2);
        P1 first(mesh);
        P0 second(mesh);
        TrialFunction u(first);
        TestFunction v(first);
        TrialFunction p(second);
        TestFunction q(second);
        using Matrix = Math::Matrix<Real>;
        using System = Math::LinearSystem<Matrix, Math::Vector<Real>>;
        using Mixed = Problem<System, decltype(u), decltype(v), decltype(p), decltype(q)>;
        using Body = typename Mixed::ProblemBodyType;
        using Input = ::Rodin::Assembly::ProblemAssemblyInput<Body, decltype(u),
          decltype(v), decltype(p), decltype(q)>;
        using FirstSolution = typename decltype(u)::SolutionType;
        using SecondSolution = typename decltype(p)::SolutionType;
        BilinearForm<FirstSolution, decltype(first), decltype(first), Matrix> uu(u, v);
        BilinearForm<SecondSolution, decltype(second), decltype(first), Matrix> pu(p, v);
        BilinearForm<FirstSolution, decltype(first), decltype(second), Matrix> up(u, q);
        BilinearForm<SecondSolution, decltype(second), decltype(second), Matrix> pp(p, q);
        uu = 2 * Integral(u, v);
        pu = 3 * Integral(p, v);
        up = 5 * Integral(u, q);
        pp = 7 * Integral(p, q);
        uu.assemble();
        pu.assemble();
        up.assemble();
        pp.assemble();
        const RealFunction exactU(
          [](const Geometry::Point& point) { return 1 + 2 * point.x() - 3 * point.y(); });
        const RealFunction exactP(2.0);
        auto loadU = Integral(2 * exactU + 3 * exactP, v);
        auto loadP = Integral(5 * exactU + 7 * exactP, q);
        // Callable coefficient degree is unknown to automatic quadrature.
        // Affine forcing times a P1/P0 test function has degree two/one.
        loadU.setOrder(2);
        loadP.setOrder(1);
        Body preassembled = Body(uu) + pu + up + pp - loadU - loadP;
        Body reference = 2 * Integral(u, v) + 3 * Integral(p, v) + 5 * Integral(u, q) +
          7 * Integral(p, q) - loadU - loadP;
        auto trials = Tuple{std::ref(u), std::ref(p)};
        auto tests = Tuple{std::ref(v), std::ref(q)};
        const size_t total = first.getSize() + second.getSize();
        std::array<size_t, 2> trialOffsets{0, first.getSize()};
        std::array<size_t, 2> testOffsets{0, first.getSize()};
        boost::bimap<FormLanguage::Base::UUID, size_t> trialMap, testMap;
        trialMap.right.insert({0, u.getUUID()});
        trialMap.right.insert({1, p.getUUID()});
        testMap.right.insert({0, v.getUUID()});
        testMap.right.insert({1, q.getUUID()});
        Input assembledInput(preassembled, trials, tests, trialOffsets, testOffsets,
          trialMap, testMap, total, total);
        Input referenceInput(reference, trials, tests, trialOffsets, testOffsets,
          trialMap, testMap, total, total);
        System actual, expected;
        Assembler<System, Mixed> assembler;
        assembler.execute(actual, assembledInput);
        assembler.execute(expected, referenceInput);
        const auto& a = actual.getOperator();
        const auto& b = actual.getVector();
        ASSERT_EQ(a.rows(), first.getSize() + second.getSize());
        ASSERT_EQ(a.cols(), a.rows());
        ASSERT_EQ(b.size(), a.rows());
        // Small reference cells and affine data permit tight quadrature-level
        // agreement; solution checks allow conditioning of the coupled system.
        constexpr Real assemblyTolerance = 1e-12;
        constexpr Real solutionTolerance = 1e-11;
        EXPECT_LT((a - expected.getOperator()).norm() /
            std::max<Real>(1, expected.getOperator().norm()),
          assemblyTolerance);
        EXPECT_LT((b - expected.getVector()).norm() /
            std::max<Real>(1, expected.getVector().norm()),
          assemblyTolerance);
        GridFunction firstExact(first), secondExact(second);
        firstExact = exactU;
        secondExact = exactP;
        Math::Vector<Real> exact(a.cols());
        exact.head(first.getSize()) = firstExact.getData();
        exact.tail(second.getSize()) = secondExact.getData();
        EXPECT_LT(
          (a * exact - b).norm() / std::max<Real>(1, b.norm()), assemblyTolerance);
        const auto factorization = a.fullPivLu();
        ASSERT_EQ(factorization.rank(), a.cols());
        const Math::Vector<Real> solution = factorization.solve(b);
        EXPECT_LT(
          (solution - exact).norm() / std::max<Real>(1, exact.norm()), solutionTolerance);

        // Targeted assembly must preserve the opposite side, including when
        // called repeatedly on a system containing unrelated previous values.
        const Matrix previousMatrix = Matrix::Constant(total, total, -17);
        const Math::Vector<Real> previousRHS = Math::Vector<Real>::Constant(total, 23);
        for (const auto target : {AssemblyTarget::LHS, AssemblyTarget::RHS})
        {
          SCOPED_TRACE(static_cast<int>(target));
          System targeted;
          targeted.getOperator() = previousMatrix;
          targeted.getVector() = previousRHS;
          for (size_t repeat = 0; repeat < 2; ++repeat)
          {
            SCOPED_TRACE(repeat);
            assembler.execute(targeted, assembledInput, target);
            if (target == AssemblyTarget::LHS)
            {
              EXPECT_LT((targeted.getOperator() - a).norm() / std::max<Real>(1, a.norm()),
                assemblyTolerance);
              EXPECT_EQ((targeted.getVector() - previousRHS).norm(), 0);
            }
            else
            {
              EXPECT_EQ((targeted.getOperator() - previousMatrix).norm(), 0);
              EXPECT_LT((targeted.getVector() - b).norm() / std::max<Real>(1, b.norm()),
                assemblyTolerance);
            }
          }
        }

        if constexpr (std::is_same_v<Assembler<System, Mixed>,
                        typename Mixed::AssemblyType>)
        {
          // Trial blocks remain P1/P0, but test blocks are P0/P1. Their sizes
          // differ, so public assembly must compute the two offset arrays independently.
          using Permuted =
            Problem<System, decltype(u), decltype(q), decltype(p), decltype(v)>;
          Matrix permutedMatrix(total, total);
          permutedMatrix.topRows(second.getSize()) = a.bottomRows(second.getSize());
          permutedMatrix.bottomRows(first.getSize()) = a.topRows(first.getSize());
          Math::Vector<Real> permutedRHS(total);
          permutedRHS.head(second.getSize()) = b.tail(second.getSize());
          permutedRHS.tail(first.getSize()) = b.head(first.getSize());
          ASSERT_NE(first.getSize(), second.getSize());
          for (const auto target : {AssemblyTarget::LHS, AssemblyTarget::RHS})
          {
            SCOPED_TRACE(static_cast<int>(target));
            Permuted problem(u, q, p, v);
            problem = preassembled;
            auto& system = problem.getLinearSystem();
            system.getOperator() = previousMatrix;
            system.getVector() = previousRHS;
            ProblemBase<System>& interface = problem;
            for (size_t repeat = 0; repeat < 2; ++repeat)
            {
              SCOPED_TRACE(repeat);
              EXPECT_EQ(&interface.assemble(target), &interface);
              EXPECT_EQ(problem.getTrialOffsets(), trialOffsets);
              EXPECT_EQ(
                problem.getTestOffsets(), (std::array<size_t, 2>{0, second.getSize()}));
              if (target == AssemblyTarget::LHS)
              {
                EXPECT_LT((system.getOperator() - permutedMatrix).norm() /
                    std::max<Real>(1, permutedMatrix.norm()),
                  assemblyTolerance);
                EXPECT_EQ((system.getVector() - previousRHS).norm(), 0);
              }
              else
              {
                EXPECT_EQ((system.getOperator() - previousMatrix).norm(), 0);
                EXPECT_LT((system.getVector() - permutedRHS).norm() /
                    std::max<Real>(1, permutedRHS.norm()),
                  assemblyTolerance);
              }
            }
          }
        }
      }
  };

  TEST_P(DensePreassembledMixed, Sequential)
  {
    check<::Rodin::Assembly::Sequential>();
  }

#ifdef RODIN_USE_OPENMP
  TEST_P(DensePreassembledMixed, OpenMP)
  {
    check<::Rodin::Assembly::OpenMP>();
  }
#endif

  INSTANTIATE_TEST_SUITE_P(TwoGeometries, DensePreassembledMixed,
    testing::Values(Polytope::Type::Triangle, Polytope::Type::Quadrilateral));
}
