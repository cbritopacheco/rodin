/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Physical and lifecycle contracts for PETSc single-field affine ties. */
#include <gtest/gtest.h>
#include <petsc.h>
#include <numeric>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/PETSc.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  /// @brief Affine component ties recover the exact field and preserve rejected targets.
  template <template <class, class> class Assembler, class Scalar = PetscScalar>
  void checkSingleFieldAffineIdentification()
  {
    const Scalar factor = [] {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return Scalar(1, -2);
      else
        return Scalar(1);
    }();
    const Scalar coefficient = [] {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return Scalar(2, 0.5);
      else
        return Scalar(2);
    }();
    for (const auto geometry : {Polytope::Type::Triangle, Polytope::Type::Quadrilateral})
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      auto mesh = Mesh<Context::Local>::UniformGrid(geometry, {2, 2});
      mesh.getConnectivity().compute(1, 0);
      mesh.getConnectivity().compute(1, 2);
      mesh.getConnectivity().compute(2, 1);
      auto check = [&](auto& space) {
        PETSc::Variational::TrialFunction u(space);
        PETSc::Variational::TestFunction v(space);
        using System = PETSc::Math::LinearSystem;
        using ProblemType = Problem<System, decltype(u), decltype(v)>;
        using Body = typename ProblemType::ProblemBodyType;
        Math::SpatialMatrix<Scalar> transform(2, 2);
        transform(0, 0) = 0;
        transform(0, 1) = coefficient;
        transform(1, 0) = 0;
        transform(1, 1) = 1;
        const auto load = VectorFunction(size_t(2), [factor, coefficient](const Point&) {
          return Math::SpatialVector<Scalar>{(coefficient + Scalar(3)) * factor, factor};
        });
        const auto offset = VectorFunction(size_t(2), [factor](const Point&) {
          return Math::SpatialVector<Scalar>{Scalar(3) * factor, Scalar(0)};
        });
        auto mass = Integral(u, v);
        auto forcing = Integral(load, v);
        mass.setOrder(8);
        forcing.setOrder(8);
        Body body =
          mass - forcing + DirichletBC(u, MatrixFunction(transform) * u, offset);
        Assembly::ProblemAssemblyInput<Body, decltype(u), decltype(v)> input(body, u, v);
        System system(PETSC_COMM_SELF);
        Assembler<System, ProblemType> assembler;
        GridFunction exact(space), direction(space);
        exact = load;
        direction = VectorFunction(size_t(2), [coefficient](const Point& point) {
          const Scalar value = Scalar(1 + point.x() + 2 * point.y());
          return Math::SpatialVector<Scalar>{coefficient * value, value};
        });
        const PetscInt size = static_cast<PetscInt>(space.getSize());
        std::vector<PetscInt> indices(static_cast<size_t>(size));
        std::iota(indices.begin(), indices.end(), PetscInt(0));
        auto matrix = [&] {
          Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> result(
            size, size);
          EXPECT_EQ(MatGetValues(system.getOperator(), size, indices.data(), size,
                      indices.data(), result.data()),
            PETSC_SUCCESS);
          return Math::Matrix<Scalar>(result);
        };
        auto vector = [&] {
          Math::Vector<Scalar> result(size);
          EXPECT_EQ(VecGetValues(system.getVector(), size, indices.data(), result.data()),
            PETSC_SUCCESS);
          return result;
        };
        for (size_t repeat = 0; repeat < 2; ++repeat)
        {
          assembler.execute(system, input);
          const auto assembled = matrix();
          const auto rhs = vector();
          EXPECT_LT((assembled * exact.getData() - rhs).norm(), 1e-11);
          const Math::Vector<Scalar> recovered = assembled.fullPivLu().solve(rhs);
          EXPECT_LT((recovered - exact.getData()).norm(), 1e-10);
          // Integral of |(coefficient*f,f)|^2, f = 1 + x + 2y, on the unit square.
          EXPECT_LT(std::abs(direction.getData().dot(assembled * direction.getData()) -
                      Scalar((std::norm(coefficient) + Real(1)) * Real(20) / Real(3))),
            1e-10);
          for (const auto target : {AssemblyTarget::LHS, AssemblyTarget::RHS})
          {
            EXPECT_THROW(assembler.execute(system, input, target), Alert::Exception);
            EXPECT_LT((matrix() - assembled).norm(), 1e-14);
            EXPECT_LT((vector() - rhs).norm(), 1e-14);
          }
        }
      };
      P1<Math::SpatialVector<Scalar>> linear(mesh, 2);
      check(linear);
      H1<2, Math::SpatialVector<Scalar>> quadratic(
        std::integral_constant<size_t, 2>{}, mesh, 2);
      check(quadratic);
    }
  }
}

TEST(Rodin_PETSc_SingleFieldConstraints, Sequential)
{
  checkSingleFieldAffineIdentification<PETSc::Assembly::Sequential>();
}

#ifdef RODIN_USE_OPENMP
TEST(Rodin_PETSc_SingleFieldConstraints, OpenMP)
{
  checkSingleFieldAffineIdentification<PETSc::Assembly::OpenMP>();
}
#endif

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  ::testing::InitGoogleTest(&argc, argv);
  const int result = RUN_ALL_TESTS();
  PetscFinalize();
  return result;
}
