/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Complex named-form insertion on PETSc local and distributed matrices.
 */
#include <gtest/gtest.h>
#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"
#include "Rodin/PETSc.h"
#include "Rodin/Solid/Linear/LinearElasticityForm.h"
#include "Rodin/Solid/Linear/LinearElasticityIntegral.h"
#ifdef RODIN_USE_MPI
#include "Rodin/Geometry/BalancedCompactPartitioner.h"
#include "Rodin/MPI/Geometry/Sharder.h"
#include "Rodin/MPI/Variational/P1.h"
#include "Rodin/MPI/Variational/H1.h"
#endif

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
#ifdef PETSC_USE_COMPLEX
  /// @brief Relative Frobenius roundoff budget for matched quadrature operators.
  constexpr PetscReal OperatorTolerance = 1e-11;

  /**
   * @brief Compares globally assembled operators without changing their entries.
   *
   * Distributed norms are intentionally collective: the regression certifies
   * the complete operator, including off-process contributions.
   *
   * @param actual Named-form matrix.
   * @param expected Explicit Integral matrix with matched quadrature.
   */
  void compare(::Mat actual, ::Mat expected)
  {
    ::Mat difference = nullptr;
    ASSERT_EQ(MatDuplicate(actual, MAT_COPY_VALUES, &difference), PETSC_SUCCESS);
    ASSERT_EQ(
      MatAXPY(difference, -1, expected, DIFFERENT_NONZERO_PATTERN), PETSC_SUCCESS);
    PetscReal error, norm;
    ASSERT_EQ(MatNorm(difference, NORM_FROBENIUS, &error), PETSC_SUCCESS);
    ASSERT_EQ(MatNorm(expected, NORM_FROBENIUS, &norm), PETSC_SUCCESS);
    EXPECT_LE(error, OperatorTolerance * std::max(PetscReal(1), norm));
    ASSERT_EQ(MatDestroy(&difference), PETSC_SUCCESS);
  }

  /**
   * @brief Exercises first insertion and repeated assembly with live coefficients.
   *
   * Mass, diffusion and Helmholtz use genuinely complex callable coefficients
   * on P1 and P2. Local meshes additionally exercise explicit sequential and
   * OpenMP assemblers. Distributed meshes use the MPI named assembler.
   * These are implemented-rule comparisons, not continuum-exactness tests.
   *
   * @tparam Distributed Whether the space has distributed ownership.
   * @tparam FES Complex scalar finite element space.
   * @param space Space whose operators are compared.
   */
  template <bool Distributed = false, class FES>
  void check(FES& space)
  {
    PETSc::Variational::TrialFunction trial(space);
    PETSc::Variational::TestFunction test(space);
    Real scale = 2;
    const ComplexFunction coefficient(
      [&](const Point& point) { return Complex(scale + point.x(), -scale / 2); });
    MassForm mass(coefficient, trial, test);
    DiffusionForm diffusion(coefficient, trial, test);
    HelmholtzForm helmholtz(coefficient, Real(3), trial, test);
    const auto cell = *space.getMesh().getCell();
    const size_t degree =
      space.getFiniteElement(cell.getDimension(), cell.getIndex()).getOrder();
    auto massIntegral = Integral(coefficient * trial, test);
    massIntegral.setOrder(2 * degree);
    auto diffusionIntegral = Integral(coefficient * Grad(trial), Grad(test));
    diffusionIntegral.setOrder(2 * degree - 2);
    auto reactionIntegral = Integral(Real(3) * trial, test);
    reactionIntegral.setOrder(2 * degree);
    const auto verify = [&](auto& form, const auto& integral) {
      BilinearForm expected(trial, test);
      expected = integral;
      expected.assemble();
      compare(form.getOperator(), expected.getOperator());
      form.assemble();
      compare(form.getOperator(), expected.getOperator());
      if constexpr (!Distributed)
      {
        ::Mat sequential = nullptr;
        Assembly::Sequential<::Mat, std::decay_t<decltype(form)>> assembler;
        for (size_t repeat = 0; repeat < 2; ++repeat)
        {
          assembler.execute(sequential, form);
          compare(sequential, expected.getOperator());
        }
        ASSERT_EQ(MatDestroy(&sequential), PETSC_SUCCESS);
#ifdef RODIN_USE_OPENMP
        ::Mat parallel = nullptr;
        Assembly::OpenMP<::Mat, std::decay_t<decltype(form)>> threaded;
        for (size_t repeat = 0; repeat < 2; ++repeat)
        {
          threaded.setThreadCount(2).execute(parallel, form);
          compare(parallel, expected.getOperator());
        }
        ASSERT_EQ(MatDestroy(&parallel), PETSC_SUCCESS);
#endif
      }
    };
    for (const Real value : {2., 5., 0.})
    {
      scale = value;
      mass.assemble();
      diffusion.assemble();
      helmholtz.assemble();
      verify(mass, massIntegral);
      verify(diffusion, diffusionIntegral);
      verify(helmholtz, diffusionIntegral + reactionIntegral);
    }
  }

  /**
   * @brief Certifies complex vector elasticity across named assembly backends.
   *
   * Nonreal Lamé parameters change between assemblies. The independent
   * integrator supplies the same quadrature rule; comparisons include first
   * insertion and reassembly, explicit sequential/OpenMP local assemblers,
   * and distributed insertion when ownership is present. Matrix norms are
   * collective for distributed spaces; no collective is added to the kernel.
   *
   * @tparam Distributed Whether the vector space has distributed ownership.
   * @tparam FES Complex vector finite element space.
   * @param space Space whose elasticity operators are compared.
   */
  template <bool Distributed = false, class FES>
  void checkElasticity(FES& space)
  {
    PETSc::Variational::TrialFunction trial(space);
    PETSc::Variational::TestFunction test(space);
    Real scale = 2;
    const ComplexFunction lambda(
      [&](const Point& point) { return Complex(scale + point.x(), scale / 2); });
    const ComplexFunction mu([&](const Point& point) {
      return Complex(2 * scale + point.x() * point.x(), -scale / 3);
    });
    LinearElasticityForm form(lambda, mu, trial, test);
    BilinearForm expected(trial, test);
    for (const Real value : {2., 5., 0.})
    {
      scale = value;
      expected = LinearElasticityIntegral(trial, test)(lambda, mu);
      expected.assemble();
      form.assemble();
      compare(form.getOperator(), expected.getOperator());
      form.assemble();
      compare(form.getOperator(), expected.getOperator());
      if constexpr (!Distributed)
      {
        ::Mat sequential = nullptr;
        Assembly::Sequential<::Mat, std::decay_t<decltype(form)>> assembler;
        for (size_t repeat = 0; repeat < 2; ++repeat)
        {
          assembler.execute(sequential, form);
          compare(sequential, expected.getOperator());
        }
        ASSERT_EQ(MatDestroy(&sequential), PETSC_SUCCESS);
#ifdef RODIN_USE_OPENMP
        ::Mat parallel = nullptr;
        Assembly::OpenMP<::Mat, std::decay_t<decltype(form)>> threaded;
        for (size_t repeat = 0; repeat < 2; ++repeat)
        {
          threaded.setThreadCount(2).execute(parallel, form);
          compare(parallel, expected.getOperator());
        }
        ASSERT_EQ(MatDestroy(&parallel), PETSC_SUCCESS);
#endif
      }
    }
  }

  /**
   * @brief Creates a small connected mesh for named assembly.
   * @param geometry Cell family.
   * @param points Grid points per axis; at least two.
   * @returns Mesh with incidence needed by both P1 and P2.
   */
  LocalMesh makeMesh(Polytope::Type geometry, size_t points = 3)
  {
    const size_t dimension = Polytope::Traits(geometry).getDimension();
    auto mesh =
      LocalMesh::UniformGrid(geometry, Array<size_t>::Constant(dimension, points));
    for (size_t d = 0; d <= dimension; ++d)
    {
      for (size_t dp = 0; dp <= dimension; ++dp)
        mesh.getConnectivity().compute(d, dp);
    }
    return mesh;
  }
#endif

  TEST(PETSc_ComplexNamedForm, Local)
  {
#ifdef PETSC_USE_COMPLEX
    for (const auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
           Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
           Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      auto mesh = makeMesh(geometry);
      P1<Complex> linear(mesh);
      check(linear);
      H1<2, Complex> quadratic(std::integral_constant<size_t, 2>{}, mesh);
      check(quadratic);
    }
#else
    GTEST_SKIP() << "Requires complex PETSc.";
#endif
  }

  /// @brief Checks P1/P2 complex vector elasticity on every affine cell family.
  TEST(PETSc_ComplexNamedForm, LocalLinearElasticity)
  {
#ifdef PETSC_USE_COMPLEX
    for (const auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
           Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
           Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      auto mesh = makeMesh(geometry, 2);
      const size_t dimension = mesh.getSpaceDimension();
      P1<Math::SpatialVector<Complex>> linear(mesh, dimension);
      checkElasticity(linear);
      H1<2, Math::SpatialVector<Complex>> quadratic(
        std::integral_constant<size_t, 2>{}, mesh, dimension);
      checkElasticity(quadratic);
    }
#else
    GTEST_SKIP() << "Requires complex PETSc.";
#endif
  }

#ifdef RODIN_USE_MPI
  TEST(PETSc_ComplexNamedForm, Distributed)
  {
#ifdef PETSC_USE_COMPLEX
    boost::mpi::environment environment;
    boost::mpi::communicator world;
    Context::MPI context(environment, world);
    Sharder<Context::MPI> sharder(context);
    if (world.rank() == 0)
    {
      auto local = makeMesh(Polytope::Type::Tetrahedron);
      BalancedCompactPartitioner partitioner(local);
      partitioner.partition(world.size());
      sharder.shard(partitioner);
    }
    sharder.scatter(0);
    auto mesh = sharder.gather(0);
    P1<Complex, decltype(mesh)> linear(mesh);
    check<true>(linear);
    H1<2, Complex, decltype(mesh)> quadratic(std::integral_constant<size_t, 2>{}, mesh);
    check<true>(quadratic);
#else
    GTEST_SKIP() << "Requires complex PETSc.";
#endif
  }

  /**
   * @brief Checks complex distributed elasticity on all seven affine cell families.
   *
   * Segments use at least one cell per rank. The test exercises nonempty
   * ownership and interface insertion, not the separate empty-shard contract.
   */
  TEST(PETSc_ComplexNamedForm, DistributedLinearElasticity)
  {
#ifdef PETSC_USE_COMPLEX
    boost::mpi::environment environment;
    boost::mpi::communicator world;
    Context::MPI context(environment, world);
    for (const auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
           Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
           Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      Sharder<Context::MPI> sharder(context);
      if (world.rank() == 0)
      {
        const size_t points = geometry == Polytope::Type::Segment
          ? std::max(size_t(3), static_cast<size_t>(world.size()) + 1)
          : size_t(3);
        auto local = makeMesh(geometry, points);
        BalancedCompactPartitioner partitioner(local);
        partitioner.partition(world.size());
        sharder.shard(partitioner);
      }
      sharder.scatter(0);
      auto mesh = sharder.gather(0);
      const size_t dimension = mesh.getSpaceDimension();
      P1<Math::SpatialVector<Complex>, decltype(mesh)> linear(mesh, dimension);
      checkElasticity<true>(linear);
      H1<2, Math::SpatialVector<Complex>, decltype(mesh)> quadratic(
        std::integral_constant<size_t, 2>{}, mesh, dimension);
      checkElasticity<true>(quadratic);
    }
#else
    GTEST_SKIP() << "Requires complex PETSc.";
#endif
  }
#endif
}

int main(int argc, char** argv)
{
#ifdef RODIN_USE_MPI
  boost::mpi::environment environment(argc, argv);
#endif
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  ::testing::InitGoogleTest(&argc, argv);
  const int result = RUN_ALL_TESTS();
  PetscFinalize();
  return result;
}
