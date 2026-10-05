/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_TESTS_UNIT_PETSC_SCALARINTEGRAL_H
#define RODIN_TESTS_UNIT_PETSC_SCALARINTEGRAL_H

#include <gtest/gtest.h>
#include <Rodin/PETSc.h>
#include <Rodin/Variational.h>

namespace Rodin::Tests::Unit
{
  /**
   * @brief Constant reproduction and scalar linearity of grade-zero integrals.
   *
   * On a unit-volume mesh, @f$\int u_h=c@f$ and
   * @f$\int z u_h=zc@f$. Complex builds deliberately use non-real values;
   * conjugate-linear form action must not leak into grade-zero integration.
   */
  class ScalarIntegral
  {
    public:
      template <class MeshType>
      static void check(const MeshType& mesh)
      {
        using namespace Variational;
        const auto space = [](const auto& fes) {
          constexpr Real tolerance = 1e-11; // Constant-reproduction roundoff budget.
#ifdef PETSC_USE_COMPLEX
          const PetscScalar value(2, 3), scale(-1, 2);
#else
          const PetscScalar value = 2, scale = -1;
#endif
          GridFunction<std::decay_t<decltype(fes)>, ::Vec> field(fes);
          field = value;
          auto integral = Integral(field);
          integral.setOrder(18);
          const auto verify = [&](PetscScalar expected) {
            const PetscScalar actual = integral.compute();
            EXPECT_NEAR(PetscRealPart(actual), PetscRealPart(expected), tolerance);
            EXPECT_NEAR(
              PetscImaginaryPart(actual), PetscImaginaryPart(expected), tolerance);
          };
          verify(value);
          field.sync();
          const PetscErrorCode ierr = VecScale(field.getData(), scale);
          assert(ierr == PETSC_SUCCESS);
          (void)ierr;
          field.sync();
          verify(scale * value);
        };
        space(P0<PetscScalar, MeshType>(mesh));
        space(P0g<PetscScalar, MeshType>(mesh));
        space(P1<PetscScalar, MeshType>(mesh));
        const auto order = [&]<size_t K>() {
          space(H1<K, PetscScalar, MeshType>(std::integral_constant<size_t, K>{}, mesh));
        };
        order.template operator()<1>();
        order.template operator()<2>();
        order.template operator()<3>();
        order.template operator()<4>();
        checkActions(mesh);
      }

    private:
      /** Non-Hermitian coefficients distinguish action from its conjugate. */
      template <class MeshType>
      static void checkActions(const MeshType& mesh)
      {
        using namespace Variational;
        constexpr Real tolerance = 1e-11; // Unit-volume assembly roundoff budget.
#ifdef PETSC_USE_COMPLEX
        const PetscScalar c(3, -4), a(2, 3), b(-1, 2), z(1, -2);
        const ComplexFunction coefficient(c);
#else
        const PetscScalar c = 3, a = 2, b = -1, z = 2;
        const RealFunction coefficient(c);
#endif
        P0g<PetscScalar, MeshType> fes(mesh);
        PETSc::Variational::GridFunction trial(fes), test(fes);
        trial = a;
        test = b;
        PETSc::Variational::TrialFunction u(fes);
        PETSc::Variational::TestFunction v(fes);
        LinearForm linear(v);
        linear = Integral(coefficient, v);
        linear.assemble();
        BilinearForm bilinear(u, v);
        bilinear = Integral(coefficient * u, v);
        bilinear.assemble();
        const auto verify = [&](PetscScalar actual, PetscScalar expected) {
          EXPECT_NEAR(PetscRealPart(actual), PetscRealPart(expected), tolerance);
          EXPECT_NEAR(
            PetscImaginaryPart(actual), PetscImaginaryPart(expected), tolerance);
        };
        verify(linear(test), c * Math::conj(b));
        verify(bilinear(trial, test), c * a * Math::conj(b));
        trial = z * a;
        verify(bilinear(trial, test), z * c * a * Math::conj(b));
        trial = a;
        test = z * b;
        verify(linear(test), Math::conj(z) * c * Math::conj(b));
        verify(bilinear(trial, test), Math::conj(z) * c * a * Math::conj(b));
      }
  };
}
#endif
