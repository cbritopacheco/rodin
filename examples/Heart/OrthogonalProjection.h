/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file OrthogonalProjection.h
 * @brief L2 projection Pi_h onto a finite element space, the building block of
 *        orthogonal subgrid scales: P_h^perp[f] = f - Pi_h[f].
 *
 * Templated on the space, so one class serves scalars (div u), vectors
 * (grad p, div sigma, (grad u)u) and symmetric tensors stored in Voigt order
 * (a vector space of dimension d(d+1)/2). A Voigt field is projected component
 * by component, so the factor 2 that the Frobenius product puts on the
 * off-diagonal entries does not enter the projection.
 */
#ifndef RODIN_EXAMPLES_HEART_ORTHOGONALPROJECTION_H
#define RODIN_EXAMPLES_HEART_ORTHOGONALPROJECTION_H

#include <cassert>
#include <string>

#include <petscsys.h>

#include <Rodin/PETSc.h>
#include <Rodin/Solver.h>
#include <Rodin/Variational.h>

namespace Rodin::Examples::Heart
{
  template <class FES>
  class OrthogonalProjection
  {
    public:
      using GridFunctionType = PETSc::Variational::GridFunction<FES>;
      using TrialFunctionType = PETSc::Variational::TrialFunction<GridFunctionType, FES>;
      using TestFunctionType = PETSc::Variational::TestFunction<FES>;
      using ProblemType = Variational::Problem<PETSc::Math::LinearSystem,
        TrialFunctionType, TestFunctionType>;

      /// @param prefix PETSc options prefix of the mass solve; CG + Jacobi
      ///        unless overridden on the command line.
      OrthogonalProjection(const FES& fes, const std::string& prefix)
        : m_trial(fes), m_test(fes), m_problem(m_trial, m_test), m_ksp(m_problem),
          m_projection(fes)
      {
        setDefault(prefix, "ksp_type", "cg");
        setDefault(prefix, "pc_type", "jacobi");
        m_ksp.setPrefix(prefix);
      }

      OrthogonalProjection(const OrthogonalProjection&) = delete;
      OrthogonalProjection& operator=(const OrthogonalProjection&) = delete;

      /// @brief Pi_h[f]: find pi in V_h with (pi, w) = (f, w) for all w in V_h.
      template <class Expression>
      OrthogonalProjection& project(const Expression& f)
      {
        m_problem = Variational::Integral(m_trial, m_test) - Variational::Integral(f, m_test);
        m_problem.solve(m_ksp);
        m_projection.setData(m_trial.getSolution().getData());
        return *this;
      }

      /// @brief The last projection, usable inside any form.
      const GridFunctionType& get() const { return m_projection; }

      GridFunctionType& get() { return m_projection; }

    private:
      static void setDefault(const std::string& prefix, const char* key, const char* value)
      {
        const std::string name = "-" + prefix + key;
        PetscBool set = PETSC_FALSE;
        PetscErrorCode ierr =
          PetscOptionsHasName(PETSC_NULLPTR, PETSC_NULLPTR, name.c_str(), &set);
        assert(ierr == PETSC_SUCCESS);
        if (!set)
        {
          ierr = PetscOptionsSetValue(PETSC_NULLPTR, name.c_str(), value);
          assert(ierr == PETSC_SUCCESS);
        }
        (void)ierr;
      }

      TrialFunctionType m_trial;
      TestFunctionType m_test;
      ProblemType m_problem;
      Solver::KSP m_ksp;
      GridFunctionType m_projection;
  };
}

#endif
