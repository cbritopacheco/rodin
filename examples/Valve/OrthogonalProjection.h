// OrthogonalProjection.h
//
// L2 projection Pi_h onto a finite element space, the building block of
// orthogonal subgrid scales (OSGS). The class is templated on the space, so the same code serves
// scalar quantities (div u, u.grad T) and vector ones ((grad u) u, grad p).
//
// Two modes: the consistent projection (mass matrix solved with the KSP of
// the given prefix, CG + Jacobi by default) and the lumped one, Pi = b ./ m_L
// with m_L the row sums of the mass matrix, computed once: one load-vector
// assembly per projection and no solve.
#ifndef EXAMPLES_VALVE_ORTHOGONALPROJECTION_H
#define EXAMPLES_VALVE_ORTHOGONALPROJECTION_H

#include <cassert>
#include <string>
#include <type_traits>

#include <Rodin/PETSc.h>
#include <Rodin/Solver.h>
#include <Rodin/Variational.h>

namespace Rodin::Examples::Valve
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
      using LinearFormType = Variational::LinearForm<FES, ::Vec>;

      /// @param prefix PETSc options prefix of the mass solve (CG + Jacobi
      ///        unless overridden on the command line).
      /// @param lumped Use the lumped mass matrix (no solve).
      OrthogonalProjection(const FES& fes, const std::string& prefix, bool lumped = false)
        : m_trial(fes), m_test(fes), m_problem(m_trial, m_test), m_ksp(m_problem),
          m_projection(fes), m_lumped(lumped), m_rhs(m_test)
      {
        m_ksp.setPrefix(prefix);
        if (m_lumped)
        {
          GridFunctionType one(fes);
          if constexpr (std::is_same_v<typename FES::RangeType, Real>)
            one = Real(1);
          else
          {
            Math::SpatialVector<Real> e(fes.getVectorDimension());
            for (size_t i = 0; i < fes.getVectorDimension(); ++i)
              e(i) = 1.0;
            one = e;
          }
          m_rhs = Variational::Integral(one, m_test);
          m_rhs.assemble();
          PetscErrorCode ierr = VecDuplicate(m_rhs.getVector(), &m_mass);
          assert(ierr == PETSC_SUCCESS);
          ierr = VecCopy(m_rhs.getVector(), m_mass);
          assert(ierr == PETSC_SUCCESS);
          (void)ierr;
        }
      }

      ~OrthogonalProjection()
      {
        if (m_mass)
        {
          PetscErrorCode ierr = VecDestroy(&m_mass);
          assert(ierr == PETSC_SUCCESS);
          (void)ierr;
        }
      }

      OrthogonalProjection(const OrthogonalProjection&) = delete;
      OrthogonalProjection& operator=(const OrthogonalProjection&) = delete;

      /// @brief Pi_h[f]: find pi in V_h with (pi, w) = (f, w) for all w in V_h
      ///        (lumped: m_L pi = (f, w)).
      template <class Expression>
      OrthogonalProjection& project(const Expression& f)
      {
        if (m_lumped)
        {
          m_rhs = Variational::Integral(f, m_test);
          m_rhs.assemble();
          ::Vec pi;
          PetscErrorCode ierr = VecDuplicate(m_rhs.getVector(), &pi);
          assert(ierr == PETSC_SUCCESS);
          ierr = VecPointwiseDivide(pi, m_rhs.getVector(), m_mass);
          assert(ierr == PETSC_SUCCESS);
          m_projection.setData(pi);
          ierr = VecDestroy(&pi);
          assert(ierr == PETSC_SUCCESS);
          (void)ierr;
          return *this;
        }
        m_problem = Variational::Integral(m_trial, m_test) - Variational::Integral(f, m_test);
        m_problem.solve(m_ksp);
        m_projection.setData(m_trial.getSolution().getData());
        return *this;
      }

      /// @brief The last projection, usable inside any form.
      const GridFunctionType& get() const { return m_projection; }

      GridFunctionType& get() { return m_projection; }

    private:
      TrialFunctionType m_trial;
      TestFunctionType m_test;
      ProblemType m_problem;
      Solver::KSP m_ksp;
      GridFunctionType m_projection;
      bool m_lumped;
      LinearFormType m_rhs;
      ::Vec m_mass = PETSC_NULLPTR;
  };
}

#endif
