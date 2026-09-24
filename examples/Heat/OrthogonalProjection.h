// OrthogonalProjection.h
//
// L2 projection Pi_h onto a finite element space, the building block of
// orthogonal subgrid scales (OSGS). The class is templated on the space, so the same code serves
// scalar quantities (div u, u.grad T) and vector ones ((grad u) u, grad p).
#ifndef EXAMPLES_HEAT_ORTHOGONALPROJECTION_H
#define EXAMPLES_HEAT_ORTHOGONALPROJECTION_H

#include <string>

#include <Rodin/PETSc.h>
#include <Rodin/Solver.h>
#include <Rodin/Variational.h>

namespace Rodin::Examples::Heat
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

      /// @param prefix PETSc options prefix of the mass solve (CG + Jacobi
      ///        unless overridden on the command line).
      OrthogonalProjection(const FES& fes, const std::string& prefix)
        : m_trial(fes), m_test(fes), m_problem(m_trial, m_test), m_ksp(m_problem),
          m_projection(fes)
      {
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
      TrialFunctionType m_trial;
      TestFunctionType m_test;
      ProblemType m_problem;
      Solver::KSP m_ksp;
      GridFunctionType m_projection;
  };
}

#endif
