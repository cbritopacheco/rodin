// OrthogonalProjection.h
//
// L2 projection Pi_h onto a finite element space, the building block of
// orthogonal subgrid scales (OSGS). The class is templated on the space, so the same code serves
// scalar quantities (div u, u.grad T) and vector ones ((grad u) u, grad p).
#ifndef EXAMPLES_HEAT_ORTHOGONALPROJECTION_H
#define EXAMPLES_HEAT_ORTHOGONALPROJECTION_H

#include <optional>
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

      /// @brief Restricts the projection to the cells of one region: (pi, w)
      ///        and (f, w) are integrated over it, plus a 1e-8 mass elsewhere
      ///        so that dofs outside it get pi = 0. For a field that only has
      ///        a meaning in that region (e.g. the pressure of a fluid next to
      ///        a solid), so that the rest does not pollute Pi_h near the
      ///        interface.
      OrthogonalProjection& restrictTo(Geometry::Attribute region)
      {
        m_region = region;
        return *this;
      }

      /// @brief Pi_h[f]: find pi in V_h with (pi, w) = (f, w) for all w in V_h.
      template <class Expression>
      OrthogonalProjection& project(const Expression& f)
      {
        if (m_region)
        {
          m_problem = Variational::Integral(m_trial, m_test).over(*m_region)
                    + 1.0e-8 * Variational::Integral(m_trial, m_test)
                    - Variational::Integral(f, m_test).over(*m_region);
        }
        else
        {
          m_problem = Variational::Integral(m_trial, m_test) - Variational::Integral(f, m_test);
        }
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
      std::optional<Geometry::Attribute> m_region;
  };
}

#endif
