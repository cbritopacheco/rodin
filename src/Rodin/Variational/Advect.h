/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_ADVECT_H
#define RODIN_VARIATIONAL_ADVECT_H

#include <mfem.hpp>

#include "ForwardDecls.h"

namespace Rodin::Variational
{
   namespace Internal
   {
      class DGSolver : public mfem::Solver
      {
         public:
            DGSolver(
                  mfem::SparseMatrix& M,
                  mfem::SparseMatrix& K,
                  const mfem::FiniteElementSpace& fes);

            void SetTimeStep(double dt);

            void SetOperator(const mfem::Operator& op) override;

            void Mult(const mfem::Vector& x, mfem::Vector& y) const override;

            virtual ~DGSolver() = default;

         private:
            mfem::SparseMatrix& m_M;
            mfem::SparseMatrix& m_K;
            mfem::GMRESSolver m_linearSolver;
            mfem::BlockILU m_prec;
            double m_dt;
      };

      class AdvectionEvolution : public mfem::TimeDependentOperator
      {
         public:
            AdvectionEvolution(mfem::BilinearForm& M, mfem::BilinearForm& K);

            virtual void Mult(
                  const mfem::Vector& x, mfem::Vector& y) const override;

            virtual void ImplicitSolve(
                  double dt, const mfem::Vector& x, mfem::Vector& k) override;

            virtual ~AdvectionEvolution() = default;

         private:
            mfem::BilinearForm& m_M;
            mfem::BilinearForm& m_K;
            mfem::CGSolver m_solver;
            std::unique_ptr<mfem::Solver> m_prec;
            std::unique_ptr<DGSolver> m_dgSolver;
            mutable mfem::Vector m_z;
      };
   }

   /**
    * @brief Advection of a function by an incompressible velocity field.
    *
    * Solves the advection equation
    * @f[
    * \left\{
    *    \begin{aligned}
    *       \dfrac{\partial u}{\partial t} + v(x) \cdot \nabla u (t, x) &= 0
    *          && \text{ in } \Omega \times (0, + \infty) \\
    *       u(x, 0) &= u_0(x) && \text{ on } \Omega \times \{ t = 0 \}
    *    \end{aligned}
    * \right.
    * @f]
    * where @f$ u_0 : \mathbb{R}^d \rightarrow \mathbb{R} @f$ is known and
    * @f$ v : \mathbb{R}^d \rightarrow \mathbb{R}^d @f$ is assumed to be
    * incompressible, i.e.  @f$ \nabla \cdot v = 0 @f$.
    */
   class Advect
   {
      public:
         template <class FEC>
         Advect(GridFunction<H1>& u, GridFunction<FEC>& velocity)
            : Advect(u, VectorCoefficient(velocity))
         {}

         /**
          * @brief Constructs the advection object.
          * @param[in,out] u Function to advect
          * @param[in] velocity Velocity field
          */
         Advect(GridFunction<H1>& u, const VectorCoefficientBase& velocity);

         void step(double dt);

         double getTime() const;

      private:
         double m_t;
         GridFunction<H1>& m_u;
         std::unique_ptr<VectorCoefficientBase> m_velocity;
         std::unique_ptr<mfem::ODESolver> m_odeSolver;
         mfem::DG_FECollection m_fec;
         mfem::FiniteElementSpace m_fes;
         mfem::BilinearForm m_M;
         mfem::BilinearForm m_K;
         Internal::AdvectionEvolution m_adv;
   };
}

#endif
