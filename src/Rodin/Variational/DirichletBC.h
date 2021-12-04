/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_DIRICHLETBC_H
#define RODIN_VARIATIONAL_DIRICHLETBC_H

#include <functional>

#include "FormLanguage/TypeTraits.h"
#include "FormLanguage/BCExpr.h"

#include "Problem.h"
#include "ScalarCoefficient.h"

namespace Rodin::Variational
{
   // template <class T>
   // struct FormLanguage::TypeTraits<DirichletBC<T>>
   // {
   //    static constexpr SyntacticConstruct Syntax = Constructor;
   //    using Rule = FormLanguage::BCExpr<DirichletBC<T>>;
   // };

   /**
    * @brief Represents a Dirichlet boundary condition.
    *
    * Imposing a Dirichlet boundary condition amounts to projecting the
    * boundary value over the function on the true boundary degrees of freedom.
    * In simpler terms, the following condition will be satisfied on the
    * segment @f$ \Gamma_D @f$ of the whole boundary @f$ \Omega @f$ specified.
    *
    * @note This class is meant to be instanced as a part of a @ref
    * Variational::FormLanguage expression.
    *
    * Example
    * -------
    *  Take the simple Poisson equation with Dirichlet boundary conditions:
    * @f[
    * \left\{
    *  \begin{aligned}
    *    - \Delta u &= 1 && \mathrm{in} \ \Omega \\
    *    u &= 0 && \mathrm{on} \ \Gamma_D \\
    *  \end{aligned}
    * \right.
    * @f]
    *
    * Integrating by parts we arrive at the variational formulation.
    * @f[
    * \text{Find } u \in H^1_0(\Omega) \text{ s.t. }
    * \forall v \in H^1_0(\Omega), \quad
    *    \int_{\Omega} \lambda \nabla u \cdot \nabla v \ dx - \int_{\Omega} v
    *    \ dx = 0
    * @f]
    *
    * In code this will look like:
    *
    * @code{.cpp}
    * Problem poisson(u, v);
    * poisson = DiffusionIntegrator(lambda)
    *         - DomainLFIntegrator(1.0)
    *         + DirichletBC(GammaD, 0.0);
    * @endcode
    *
    */
   template <class T>
   class DirichletBC
      :  public FormLanguage::BCExpr<DirichletBC<T>>
   {
      public:
         /**
          * @brief Constructs a Dirichlet boundary condition on the part of the
          * boundary specified by the boundary attribute.
          *
          * @param[in] bdrAttr Attribute corresponding to the segment @f$
          * \Gamma_D @f$ where the Dirichlet boundary condition is imposed.
          * @param[in] value Value of the trial function @f$ u @f$ at the
          * boundary @f$ \Gamma_D @f$ specified by `bdrAttr`.
          */
         DirichletBC(int bdrAttr, const T& value);

         DirichletBC(const DirichletBC& other);
         virtual ~DirichletBC() = default;

         DirichletBC& setProblem(ProblemBase& problem) override;

         void eval() override;

         template <class ... Args>
         static DirichletBC* create(Args&&... args) noexcept;
         virtual DirichletBC* copy() const noexcept override;

      private:
         int m_bdrAttr;
         ScalarCoefficient<T> m_value;
         mfem::Array<int> m_essBdr;
         std::optional<std::reference_wrapper<ProblemBase>> m_problem;
   };

   template <class ... Values>
   class DirichletBC<VectorCoefficient<Values...>>
      :  public FormLanguage::BCExpr<DirichletBC<VectorCoefficient<Values...>>>
   {
      public:
         DirichletBC(int bdrAttr, const VectorCoefficient<Values...>& value);
         DirichletBC(const DirichletBC& other);

         DirichletBC& setProblem(ProblemBase& problem) override;

         void eval() override;

         template <class ... Args>
         static DirichletBC* create(Args&&... args) noexcept;
         virtual DirichletBC* copy() const noexcept override;

      private:
         int m_bdrAttr;
         VectorCoefficient<Values...> m_value;
         mfem::Array<int> m_essBdr;
         std::optional<std::reference_wrapper<ProblemBase>> m_problem;
   };
}

#include "DirichletBC.hpp"

#endif
