/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Class @ref NeumannBC
 */
#ifndef RODIN_VARIATIONAL_NEUMANNBC_H
#define RODIN_VARIATIONAL_NEUMANNBC_H

#include <variant>
#include <functional>

#include "ForwardDecls.h"

#include "FormLanguage/TypeTraits.h"
#include "FormLanguage/BCExpr.h"

#include "Problem.h"
#include "ScalarCoefficient.h"

namespace Rodin::Variational
{
   // template <class T>
   // struct FormLanguage::TypeTraits<NeumannBC<T>>
   // {
   //    static constexpr SyntacticConstruct Syntax = Constructor;
   //    using Rule = FormLanguage::BCExpr<NeumannBC<T>>;
   // };

   /**
    * @brief Represents a Neumann boundary condition.
    *
    * The usage of this class assumes the proper coefficients are present
    * in the boundary conditions of the equation. In particular, when utilized
    * in a @ref Problem declaration, it will add a term of the form:
    * @f[
    *    \int_{\Gamma_N} g v \ dx
    * @f]
    * to the right hand side of the variational formulation.
    *
    * | Detail                | Description                                  |
    * |-----------------------|----------------------------------------------|
    * |  Spaces supported     | L2, H1                                       |
    * |  Dimensions supported | 1D, 2D, 3D                                   |
    *
    * @note This class is meant to be instanced as a part of a @ref
    * Variational::FormLanguage expression.
    *
    * Example
    * -------
    *  Take the simple heat conduction equation:
    * @f[
    * \left\{
    *  \begin{aligned}
    *    - \nabla \cdot (\lambda \nabla u) &= f && \mathrm{in} \ \Omega \\
    *    u &= 0 && \mathrm{on} \ \Gamma_D \\
    *    \lambda \partial_n u &= g && \mathrm{on} \ \Gamma_N
    *  \end{aligned}
    * \right.
    * @f]
    *
    * Integrating by parts we arrive at the variational formulation.
    * @f[
    * \text{Find } u \in H^1_{\Gamma_D}(\Omega) \text{ s.t. }
    * \forall v \in H^1_0(\Omega), \quad
    *    \int_{\Omega} \lambda \nabla u \cdot \nabla v \ dx - \int_{\Omega} fv
    *    \ dx = \int_{\Gamma_N} \underbrace{g}_{\lambda \partial_n u} v \ dx
    * @f]
    *
    * In code this will look like:
    *
    * @code{.cpp}
    * Problem poisson(u, v);
    * poisson = DiffusionIntegrator(lambda)
    *         - DomainLFIntegrator(f)
    *         + DirichletBC(GammaD, 0.0)
    *         + NeumannBC(GammaN, g);
    * @endcode
    *
    */
   template <class T>
   class NeumannBC
      :  public FormLanguage::BCExpr<NeumannBC<T>>
   {
      public:
         /**
          * @brief Constructs a Neumann boundary condition on the part of the
          * boundary specified by the boundary attribute.
          *
          * @param[in] bdrAttr Attribute corresponding to the segment @f$
          * \Gamma_N @f$ where the Neumann boundary condition is imposed.
          * @param[in] value Value of the normal derivative at the
          * boundary @f$ \Gamma_N @f$ specified by `bdrAttr`.
          */
         NeumannBC(int bdrAttr, const T& value);
         NeumannBC(const NeumannBC& other);

         NeumannBC& setProblem(ProblemBase& problem) override;

         void eval() override;

         template <class ... Args>
         static NeumannBC* create(Args&&... args) noexcept;
         virtual NeumannBC* copy() const noexcept override;

      private:
         int m_bdrAttr;
         ScalarCoefficient<T> m_value;
         mfem::Array<int> m_nbcBdr;
         std::optional<std::reference_wrapper<ProblemBase>> m_problem;
   };

   template <class ... Values>
   class NeumannBC<VectorCoefficient<Values...>>
      :  public FormLanguage::BCExpr<NeumannBC<VectorCoefficient<Values...>>>
   {
      public:
         NeumannBC(int bdrAttr, const VectorCoefficient<Values...>& value);
         NeumannBC(const NeumannBC& other);

         NeumannBC& setProblem(ProblemBase& problem) override;

         void eval() override;

         template <class ... Args>
         static NeumannBC* create(Args&&... args) noexcept;
         virtual NeumannBC* copy() const noexcept override;

      private:
         int m_bdrAttr;
         VectorCoefficient<Values...> m_value;
         mfem::Array<int> m_nbcBdr;
         std::optional<std::reference_wrapper<ProblemBase>> m_problem;
   };
}

#include "NeumannBC.hpp"

#endif
