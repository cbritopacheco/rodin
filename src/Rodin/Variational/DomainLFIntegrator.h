/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_DOMAINLFINTEGRATOR_H
#define RODIN_VARIATIONAL_DOMAINLFINTEGRATOR_H

#include <functional>

#include "ForwardDecls.h"
#include "ScalarCoefficient.h"

#include "FormLanguage/LinearFormExpr.h"

namespace Rodin::Variational
{
   // template <class T>
   // struct FormLanguage::TypeTraits<DomainLFIntegrator<T>>
   // {
   //    static constexpr SyntacticConstruct Syntax = Constructor;
   //    using Rule = FormLanguage::LinearFormExpr<DiffusionIntegrator<T>>;
   // };

   /**
    *
    * @brief Object used to integrate a coefficient against a fixed coefficient.
    *
    * Represents the integration of the linear form
    * @f[
    *    L(v) = \int_{\Omega} f v \ dx
    * @f]
    * where @f$ f @f$ is a scalar coefficient.
    *
    * ----
    *
    * | Detail                | Description                                  |
    * |-----------------------|----------------------------------------------|
    * |  Spaces supported     | L2, H1                                       |
    * |  Dimensions supported | 1D, 2D, 3D                                   |
    * |  Continuous operator  | @f$ f @f$                                    |
    * |  @f$ \lambda @f$      | Scalar                                       |
    *
    * ----
    */
   template <class T = double>
   class DomainLFIntegrator
      :  public FormLanguage::LinearFormExpr<DomainLFIntegrator<T>>
   {
      public:
         /**
          * @brief Constructs a DomainLFIntegrator with scalar coefficient
          * `f`.
          *
          * @param[in] f Coefficient to integrate.
          */
         DomainLFIntegrator(const T& f = 1.0);

         DomainLFIntegrator(const DomainLFIntegrator& other);

         virtual DomainLFIntegrator& setLinearForm(LinearFormBase& lf) override;

         void eval() override;

         DomainLFIntegrator& toggleSign() override;

         template <class ... Args>
         static DomainLFIntegrator* create(Args&&... args) noexcept;

         virtual DomainLFIntegrator* copy() const noexcept override;

      private:
         ScalarCoefficient<T> m_f;
         std::optional<std::reference_wrapper<LinearFormBase>> m_lf;
   };
}

#include "DomainLFIntegrator.hpp"

#endif

