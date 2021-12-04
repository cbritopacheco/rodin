/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FORMLANGUAGE_LINEARFORMEXPR_H
#define RODIN_VARIATIONAL_FORMLANGUAGE_LINEARFORMEXPR_H

#include <variant>

#include "TypeTraits.h"
#include "ForwardDecls.h"

#include "ProblemBody.h"

#include "RodinBase.h"

namespace Rodin::Variational::FormLanguage
{
   // template <class DerivedType>
   // struct TypeTraits<LinearFormExpr<DerivedType>>
   // {
   //    static constexpr SyntacticConstruct Syntax = Rule;
   //    using Derived = DerivedType;
   // };

   /**
    * @brief Represents an expression of bilinear forms.
    *
    * @code{.unparsed}
    * BilinearFormExpr ::= LinearFormIntegrator
    *                    | LinearFormIntegrator + LinearFormIntegrator
    *                    | LinearFormIntegrator - LinearFormIntegrator
    *                    | - LinearFormIntegrator
    * @endcode
    */
   template <class Derived>
   class LinearFormExpr : public RodinBase
   {
      public:
         LinearFormExpr() = default;
         virtual ~LinearFormExpr() = default;

         LinearFormExpr(const LinearFormExpr&) = delete;
         void operator=(const LinearFormExpr& other) = delete;

         virtual LinearFormExpr& setLinearForm(LinearFormBase& bf)
         {
            static_cast<Derived*>(this)->setLinearForm(bf);
            return *this;
         }

         virtual void eval() override
         {
            static_cast<Derived*>(this)->eval();
         }

         virtual LinearFormExpr& toggleSign()
         {
            static_cast<Derived*>(this)->toggleSign();
            return *this;
         }

         template <class ... Args>
         static LinearFormExpr* create(Args&&... args) noexcept
         {
            return Derived::create(std::forward<Args>(args)...);
         }

         virtual LinearFormExpr* copy() const noexcept override
         {
            return static_cast<const Derived*>(this)->copy();
         }
   };
}

#endif
