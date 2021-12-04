/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FORMLANGUAGE_BILINEARFORMEXPR_H
#define RODIN_VARIATIONAL_FORMLANGUAGE_BILINEARFORMEXPR_H

#include <variant>

#include "TypeTraits.h"
#include "ForwardDecls.h"

#include "ProblemBody.h"

namespace Rodin::Variational::FormLanguage
{
   // template <class DerivedType>
   // struct TypeTraits<BilinearFormExpr<DerivedType>>
   // {
   //    static constexpr SyntacticConstruct Syntax = Rule;
   //    using Derived = DerivedType;
   // };

   /**
    * @brief Represents an expression of bilinear forms.
    *
    * @code{.unparsed}
    * BilinearFormExpr ::= BilinearFormIntegrator
    *                    | BilinearFormExpr + BilinearFormExpr
    *                    | BilinearFormExpr - BilinearFormExpr
    *                    | - BilinearFormExpr
    * @endcode
    */
   template <class Derived>
   class BilinearFormExpr : public ProblemBody<BilinearFormExpr<Derived>>
   {
      public:
         BilinearFormExpr() = default;
         virtual ~BilinearFormExpr() = default;

         BilinearFormExpr(const BilinearFormExpr&) = delete;
         void operator=(const BilinearFormExpr& other) = delete;

         virtual BilinearFormExpr& setBilinearForm(BilinearFormBase& bf)
         {
            static_cast<Derived*>(this)->setBilinearForm(bf);
            return *this;
         }

         virtual void eval() override
         {
            static_cast<Derived*>(this)->eval();
         }

         virtual BilinearFormExpr& toggleSign()
         {
            static_cast<Derived*>(this)->toggleSign();
            return *this;
         }

         template <class ... Args>
         static BilinearFormExpr* create(Args&&... args) noexcept
         {
            return Derived::create(std::forward<Args>(args)...);
         }

         virtual BilinearFormExpr* copy() const noexcept override
         {
            return static_cast<const Derived*>(this)->copy();
         }
   };
}

#endif
