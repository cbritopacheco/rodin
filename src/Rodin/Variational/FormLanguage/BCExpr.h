/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FORMLANGUAGE_BCEXPR_H
#define RODIN_VARIATIONAL_FORMLANGUAGE_BCEXPR_H

#include <memory>
#include <variant>

#include "TypeTraits.h"
#include "ForwardDecls.h"

#include "BCExprList.h"

namespace Rodin::Variational::FormLanguage
{
   // template <class DerivedType>
   // struct TypeTraits<BCExpr<DerivedType>>
   // {
   //    static constexpr SyntacticConstruct Syntax = Rule;
   //    using Derived = DerivedType;
   // };

   template <class Derived>
   class BCExpr : public BCExprList<BCExpr<Derived>>
   {
      public:
         BCExpr() = default;
         virtual ~BCExpr() = default;

         BCExpr(const BCExpr&) = delete;
         void operator=(const BCExpr& other) = delete;

         virtual BCExpr& setProblem(ProblemBase& problem) override
         {
            return static_cast<Derived*>(this)->setProblem(problem);
         }

         virtual void eval() override
         {
            static_cast<Derived*>(this)->eval();
         }

         template <class ... Args>
         static BCExpr* create(Args&&... args) noexcept
         {
            return Derived::create(std::forward<Args>(args)...);
         }

         virtual BCExpr* copy() const noexcept override
         {
            return static_cast<const Derived*>(this)->copy();
         }
   };
}

#endif
