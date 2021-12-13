/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FORMLANGUAGE_BCEXPRLIST_H
#define RODIN_VARIATIONAL_FORMLANGUAGE_BCEXPRLIST_H

#include <memory>
#include <variant>

#include "TypeTraits.h"
#include "ForwardDecls.h"

#include "Base.h"

namespace Rodin::Variational::FormLanguage
{
   template <class Derived>
   class BCExprList : public Base
   {
      public:
         BCExprList() = default;
         virtual ~BCExprList() = default;

         BCExprList(const BCExprList&) = delete;
         void operator=(const BCExprList& other) = delete;

         virtual BCExprList& setProblem(ProblemBase& problem)
         {
            return static_cast<Derived*>(this)->setProblem(problem);
         }

         virtual void eval()
         {
            static_cast<Derived*>(this)->eval();
         }

         virtual BCExprList* copy() const noexcept override
         {
            return static_cast<const Derived*>(this)->copy();
         }
   };
}

#endif

