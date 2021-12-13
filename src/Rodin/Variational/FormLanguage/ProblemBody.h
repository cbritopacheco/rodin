/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FORMLANGUAGE_PROBLEMBODY_H
#define RODIN_VARIATIONAL_FORMLANGUAGE_PROBLEMBODY_H

#include <variant>

#include "TypeTraits.h"
#include "ForwardDecls.h"

#include "Base.h"

namespace Rodin::Variational::FormLanguage
{
   class ProblemBodyBase : public Base
   {
      public:
         virtual ProblemBodyBase& setProblem(ProblemBase& problem) = 0;
         virtual void eval() = 0;
         virtual ProblemBodyBase* copy() const noexcept = 0;
   };

   /**
    * @brief Represents the body of a variational problem.
    *
    * @code{.unparsed}
    *    ProblemBody ::= BilinearFormExpr
    *                  | BilinearFormExpr [+/-] LinearFormExpr
    *                  | BilinearFormExpr   +   BCExpr
    *                  | BilinearFormExpr [+/-] LinearFormExpr + BCExpr
    * @endcode
    */
   template <class Derived>
   class ProblemBody : public ProblemBodyBase
   {
      public:
         ProblemBody() = default;
         virtual ~ProblemBody() = default;

         ProblemBody(const ProblemBody&) = delete;
         void operator=(const ProblemBody& other) = delete;

         virtual ProblemBody& setProblem(ProblemBase& problem) override
         {
            static_cast<Derived*>(this)->setProblem(problem);
            return *this;
         }

         virtual void eval()
         {
            static_cast<Derived*>(this)->eval();
         }

         virtual ProblemBody* copy() const noexcept override
         {
            return static_cast<const Derived*>(this)->copy();
         }
   };
}

#endif
