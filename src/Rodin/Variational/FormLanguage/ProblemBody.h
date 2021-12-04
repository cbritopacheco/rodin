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

#include "RodinBase.h"

namespace Rodin::Variational::FormLanguage
{
   // template <class DerivedType>
   // struct TypeTraits<ProblemBody<DerivedType>>
   // {
   //    static constexpr SyntacticConstruct Syntax = Rule;
   //    using Derived = DerivedType;
   // };

   class ProblemBodyBase : public RodinBase
   {
      public:
         virtual ProblemBodyBase& setProblem(ProblemBase& problem) = 0;
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

         virtual void eval() override
         {
            static_cast<Derived*>(this)->eval();
         }

         template <class ... Args>
         static ProblemBody* create(Args&&... args) noexcept
         {
            return Derived::create(std::forward<Args>(args)...);
         }

         virtual ProblemBody* copy() const noexcept override
         {
            return static_cast<const Derived*>(this)->copy();
         }
   };
}

#endif
