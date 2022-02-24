/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FORMLANGUAGE_FORWARDDECLS_H
#define RODIN_VARIATIONAL_FORMLANGUAGE_FORWARDDECLS_H

#include "../ForwardDecls.h"

namespace Rodin::Variational::FormLanguage
{
   class Base;

   template <class InternalValue>
   class Buildable;

   class ProblemBody;

   class BilinearFormIntegratorSum;

   template <class T>
   class BilinearFormIntegratorUnaryMinus;

   class LinearFormIntegratorSum;

   template <class Lhs, class Rhs>
   class Sum;

   template <class Lhs, class Rhs>
   class Mult;

   template <class Operand>
   class UnaryMinus;
}

#endif
