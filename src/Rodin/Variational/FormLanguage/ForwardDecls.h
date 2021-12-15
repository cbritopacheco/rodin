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

   // ---- FormExprBCExprSum -------------------------------------------------
   template <class Derived>
   class ProblemBody;

   template <class Lhs, class Rhs>
   class BFExprLFExprSum;

   template <class Lhs, class Rhs>
   class BFExprBCExprSum;

   template <class Lhs, class Rhs>
   class FormExprBCExprSum;

   // ---- BilinearFormExpr --------------------------------------------------
   template <class Derived>
   class BilinearFormExpr;

   template <class Lhs, class Rhs>
   class BilinearFormExprSum;

   template <class T>
   class BilinearFormExprUnaryMinus;

   // ---- LinearFormExpr ----------------------------------------------------
   template <class Derived>
   class LinearFormExpr;

   template <class Lhs, class Rhs>
   class LinearFormExprSum;

   template <class T>
   class LinearFormExprUnaryMinus;

   // ---- BCExpr ------------------------------------------------------------
   template <class Derived>
   class BCExpr;

   template <class Derived>
   class BCExprList;

   template <class HeadNested, class Tail = void>
   class BCExprListCons;

   template <class Lhs, class Rhs>
   class FormBCExprSum;

   class ScalarSum;
   class ScalarCoefficientUnaryMinus;
}

#endif
