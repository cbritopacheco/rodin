/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Mult.h"

namespace Rodin::Variational
{
   Mult<ScalarCoefficientBase, ScalarCoefficientBase>
   operator*(const ScalarCoefficientBase& lhs, const ScalarCoefficientBase& rhs)
   {
      return Mult(lhs, rhs);
   }

   Mult<ScalarCoefficientBase, MatrixCoefficientBase>
   operator*(const ScalarCoefficientBase& lhs, const MatrixCoefficientBase& rhs)
   {
      return Mult(lhs, rhs);
   }

   Mult<ScalarCoefficientBase, MatrixCoefficientBase>
   operator*(const MatrixCoefficientBase& lhs, const ScalarCoefficientBase& rhs)
   {
      return Mult(rhs, lhs);
   }

   Mult<ScalarCoefficientBase, VectorCoefficientBase>
   operator*(const ScalarCoefficientBase& lhs, const VectorCoefficientBase& rhs)
   {
      return Mult(lhs, rhs);
   }

   Mult<ScalarCoefficientBase, VectorCoefficientBase>
   operator*(const VectorCoefficientBase& lhs, const ScalarCoefficientBase& rhs)
   {
      return Mult(rhs, lhs);
   }
}
