/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "ScalarFunction.h"

#include "Rodin/Mesh/SubMesh.h"

#include "Utility.h"
#include "Restriction.h"
#include "Exceptions.h"

namespace Rodin::Variational
{
   // ---- FunctionBase ------------------------------------------------------
   ScalarFunction<FunctionBase>::ScalarFunction(const FunctionBase& nested)
      : m_nested(nested.copy())
   {
      if (nested.getRangeType() != RangeType::Scalar)
         UnexpectedRangeTypeException(RangeType::Scalar, nested.getRangeType()).raise();
   }

   ScalarFunction<FunctionBase>::ScalarFunction(const ScalarFunction& other)
      :  ScalarFunctionBase(other),
         m_nested(other.m_nested->copy())
   {}

   ScalarFunction<FunctionBase>::ScalarFunction(ScalarFunction&& other)
      : ScalarFunctionBase(std::move(other)),
        m_nested(std::move(other.m_nested))
   {}

   double ScalarFunction<FunctionBase>::getValue(
         mfem::ElementTransformation& trans,
         const mfem::IntegrationPoint& ip) const
   {
      mfem::DenseMatrix v;
      m_nested->getValue(v, trans, ip);
      return v(0, 0);
   }
}
