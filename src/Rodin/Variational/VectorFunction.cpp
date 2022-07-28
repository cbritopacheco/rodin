/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "VectorFunction.h"

#include "Component.h"
#include "GridFunction.h"

namespace Rodin::Variational
{
   Component<FunctionBase> VectorFunctionBase::operator()(int i) const
   {
      assert(0 <= i);
      assert(i < getDimension());
      return Component(*this, i);
   }

   Component<FunctionBase> VectorFunctionBase::x() const
   {
      assert(getDimension() >= 1);
      return Component<FunctionBase>(*this, 0);
   }

   Component<FunctionBase> VectorFunctionBase::y() const
   {
      assert(getDimension() >= 2);
      return Component<FunctionBase>(*this, 1);
   }

   Component<FunctionBase> VectorFunctionBase::z() const
   {
      assert(getDimension() >= 3);
      return Component<FunctionBase>(*this, 2);
   }
}

