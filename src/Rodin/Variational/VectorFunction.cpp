#include "Component.h"

#include "VectorFunction.h"

namespace Rodin::Variational
{
   Component<VectorFunctionBase> VectorFunctionBase::operator()(int i) const
   {
      assert(0 <= i);
      assert(i < getDimension());
      return Component<VectorFunctionBase>(*this, i);
   }

   Component<VectorFunctionBase> VectorFunctionBase::x() const
   {
      assert(getDimension() >= 1);
      return Component<VectorFunctionBase>(*this, 0);
   }

   Component<VectorFunctionBase> VectorFunctionBase::y() const
   {
      assert(getDimension() >= 2);
      return Component<VectorFunctionBase>(*this, 1);
   }

   Component<VectorFunctionBase> VectorFunctionBase::z() const
   {
      assert(getDimension() >= 3);
      return Component<VectorFunctionBase>(*this, 2);
   }

  std::unique_ptr<mfem::VectorCoefficient> VectorFunctionBase::build() const
  {
     return std::unique_ptr<mfem::VectorCoefficient>(new Internal::ProxyVectorFunction(*this));
  }
}

