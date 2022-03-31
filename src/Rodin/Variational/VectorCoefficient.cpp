#include "Component.h"

#include "VectorCoefficient.h"

namespace Rodin::Variational
{
   Component<VectorCoefficientBase> VectorCoefficientBase::operator()(int i) const
   {
      assert(0 <= i);
      assert(i < getDimension());
      return Component<VectorCoefficientBase>(*this, i);
   }

   Component<VectorCoefficientBase> VectorCoefficientBase::x() const
   {
      assert(getDimension() >= 1);
      return Component<VectorCoefficientBase>(*this, 0);
   }

   Component<VectorCoefficientBase> VectorCoefficientBase::y() const
   {
      assert(getDimension() >= 2);
      return Component<VectorCoefficientBase>(*this, 1);
   }

   Component<VectorCoefficientBase> VectorCoefficientBase::z() const
   {
      assert(getDimension() >= 3);
      return Component<VectorCoefficientBase>(*this, 2);
   }

  std::unique_ptr<mfem::VectorCoefficient> VectorCoefficientBase::build() const
  {
     return std::unique_ptr<mfem::VectorCoefficient>(new Internal::ProxyVectorCoefficient(*this));
  }
}

