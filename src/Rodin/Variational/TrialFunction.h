#ifndef RODIN_VARIATIONAL_TRIALFUNCTION_H
#define RODIN_VARIATIONAL_TRIALFUNCTION_H

#include "Component.h"
#include "ShapeFunction.h"

namespace Rodin::FormLanguage
{
  template <class FESType>
  struct Traits<Variational::TrialFunction<FESType>>
  {
    using FES = FESType;
    static constexpr Variational::ShapeFunctionSpaceType SpaceType = Variational::TrialSpace;
  };
}

namespace Rodin::Variational
{
  template <class FESType>
  class TrialFunction final
    : public ShapeFunction<TrialFunction<FESType>, FESType, TrialSpace>
  {
    public:
      using FES = FESType;
      static constexpr ShapeFunctionSpaceType Space = TrialSpace;

      using Parent = ShapeFunction<TrialFunction<FESType>, FESType, TrialSpace>;

      static_assert(std::is_base_of_v<FiniteElementSpaceBase, FES>,
          "FES is not a finite element space.");

      constexpr
      TrialFunction(const FES& fes)
        : Parent(fes)
      {}

      constexpr
      TrialFunction(const TrialFunction& other)
        : Parent(other)
      {}

      constexpr
      TrialFunction(TrialFunction&& other)
        : Parent(std::move(other))
      {}

      void operator=(const TrialFunction&) = delete;

      void operator=(TrialFunction&&) = delete;

      inline
      constexpr
      auto x() const
      {
        assert(this->getFiniteElementSpace().getVectorDimension() >= 1);
        return Component(*this, 0);
      }

      inline
      constexpr
      auto y() const
      {
        assert(this->getFiniteElementSpace().getVectorDimension() >= 2);
        return Component(*this, 1);
      }

      inline
      constexpr
      auto z() const
      {
        assert(this->getFiniteElementSpace().getVectorDimension() >= 3);
        return Component(*this, 2);
      }

      inline
      constexpr
      const TrialFunction& getLeaf() const
      {
        return *this;
      }

      inline
      TrialFunction* copy() const noexcept override
      {
        return new TrialFunction(*this);
      }
  };

  template <class FES>
  TrialFunction(const FES&) -> TrialFunction<FES>;
}
#endif

