#ifndef RODIN_VARIATIONAL_TRIALFUNCTION_H
#define RODIN_VARIATIONAL_TRIALFUNCTION_H

#include "Component.h"
#include "GridFunction.h"
#include "ShapeFunction.h"

namespace Rodin::Variational
{
  template <class FES>
  class TrialFunction : public ShapeFunction<FES, TrialSpace>
  {
    public:
      constexpr
      TrialFunction(FES& fes)
        : ShapeFunction<FES, TrialSpace>(fes)
      {}

      constexpr
      TrialFunction(const TrialFunction& other)
        : ShapeFunction<FES, TrialSpace>(other)
      {}

      constexpr
      TrialFunction(TrialFunction&& other)
        : ShapeFunction<FES, TrialSpace>(std::move(other))
      {}

      void operator=(const TrialFunction&) = delete;

      void operator=(TrialFunction&&) = delete;

      constexpr
      TrialFunction& emplace()
      {
        m_gf.emplace(this->getFiniteElementSpace());
        return *this;
      }

      constexpr
      GridFunction<FES>& getSolution()
      {
        assert(m_gf);
        return *m_gf;
      }

      constexpr
      const GridFunction<FES>& getSolution() const
      {
        assert(m_gf);
        return *m_gf;
      }

      constexpr
      Component<TrialFunction<FES>> x() const
      {
        assert(this->getFiniteElementSpace().getVectorDimension() >= 1);
        return Component<TrialFunction<FES>>(*this, 0);
      }

      constexpr
      Component<TrialFunction<FES>> y() const
      {
        assert(this->getFiniteElementSpace().getVectorDimension() >= 2);
        return Component<TrialFunction<FES>>(*this, 1);
      }

      constexpr
      Component<TrialFunction<FES>> z() const
      {
        assert(this->getFiniteElementSpace().getVectorDimension() >= 3);
        return Component<TrialFunction<FES>>(*this, 2);
      }

      const TrialFunction& getLeaf() const override
      {
        return *this;
      }

      TrialFunction* copy() const noexcept override
      {
        return new TrialFunction(*this);
      }
    private:
      std::optional<GridFunction<FES>> m_gf;
  };

  template <class FES>
  TrialFunction(FES&) -> TrialFunction<FES>;
}
#endif

