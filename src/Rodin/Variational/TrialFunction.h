#ifndef RODIN_VARIATIONAL_TRIALFUNCTION_H
#define RODIN_VARIATIONAL_TRIALFUNCTION_H

#include "Component.h"
#include "ShapeFunction.h"

namespace Rodin::FormLanguage
{
  template <class Solution, class FES>
  struct Traits<Variational::TrialFunction<Solution, FES>>
  {
    using FESType = FES;
    static constexpr Variational::ShapeFunctionSpaceType SpaceType = Variational::TrialSpace;

    using SolutionType = Solution;
  };
}

namespace Rodin::Variational
{
  template <class Data, class FES>
  class TrialFunction final
    : public ShapeFunction<TrialFunction<Data, FES>, FES, TrialSpace>
  {
    public:
      using FESType = FES;
      static constexpr ShapeFunctionSpaceType Space = TrialSpace;

      using DataType = Data;

      using SolutionType = GridFunction<FES, DataType>;

      using Parent =
        ShapeFunction<TrialFunction<DataType, FESType>, FESType, TrialSpace>;

      static_assert(std::is_base_of_v<FiniteElementSpaceBase, FES>);

      constexpr
      TrialFunction(const FES& fes)
        : Parent(fes),
          m_gf(fes)
      {}

      constexpr
      TrialFunction(const FES& fes, DataType& data)
        : Parent(fes),
          m_gf(fes, data)
      {}

      constexpr
      TrialFunction(const FES& fes, DataType&& data)
        : Parent(fes),
          m_gf(fes, std::move(data))
      {}

      constexpr
      TrialFunction(const TrialFunction& other)
        : Parent(other),
          m_gf(other.m_gf)
      {}

      constexpr
      TrialFunction(TrialFunction&& other)
        : Parent(std::move(other)),
          m_gf(std::move(other.m_gf))
      {}

      void operator=(const TrialFunction&) = delete;

      void operator=(TrialFunction&&) = delete;

      constexpr
      auto x() const
      {
        assert(this->getFiniteElementSpace().getVectorDimension() >= 1);
        return Component(*this, 0);
      }

      constexpr
      auto y() const
      {
        assert(this->getFiniteElementSpace().getVectorDimension() >= 2);
        return Component(*this, 1);
      }

      constexpr
      auto z() const
      {
        assert(this->getFiniteElementSpace().getVectorDimension() >= 3);
        return Component(*this, 2);
      }

      constexpr
      const TrialFunction& getLeaf() const
      {
        return *this;
      }

      constexpr
      SolutionType& getSolution()
      {
        return m_gf;
      }

      constexpr
      const SolutionType& getSolution() const
      {
        return m_gf;
      }

      TrialFunction* copy() const noexcept override
      {
        // No data is copied.
        return new TrialFunction(this->getFiniteElementSpace());
      }

    private:
      SolutionType m_gf;
  };

  template <class FES>
  TrialFunction(const FES&)
    -> TrialFunction<Math::Vector<typename FormLanguage::Traits<FES>::ScalarType>, FES>;

  template <class Data, class FES>
  TrialFunction(const FES&, Data&) -> TrialFunction<Data, FES>;

  template <class Data, class FES>
  TrialFunction(const FES&, Data&&) -> TrialFunction<Data, FES>;
}
#endif

