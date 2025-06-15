#ifndef RODIN_VARIATIONAL_TRIALFUNCTION_H
#define RODIN_VARIATIONAL_TRIALFUNCTION_H

#include "Component.h"
#include "ShapeFunction.h"

namespace Rodin::FormLanguage
{
  template <class FES, class Solution>
  struct Traits<Variational::TrialFunction<FES, Solution>>
  {
    using FESType = FES;
    static constexpr Variational::ShapeFunctionSpaceType SpaceType = Variational::TrialSpace;
  };
}

namespace Rodin::Variational
{
  template <class Solution, class FES>
  class TrialFunction final
    : public ShapeFunction<TrialFunction<FES, Solution>, FES, TrialSpace>
  {
    public:
      using FESType = FES;
      static constexpr ShapeFunctionSpaceType Space = TrialSpace;

      using SolutionType = Solution;

      using Parent =
        ShapeFunction<TrialFunction<FESType, Solution>, FESType, TrialSpace>;

      static_assert(std::is_base_of_v<FiniteElementSpaceBase, FES>,
          "FES is not a finite element space.");

      constexpr
      TrialFunction(const FES& fes, SolutionType& gf)
        : Parent(fes),
          m_gf(std::ref(gf))
      {}

      constexpr
      TrialFunction(const FES& fes, SolutionType&& gf)
        : Parent(fes),
          m_gf(std::move(gf))
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
        assert(m_gf.has_value());
        return m_gf.value();
      }

      constexpr
      const SolutionType& getSolution() const
      {
        assert(m_gf.has_value());
        return m_gf.value();
      }

      constexpr
      auto& emplace()
      {
        m_gf.emplace(SolutionType(this->getFiniteElementSpace()));
        return *this;
      }

      TrialFunction* copy() const noexcept override
      {
        return new TrialFunction(*this);
      }

    private:
      std::optional<SolutionType> m_gf;
  };

  template <class Solution, class FES>
  TrialFunction(const FES&, Solution&)
    -> TrialFunction<Solution, FES>;

  template <class Solution, class FES>
  TrialFunction(const FES&, Solution&&)
    -> TrialFunction<Solution, FES>;

  template <class FES>
  TrialFunction(const FES&)
    -> TrialFunction<
        GridFunction<FES, Math::Vector<typename FormLanguage::Traits<FES>::ScalarType>>,
        FES>;
}
#endif

