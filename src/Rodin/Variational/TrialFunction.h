#ifndef RODIN_VARIATIONAL_TRIALFUNCTION_H
#define RODIN_VARIATIONAL_TRIALFUNCTION_H

#include <functional>

#include "Rodin/Variational/ForwardDecls.h"
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
  template <class Solution, class FES>
  class TrialFunctionReference
    : public ShapeFunction<TrialFunctionReference<Solution, FES>, FES, ShapeFunctionSpaceType::Trial>
  {
    public:
      using FESType = FES;
      static constexpr ShapeFunctionSpaceType SpaceType = ShapeFunctionSpaceType::Trial;

      using Parent =
        ShapeFunction<TrialFunctionReference<Solution, FESType>, FESType, SpaceType>;

      explicit
      TrialFunctionReference(const TrialFunction<Solution, FESType>& ref, const FESType& fes)
        : Parent(fes),
          m_ref(std::cref(ref))
      {}

      TrialFunctionReference(const TrialFunctionReference& other)
        : Parent(other),
          m_ref(other.m_ref)
      {}

      TrialFunctionReference(TrialFunctionReference&& other)
        : Parent(std::move(other)),
          m_ref(std::move(other.m_ref))
      {}

      TrialFunctionReference& setPoint(const Geometry::Point& p)
      {
        m_ref.get().setPoint(p);
        return *this;
      }

      constexpr
      auto x() const
      {
        return m_ref.get().x();
      }

      constexpr
      auto y() const
      {
        return m_ref.get().y();
      }

      constexpr
      auto z() const
      {
        return m_ref.get().z();
      }

      constexpr
      const auto& getLeaf() const
      {
        return m_ref.get().getLeaf();
      }

      constexpr
      auto& getSolution()
      {
        return m_ref.get().getSolution();
      }

      constexpr
      const auto& getSolution() const
      {
        return m_ref.get().getSolution();
      }

      TrialFunctionReference* copy() const noexcept final override
      {
        return new TrialFunctionReference(*this);
      }

    private:
      std::reference_wrapper<const TrialFunction<Solution, FESType>> m_ref;
  };

  template <class Solution, class FES>
  class TrialFunction : public TrialFunctionReference<Solution, FES>
  {
    public:
      using FESType = FES;
      static constexpr ShapeFunctionSpaceType Space = TrialSpace;

      using SolutionType = Solution;

      using Parent = TrialFunctionReference<SolutionType, FESType>;

      static_assert(std::is_base_of_v<FiniteElementSpaceBase, FES>);

      constexpr
      TrialFunction(const FES& fes)
        : Parent(*this, fes),
          m_gf(fes)
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

    private:
      SolutionType m_gf;
  };

  template <class FES>
  TrialFunction(const FES&)
    -> TrialFunction<
        GridFunction<FES, Math::Vector<
          typename FormLanguage::Traits<FES>::ScalarType>>, FES>;
}
#endif

