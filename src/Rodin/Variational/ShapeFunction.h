/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2024.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_SHAPEFUNCTION_H
#define RODIN_VARIATIONAL_SHAPEFUNCTION_H

#include "ForwardDecls.h"

#include "Rodin/Geometry/Point.h"
#include "Rodin/FormLanguage/Base.h"
#include "Rodin/FormLanguage/Traits.h"
#include "Rodin/Variational/Traits.h"

namespace Rodin::FormLanguage
{
  template <class Derived, class FES, Variational::ShapeFunctionSpaceType Space>
  struct Traits<Variational::ShapeFunctionBase<Derived, FES, Space>>
  {
    using DerivedType = Derived;

    using FESType = FES;
    static constexpr const Variational::ShapeFunctionSpaceType SpaceType = Space;

    using ResultType =
      typename ResultOf<Variational::ShapeFunctionBase<Derived, FES, SpaceType>>::Type;

    using RangeType =
      typename RangeOf<Variational::ShapeFunctionBase<Derived, FES, SpaceType>>::Type;

    using ScalarType = typename FormLanguage::Traits<RangeType>::ScalarType;
  };

  template <class Derived, class FES, Variational::ShapeFunctionSpaceType Space>
  struct Traits<Variational::ShapeFunction<Derived, FES, Space>>
  {
    using DerivedType = Derived;

    using FESType = FES;
    static constexpr const Variational::ShapeFunctionSpaceType SpaceType = Space;

    using ResultType =
      typename ResultOf<
        Variational::ShapeFunctionBase<
          Variational::ShapeFunction<Derived, FES, SpaceType>, FES, SpaceType>>::Type;

    using RangeType =
      typename RangeOf<
        Variational::ShapeFunctionBase<
          Variational::ShapeFunction<Derived, FES, SpaceType>, FES, SpaceType>>::Type;

    using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;
  };
}

namespace Rodin::Variational
{
  /**
  * @defgroup ShapeFunctionSpecializations ShapeFunction Template Specializations
  * @brief Template specializations of the ShapeFunction class.
  * @see ShapeFunction
  */

  template <class T>
  struct IsTrialFunction
  {
    static constexpr Boolean Value = false;
  };

  template <class Solution, class FES>
  struct IsTrialFunction<TrialFunction<Solution, FES>>
  {
    static constexpr Boolean Value = true;
  };

  template <class T>
  struct IsTestFunction
  {
    static constexpr Boolean Value = false;
  };

  template <class FES>
  struct IsTestFunction<TestFunction<FES>>
  {
    static constexpr Boolean Value = true;
  };

  template <
    class Derived,
    class FES = typename FormLanguage::Traits<Derived>::FESType,
    ShapeFunctionSpaceType SpaceType = FormLanguage::Traits<Derived>::SpaceType>
  class ShapeFunctionBase : public FormLanguage::Base
  {
    public:
      using FESType = FES;
      static constexpr ShapeFunctionSpaceType Space = SpaceType;

      using ScalarType =
        typename FormLanguage::Traits<FESType>::ScalarType;

      using Parent =
        FormLanguage::Base;

      constexpr
      ShapeFunctionBase(const FES& fes)
        : m_fes(fes)
      {}

      constexpr
      ShapeFunctionBase(const ShapeFunctionBase& other)
        : Parent(other),
          m_fes(other.m_fes)
      {}

      constexpr
      ShapeFunctionBase(ShapeFunctionBase&& other)
        : Parent(std::move(other)),
          m_fes(std::move(other.m_fes))
      {}

      /**
       * @brief Indicates whether the shape function is part of a %Trial or %Test
       * function expression.
       */
      constexpr
      ShapeFunctionSpaceType getSpaceType() const
      {
        return Space;
      }

      auto x() const
      {
        return Component(*this, 0);
      }

      auto y() const
      {
        return Component(*this, 1);
      }

      auto z() const
      {
        return Component(*this, 2);
      }

      constexpr
      auto T() const
      {
        return Transpose(*this);
      }

      /**
       * @brief Gets the operand in the shape function expression.
       * @note CRTP function to be overriden in the Derived class.
       */
      constexpr
      const auto& getLeaf() const
      {
        return static_cast<const Derived&>(*this).getLeaf();
      }

      /**
       * @brief Gets the number of degrees of freedom for the given polytope.
       * @param[in] polytope Polytope
       * @note CRTP function to be overriden in the Derived class.
       */
      constexpr
      size_t getDOFs(const Geometry::Polytope& polytope) const
      {
        return static_cast<const Derived&>(*this).getDOFs(polytope);
      }

      constexpr
      const Geometry::Point& getPoint() const
      {
        return static_cast<const Derived&>(*this).getPoint();
      }

      constexpr
      Derived& setPoint(const Geometry::Point& p)
      {
        return static_cast<Derived&>(*this).setPoint(p);
      }

      /**
       * @brief Gets an expression which yields the shape function basis at the
       * given point.
       * @param[in] p Point where the shape function basis will be calculated
       * @note CRTP function to be overriden in the Derived class.
       */
      constexpr
      decltype(auto) getBasis(size_t local) const
      {
        return static_cast<const Derived&>(*this).getBasis(local);
      }

      /**
       * @brief Call operator to get an expression which yields the shape
       * function basis at the given point.
       *
       * Synonym to getBasis(size_t).
       */
      constexpr
      decltype(auto) operator()(size_t local) const
      {
        return static_cast<const Derived&>(*this).getBasis(local);
      }

      /**
       * @brief Gets the finite element space to which the shape function
       * belongs to.
       */
      constexpr
      const FES& getFiniteElementSpace() const
      {
        return m_fes.get();
      }

      Derived& getDerived() noexcept
      {
        return static_cast<Derived&>(*this);
      }

      const Derived& getDerived() const noexcept
      {
        return static_cast<const Derived&>(*this);
      }

      virtual ShapeFunctionBase* copy() const noexcept override
      {
        return static_cast<const Derived&>(*this).copy();
      }

    private:
      std::reference_wrapper<const FES> m_fes;
  };

  /**
  * @ingroup ShapeFunctionSpecializations
  * @brief ShapeFunction
  */
  template <class Derived, class FES, ShapeFunctionSpaceType Space>
  class ShapeFunction
    : public ShapeFunctionBase<ShapeFunction<Derived, FES, Space>, FES, Space>
  {
    public:
      using FESType = FES;
      static constexpr ShapeFunctionSpaceType SpaceType = Space;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using RangeType = typename FormLanguage::Traits<FESType>::RangeType;

      using Parent =
        ShapeFunctionBase<
          ShapeFunction<Derived, FESType, SpaceType>, FESType, SpaceType>;

      ShapeFunction() = delete;

      constexpr
      ShapeFunction(const FESType& fes)
        : Parent(fes),
          m_p(nullptr)
      {}

      constexpr
      ShapeFunction(const ShapeFunction& other)
        : Parent(other),
          m_basis(other.m_basis),
          m_p(nullptr)
      {}

      constexpr
      ShapeFunction(ShapeFunction&& other)
        : Parent(std::move(other)),
          m_basis(std::move(other.m_basis)),
          m_p(std::exchange(other.m_p, nullptr))
      {}

      constexpr
      size_t getDOFs(const Geometry::Polytope& polytope) const
      {
        const size_t d = polytope.getDimension();
        const size_t i = polytope.getIndex();
        return this->getFiniteElementSpace().getFiniteElement(d, i).getCount();
      }

      constexpr
      const Geometry::Point& getPoint() const
      {
        assert(m_p);
        return *m_p;
      }

      ShapeFunction& setPoint(const Geometry::Point& p)
      {
        if (m_p == &p)
          return *this;
        m_p = &p;
        const auto& polytope = p.getPolytope();
        const size_t d = polytope.getDimension();
        const Index i = polytope.getIndex();
        const auto& fes = this->getFiniteElementSpace();
        const auto& fe = fes.getFiniteElement(d, i);
        const size_t count = fe.getCount();
        m_basis.resize(count);
        for (size_t local = 0; local < count; local++)
          m_basis[local] = fes.getPushforward({ d, i }, fe.getBasis(local))(p);
        return *this;
      }

      constexpr
      const RangeType& getBasis(size_t local) const
      {
        assert(local < m_basis.size());
        return m_basis[local];
      }

      constexpr
      const auto& getLeaf() const
      {
        return static_cast<const Derived&>(*this).getLeaf();
      }

      virtual ShapeFunction* copy() const noexcept override
      {
        return static_cast<const Derived&>(*this).copy();
      }

    private:
      std::vector<RangeType> m_basis;

      const Geometry::Point* m_p;
  };
}

#endif
