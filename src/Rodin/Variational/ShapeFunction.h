/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2024.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_SHAPEFUNCTION_H
#define RODIN_VARIATIONAL_SHAPEFUNCTION_H

#include "Rodin/Alert/Exception.h"
#include "Rodin/FormLanguage/Base.h"
#include "Rodin/FormLanguage/Traits.h"

#include "ForwardDecls.h"

#include "RangeType.h"
#include "RangeShape.h"
#include "FiniteElementSpace.h"

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

  template <class FES>
  struct IsTrialFunction<TrialFunction<FES>>
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

      using Parent = FormLanguage::Base;

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

      inline
      Derived& getDerived()
      {
        return static_cast<Derived&>(*this);
      }

      inline
      const Derived& getDerived() const
      {
        return static_cast<const Derived&>(*this);
      }

      /**
       * @brief Indicates whether the shape function is part of a %Trial or %Test
       * function expression.
       */
      inline
      constexpr
      ShapeFunctionSpaceType getSpaceType() const
      {
        return Space;
      }

      /**
       * @brief Gets the shape of the range space.
       * @note CRTP function to be overriden in the Derived class.
       */
      inline
      constexpr
      RangeShape getRangeShape() const
      {
        return static_cast<const Derived&>(*this).getRangeShape();
      }

      inline
      constexpr
      RangeType getRangeType() const
      {
        using R = typename FormLanguage::Traits<ShapeFunctionBase<Derived, FES, Space>>::RangeType;
        if constexpr (std::is_same_v<R, Boolean>)
        {
          return RangeType::Boolean;
        }
        else if constexpr (std::is_same_v<R, Integer>)
        {
          return RangeType::Integer;
        }
        else if constexpr (std::is_same_v<R, Real>)
        {
          return RangeType::Real;
        }
        else if constexpr (std::is_same_v<R, Math::Vector<Real>>)
        {
          return RangeType::Vector;
        }
        else if constexpr (std::is_same_v<R, Math::Matrix<Real>>)
        {
          return RangeType::Matrix;
        }
        else
        {
          assert(false);
          static_assert(Utility::DependentFalse<R>::Value);
        }
      }

      inline
      auto x() const
      {
        return Component(*this, 0);
      }

      inline
      auto y() const
      {
        return Component(*this, 1);
      }

      inline
      auto z() const
      {
        return Component(*this, 2);
      }

      inline
      constexpr
      auto T() const
      {
        return Transpose(*this);
      }

      /**
       * @brief Gets the operand in the shape function expression.
       * @note CRTP function to be overriden in the Derived class.
       */
      inline
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
      inline
      constexpr
      size_t getDOFs(const Geometry::Polytope& polytope) const
      {
        return static_cast<const Derived&>(*this).getDOFs(polytope);
      }

      inline
      const Geometry::Point& getPoint() const
      {
        return static_cast<const Derived&>(*this).getPoint();
      }

      inline
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
      inline
      constexpr
      auto getBasis(size_t local) const
      {
        return static_cast<const Derived&>(*this).getBasis(local);
      }

      /**
       * @brief Call operator to get an expression which yields the shape
       * function basis at the given point.
       *
       * Synonym to getBasis().
       */
      template <class ... Args>
      inline
      constexpr
      auto operator()(Args&&... args) const
      {
        return this->getBasis(std::forward<Args>(args)...);
      }

      /**
       * @brief Gets the finite element space to which the shape function
       * belongs to.
       */
      inline
      constexpr
      const FES& getFiniteElementSpace() const
      {
        return m_fes.get();
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

      using RangeType = typename FormLanguage::Traits<FESType>::RangeType;

      using Parent = ShapeFunctionBase<ShapeFunction<Derived, FESType, SpaceType>, FESType, SpaceType>;

      ShapeFunction() = delete;

      constexpr
      ShapeFunction(const FESType& fes)
        : Parent(fes)
      {}

      constexpr
      ShapeFunction(const ShapeFunction& other)
        : Parent(other)
      {}

      constexpr
      ShapeFunction(ShapeFunction&& other)
        : Parent(std::move(other))
      {}

      inline
      constexpr
      auto& emplace()
      {
        m_gf.emplace(this->getFiniteElementSpace());
        return *this;
      }

      inline
      constexpr
      RangeShape getRangeShape() const
      {
        return { this->getFiniteElementSpace().getVectorDimension(), 1 };
      }

      inline
      constexpr
      GridFunction<FES>& getSolution()
      {
        assert(m_gf.has_value());
        return m_gf.value();
      }

      inline
      constexpr
      const GridFunction<FES>& getSolution() const
      {
        assert(m_gf.has_value());
        return m_gf.value();
      }

      inline
      constexpr
      size_t getDOFs(const Geometry::Polytope& polytope) const
      {
        const size_t d = polytope.getDimension();
        const size_t i = polytope.getIndex();
        return this->getFiniteElementSpace().getFiniteElement(d, i).getCount();
      }

      inline
      const Geometry::Point& getPoint() const
      {
        assert(m_p.has_value());
        return m_p.value().get();
      }

      inline
      ShapeFunction& setPoint(const Geometry::Point& p)
      {
        m_p = p;
        return *this;
      }

      inline
      constexpr
      auto getBasis(size_t local) const
      {
        const auto& p = m_p.value().get();
        const size_t d = p.getPolytope().getDimension();
        const Index i = p.getPolytope().getIndex();
        const auto& fes = this->getFiniteElementSpace();
        const auto& fe = fes.getFiniteElement(d, i);
        if constexpr (std::is_same_v<RangeType, Real>)
        {
          return fes.getInverseMapping({ d, i }, fe.getBasis(local))(p);
        }
        else if constexpr (std::is_same_v<RangeType, Math::Vector<Real>>)
        {
          return this->object(fes.getInverseMapping({ d, i }, fe.getBasis(local))(p));
        }
        else
        {
          assert(false);
          return void();
        }
      }

      inline
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
      std::optional<GridFunction<FES>> m_gf;

      std::optional<std::reference_wrapper<const Geometry::Point>> m_p;
  };
}

#endif
