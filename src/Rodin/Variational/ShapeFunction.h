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
  namespace Internal
  {
    template <typename T, class ... Args>
    struct HasGetBasisMethod
    {
      template<typename U, typename = decltype(std::declval<U>().getBasis(std::declval<Args>()...))>
      static std::true_type Test(int);

      template<typename U>
      static std::false_type Test(...);

      using Type = decltype(Test<T>(0));
      static constexpr bool Value = Type::value;
    };

    template <typename T, typename... Args>
    struct HasGetBasisMethod<T, Args&...>
    {
      template <typename U, typename = decltype(std::declval<U>().getBasis(std::declval<Args&>()...))>
      static std::true_type Test(int);

      template <typename U>
      static std::false_type Test(...);

      using Type = decltype(Test<T>(0));
      static constexpr bool Value = Type::value;
    };

    template <typename T, class... Args>
    struct HasGetBasisMethodR
    {
        template<typename U, typename = decltype(std::declval<U>().getBasis(std::declval<Args>()...))>
        static auto Test(int) ->
          decltype(std::is_same<typename std::invoke_result<decltype(&U::getBasis)(U, Args...)>::type, T>::value, std::true_type{});

        template<typename U>
        static std::false_type Test(...);

        using Type = decltype(Test<T>(0));
        static constexpr bool Value = Type::value;
    };
  }

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

      Derived& getDerived()
      {
        return static_cast<Derived&>(*this);
      }

      const Derived& getDerived() const
      {
        return static_cast<const Derived&>(*this);
      }

      /**
       * @brief Indicates whether the shape function is part of a %Trial or %Test
       * function expression.
       */
      constexpr
      ShapeFunctionSpaceType getSpaceType() const
      {
        return Space;
      }

      /**
       * @brief Gets the shape of the range space.
       * @note CRTP function to be overriden in the Derived class.
       */
      constexpr
      RangeShape getRangeShape() const
      {
        return static_cast<const Derived&>(*this).getRangeShape();
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
      auto getBasis(size_t local) const
      {
        return static_cast<const Derived&>(*this).getBasis(local);
      }

      template <class T>
      constexpr
      void getBasis(T& basis, size_t local) const
      {
        if constexpr (Internal::HasGetBasisMethod<Derived, T&, size_t>::Value)
        {
          static_cast<const Derived&>(*this).getBasis(basis, local);
        }
        else
        {
          basis = getBasis(local);
        }
      }

      /**
       * @brief Call operator to get an expression which yields the shape
       * function basis at the given point.
       *
       * Synonym to getBasis(size_t).
       */
      constexpr
      auto operator()(size_t local) const
      {
        return getBasis(local);
      }

      template <class T>
      constexpr
      void operator()(T& res, size_t local) const
      {
        getBasis(res, local);
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

      constexpr
      auto& emplace()
      {
        m_gf.emplace(this->getFiniteElementSpace());
        return *this;
      }

      constexpr
      RangeShape getRangeShape() const
      {
        return { this->getFiniteElementSpace().getVectorDimension(), 1 };
      }

      constexpr
      GridFunction<FES>& getSolution()
      {
        assert(m_gf.has_value());
        return m_gf.value();
      }

      constexpr
      const GridFunction<FES>& getSolution() const
      {
        assert(m_gf.has_value());
        return m_gf.value();
      }

      constexpr
      size_t getDOFs(const Geometry::Polytope& polytope) const
      {
        const size_t d = polytope.getDimension();
        const size_t i = polytope.getIndex();
        return this->getFiniteElementSpace().getFiniteElement(d, i).getCount();
      }

      const Geometry::Point& getPoint() const
      {
        assert(m_p.has_value());
        return m_p.value().get();
      }

      ShapeFunction& setPoint(const Geometry::Point& p)
      {
        m_p = p;
        return *this;
      }

      constexpr
      auto getBasis(size_t local) const
      {
        const auto& p = m_p.value().get();
        const size_t d = p.getPolytope().getDimension();
        const Index i = p.getPolytope().getIndex();
        const auto& fes = this->getFiniteElementSpace();
        const auto& fe = fes.getFiniteElement(d, i);
        return this->object(fes.getInverseMapping({ d, i }, fe.getBasis(local))(p));
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
      std::optional<GridFunction<FES>> m_gf;

      std::optional<std::reference_wrapper<const Geometry::Point>> m_p;
  };
}

#endif
