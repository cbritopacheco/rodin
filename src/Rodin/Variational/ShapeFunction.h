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
#include "TensorBasis.h"
#include "FiniteElementSpace.h"

namespace Rodin::FormLanguage
{
  template <class Derived, class FESType, Variational::ShapeFunctionSpaceType SpaceType>
  struct Traits<Variational::ShapeFunctionBase<Derived, FESType, SpaceType>>
  {
    using FES = FESType;

    static constexpr const Variational::ShapeFunctionSpaceType Space = SpaceType;

    using ResultType =
      typename ResultOf<Variational::ShapeFunctionBase<Derived, FES, SpaceType>>::Type;

    using RangeType =
      typename RangeOf<Variational::ShapeFunctionBase<Derived, FES, Space>>::Type;

  };

  template <class Derived, class FESType, Variational::ShapeFunctionSpaceType SpaceType>
  struct Traits<Variational::ShapeFunction<Derived, FESType, SpaceType>>
  {
    using FES = FESType;
    static constexpr const Variational::ShapeFunctionSpaceType Space = SpaceType;
  };
}

namespace Rodin::Variational
{
  /**
  * @defgroup ShapeFunctionSpecializations ShapeFunction Template Specializations
  * @brief Template specializations of the ShapeFunction class.
  * @see ShapeFunction
  */

  template <ShapeFunctionSpaceType Space>
  struct DualSpaceType;

  template <>
  struct DualSpaceType<TrialSpace>
  {
    static constexpr ShapeFunctionSpaceType Value = ShapeFunctionSpaceType::Test;
  };

  template <>
  struct DualSpaceType<TestSpace>
  {
    static constexpr ShapeFunctionSpaceType Value = ShapeFunctionSpaceType::Trial;
  };

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
    class FESType = typename FormLanguage::Traits<Derived>::FES,
    ShapeFunctionSpaceType SpaceType = FormLanguage::Traits<Derived>::Space>
  class ShapeFunctionBase : public FormLanguage::Base
  {
    public:
      using FES = FESType;
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
        else if constexpr (std::is_same_v<R, Scalar>)
        {
          return RangeType::Scalar;
        }
        else if constexpr (std::is_same_v<R, Math::Vector>)
        {
          return RangeType::Vector;
        }
        else if constexpr (std::is_same_v<R, Math::Matrix>)
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

      /**
       * @brief Gets an expression which yields the shape function basis at the
       * given point.
       * @param[in] p Point where the shape function basis will be calculated
       * @note CRTP function to be overriden in the Derived class.
       */
      inline
      constexpr
      auto getTensorBasis(const Geometry::Point& p) const
      {
        return static_cast<const Derived&>(*this).getTensorBasis(p);
      }

      /**
       * @brief Call operator to get an expression which yields the shape
       * function basis at the given point.
       *
       * Synonym to getTensorBasis().
       */
      template <class ... Args>
      inline
      constexpr
      auto operator()(Args&&... args) const
      {
        return this->getTensorBasis(std::forward<Args>(args)...);
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
  template <class Derived, class FESType, ShapeFunctionSpaceType SpaceType>
  class ShapeFunction
    : public ShapeFunctionBase<ShapeFunction<Derived, FESType, SpaceType>, FESType, SpaceType>
  {
    public:
      using FES = FESType;

      static constexpr ShapeFunctionSpaceType Space = SpaceType;

      using Parent = ShapeFunctionBase<ShapeFunction<Derived, FES, Space>, FES, Space>;

      ShapeFunction() = delete;

      constexpr
      ShapeFunction(const FES& fes)
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
      auto& getSolution()
      {
        assert(m_gf);
        return *m_gf;
      }

      inline
      constexpr
      const auto& getSolution() const
      {
        assert(m_gf);
        return *m_gf;
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
      constexpr
      auto getTensorBasis(const Geometry::Point& p) const
      {
        using RangeType = typename FES::RangeType;
        const size_t d = p.getPolytope().getDimension();
        const Index i = p.getPolytope().getIndex();
        const auto& fes = this->getFiniteElementSpace();
        const auto& fe = fes.getFiniteElement(d, i);
        if constexpr (std::is_same_v<RangeType, Scalar>)
        {
          return TensorBasis(fe.getCount(),
              [&](size_t local)
              {
                return fes.getInverseMapping({ d, i }, fe.getBasis(local))(p);
              });
        }
        else if constexpr (std::is_same_v<RangeType, Math::Vector>)
        {
          return TensorBasis(fe.getCount(),
              [&](size_t local)
              {
                return this->object(fes.getInverseMapping({ d, i }, fe.getBasis(local))(p));
              });
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
  };
}

#endif
