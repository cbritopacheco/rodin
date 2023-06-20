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

  template <class Derived, class FESType, ShapeFunctionSpaceType Space>
  class ShapeFunctionBase : public FormLanguage::Base
  {
    public:
      using Parent = FormLanguage::Base;
      using FES = FESType;

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
      constexpr
      ShapeFunctionSpaceType getSpaceType() const
      {
        return Space;
      }

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

      inline
      constexpr
      const auto& getLeaf() const
      {
        return static_cast<const Derived&>(*this).getLeaf();
      }

      inline
      constexpr
      size_t getDOFs(const Geometry::Polytope& element) const
      {
        return static_cast<const Derived&>(*this).getDOFs(element);
      }

      inline
      constexpr
      auto getTensorBasis(const Geometry::Point& p) const
      {
        return static_cast<const Derived&>(*this).getTensorBasis(p);
      }

      template <class ... Args>
      inline
      constexpr
      auto operator()(Args&&... args) const
      {
        return getTensorBasis(std::forward<Args>(args)...);
      }

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

  // /**
  // * @ingroup ShapeFunctionSpecializations
  // * @brief H1 ShapeFunction
  // */
  // template <class Derived, class ... Ps, ShapeFunctionSpaceType Space>
  // class ShapeFunction<Derived, H1<Scalar, Ps...>, Space>
  //   : public ShapeFunctionBase<ShapeFunction<Derived, H1<Scalar, Ps...>, Space>, H1<Scalar, Ps...>, Space>
  // {
  //   public:
  //     using FES = H1<Scalar, Ps...>;
  //     using Parent = ShapeFunctionBase<ShapeFunction<Derived, FES, Space>, FES, Space>;

  //     ShapeFunction() = delete;

  //     constexpr
  //     ShapeFunction(const FES& fes)
  //       : Parent(fes)
  //     {
  //       assert(fes.getVectorDimension() == 1);
  //     }

  //     constexpr
  //     ShapeFunction(const ShapeFunction& other)
  //       : Parent(other),
  //         m_gf(other.m_gf)
  //     {}

  //     constexpr
  //     ShapeFunction(ShapeFunction&& other)
  //       : Parent(std::move(other)),
  //         m_gf(std::move(other.m_gf))
  //     {}

  //     inline
  //     constexpr
  //     auto& emplace()
  //     {
  //       m_gf.emplace(this->getFiniteElementSpace());
  //       return *this;
  //     }

  //     inline
  //     constexpr
  //     RangeShape getRangeShape() const
  //     {
  //       return { 1, 1 };
  //     }

  //     inline
  //     constexpr
  //     GridFunction<FES>& getSolution()
  //     {
  //       assert(m_gf);
  //       return *m_gf;
  //     }

  //     inline
  //     constexpr
  //     const GridFunction<FES>& getSolution() const
  //     {
  //       assert(m_gf);
  //       return *m_gf;
  //     }

  //     inline
  //     constexpr
  //     size_t getDOFs(const Geometry::Polytope& element) const
  //     {
  //       return this->getFiniteElementSpace().getFiniteElement(element).getDOFs();
  //     }

  //     inline
  //     TensorBasis<Scalar> getTensorBasis(const Geometry::Point& p) const
  //     {
  //       const auto& fe = this->getFiniteElementSpace().getFiniteElement(p.getSimplex());
  //       return fe.getBasis(p.getCoordinates(Geometry::Point::Coordinates::Reference));
  //     }

  //     inline
  //     constexpr
  //     const auto& getLeaf() const
  //     {
  //       return static_cast<const Derived&>(*this).getLeaf();
  //     }

  //     virtual ShapeFunction* copy() const noexcept override
  //     {
  //       return static_cast<const Derived&>(*this).copy();
  //     }

  //   private:
  //     std::optional<GridFunction<FES>> m_gf;
  // };

  // /**
  // * @ingroup ShapeFunctionSpecializations
  // * @brief H1 ShapeFunction
  // */
  // template <class Derived, class ... Ps, ShapeFunctionSpaceType Space>
  // class ShapeFunction<Derived, H1<Math::Vector, Ps...>, Space>
  //   : public ShapeFunctionBase<ShapeFunction<Derived, H1<Math::Vector, Ps...>, Space>, H1<Math::Vector, Ps...>, Space>
  // {
  //   public:
  //     using FES = H1<Math::Vector, Ps...>;
  //     using Parent = ShapeFunctionBase<ShapeFunction<Derived, FES, Space>, FES, Space>;

  //     ShapeFunction() = delete;

  //     constexpr
  //     ShapeFunction(const FES& fes)
  //       : Parent(fes)
  //     {}

  //     constexpr
  //     ShapeFunction(const ShapeFunction& other)
  //       : Parent(other),
  //         m_gf(other.m_gf)
  //     {}

  //     constexpr
  //     ShapeFunction(ShapeFunction&& other)
  //       : Parent(std::move(other)),
  //         m_gf(std::move(other.m_gf))
  //     {}

  //     inline
  //     constexpr
  //     auto& emplace()
  //     {
  //       m_gf.emplace(this->getFiniteElementSpace());
  //       return *this;
  //     }

  //     inline
  //     constexpr
  //     GridFunction<FES>& getSolution()
  //     {
  //       assert(m_gf);
  //       return *m_gf;
  //     }

  //     inline
  //     constexpr
  //     const GridFunction<FES>& getSolution() const
  //     {
  //       assert(m_gf);
  //       return *m_gf;
  //     }

  //     inline
  //     constexpr
  //     RangeShape getRangeShape() const
  //     {
  //       return { this->getFiniteElementSpace().getVectorDimension(), 1 };
  //     }

  //     inline
  //     constexpr
  //     size_t getDOFs(const Geometry::Polytope& element) const
  //     {
  //       return this->getFiniteElementSpace().getFiniteElement(element).getDOFs();
  //     }

  //     inline
  //     TensorBasis<Math::Vector> getTensorBasis(const Geometry::Point& p) const
  //     {
  //       const auto& fe = this->getFiniteElementSpace().getFiniteElement(p.getSimplex());
  //       const Math::Vector& coords = p.getCoordinates(Geometry::Point::Coordinates::Reference);
  //       return fe.getBasis(coords).transpose();
  //     }

  //     inline
  //     constexpr
  //     const auto& getLeaf() const
  //     {
  //       return static_cast<const Derived&>(*this).getLeaf();
  //     }

  //     virtual ShapeFunction* copy() const noexcept override
  //     {
  //       return static_cast<const Derived&>(*this).copy();
  //     }

  //   private:
  //     std::optional<GridFunction<FES>> m_gf;
  // };

  /**
  * @ingroup ShapeFunctionSpecializations
  * @brief ShapeFunction
  */
  template <class Derived, class FESType, ShapeFunctionSpaceType Space>
  class ShapeFunction : public ShapeFunctionBase<ShapeFunction<Derived, FESType, Space>, FESType, Space>
  {
    public:
      using FES = FESType;
      using Parent = ShapeFunctionBase<ShapeFunction<Derived, FES, Space>, FES, Space>;

      ShapeFunction() = delete;

      constexpr
      ShapeFunction(const FES& fes)
        : Parent(fes)
      {
        assert(fes.getVectorDimension() == 1);
      }

      constexpr
      ShapeFunction(const ShapeFunction& other)
        : Parent(other),
          m_gf(other.m_gf)
      {}

      constexpr
      ShapeFunction(ShapeFunction&& other)
        : Parent(std::move(other)),
          m_gf(std::move(other.m_gf))
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
        return { 1, 1 };
      }

      inline
      constexpr
      GridFunction<FES>& getSolution()
      {
        assert(m_gf);
        return *m_gf;
      }

      inline
      constexpr
      const GridFunction<FES>& getSolution() const
      {
        assert(m_gf);
        return *m_gf;
      }

      inline
      constexpr
      size_t getDOFs(const Geometry::Polytope& element) const
      {
        const size_t d = element.getDimension();
        const size_t i = element.getIndex();
        return this->getFiniteElementSpace().getFiniteElement(d, i).getCount();
      }

      inline
      constexpr
      auto getTensorBasis(const Geometry::Point& p) const
      {
        const size_t d = p.getPolytope().getDimension();
        const Index i = p.getPolytope().getIndex();
        const auto& fe = this->getFiniteElementSpace().getFiniteElement(d, i);
        const auto& rc = p.getCoordinates(Geometry::Point::Coordinates::Reference);
        return TensorBasis(fe.getCount(),
            [&](size_t local){ return this->object(fe.getBasis(local)(rc)); });
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
