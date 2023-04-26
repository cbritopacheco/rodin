#ifndef RODIN_VARIATIONAL_SHAPEFUNCTION_H
#define RODIN_VARIATIONAL_SHAPEFUNCTION_H

#include <mfem.hpp>

#include "Rodin/Alert/Exception.h"
#include "Rodin/FormLanguage/Base.h"
#include "Rodin/FormLanguage/Traits.h"

#include "ForwardDecls.h"

#include "H1.h"
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
  // * @brief L2 ShapeFunction
  // */
  // template <class Derived, class ... Ts, ShapeFunctionSpaceType Space>
  // class ShapeFunction<Derived, L2<Ts...>, Space>
  //   : public ShapeFunctionBase<ShapeFunction<Derived, L2<Ts...>, Space>, L2<Ts...>, Space>
  // {
  //   public:
  //     using FES = L2<Ts...>;
  //     using Parent = ShapeFunctionBase<ShapeFunction<Derived, L2<Ts...>, Space>, L2<Ts...>, Space>;

  //     constexpr
  //     ShapeFunction(FES& fes)
  //       : Parent(fes)
  //     {}

  //     constexpr
  //     ShapeFunction(const ShapeFunction& other)
  //       : Parent(other)
  //     {}

  //     constexpr
  //     ShapeFunction(ShapeFunction&& other)
  //       : Parent(std::move(other))
  //     {}

  //     inline
  //     constexpr
  //     auto& emplace()
  //     {
  //       m_gf.emplace(this->getFiniteElementSpace().get());
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
  //     auto x() const
  //     {
  //       assert(this->getFiniteElementSpace().getVectorDimension() >= 1);
  //       return static_cast<const Derived&>(*this).x();
  //     }

  //     inline
  //     constexpr
  //     auto y() const
  //     {
  //       assert(this->getFiniteElementSpace().getVectorDimension() >= 2);
  //       return static_cast<const Derived&>(*this).y();
  //     }

  //     inline
  //     constexpr
  //     auto z() const
  //     {
  //       assert(this->getFiniteElementSpace().getVectorDimension() >= 3);
  //       return static_cast<const Derived&>(*this).z();
  //     }

  //     inline
  //     constexpr
  //     size_t getDOFs(const Geometry::Simplex& element) const
  //     {
  //       const auto& fe = this->getFiniteElementSpace().getFiniteElement(element);
  //       return fe.GetDof() * this->getFiniteElementSpace().getVectorDimension();
  //     }

  //     inline
  //     constexpr
  //     const auto& getLeaf() const
  //     {
  //       return static_cast<const Derived&>(*this).getLeaf();
  //     }

  //     inline
  //     constexpr
  //     auto getTensorBasis(const Geometry::Point& p) const
  //     {
  //       return static_cast<const Derived&>(*this).getTensorBasis(p);
  //     }

  //     virtual ShapeFunction* copy() const noexcept override
  //     {
  //       return static_cast<const Derived&>(*this).copy();
  //     }

  //     // void getTensorBasis(
  //     //    DenseBasisOperator& op,
  //     //    ShapeComputator& compute,
  //     //    const Geometry::Point& point,
  //     //    const Geometry::Element& element) const override
  //     // {
  //     //   const auto& shape =
  //     //    compute.getPhysicalShape(
  //     //       getFiniteElementSpace().getFiniteElement(element),
  //     //       element.getTransformation(),
  //     //       element.getTransformation().GetIntPoint());
  //     //   const int n = shape.Size();
  //     //   const int vdim = getFiniteElementSpace().getVectorDimension();
  //     //   op.setSize(vdim, 1, vdim * n);
  //     //   op = 0.0;
  //     //   for (int i = 0; i < vdim; i++)
  //     //    for (int j = 0; j < n; j++)
  //     //     op(i, 0, j + i * n) = shape(j);
  //     // }

  //   private:
  //     std::optional<GridFunction<FES>> m_gf;
  // };

  /**
  * @ingroup ShapeFunctionSpecializations
  * @brief H1 ShapeFunction
  */
  template <class Derived, class ... Ps, ShapeFunctionSpaceType Space>
  class ShapeFunction<Derived, H1<Scalar, Ps...>, Space>
    : public ShapeFunctionBase<ShapeFunction<Derived, H1<Scalar, Ps...>, Space>, H1<Scalar, Ps...>, Space>
  {
    public:
      using FES = H1<Scalar, Ps...>;
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
        return this->getFiniteElementSpace().getFiniteElement(element).getDOFs();
      }

      inline
      TensorBasis<Scalar> getTensorBasis(const Geometry::Point& p) const
      {
        const auto& fe = this->getFiniteElementSpace().getFiniteElement(p.getSimplex());
        return fe.getBasis(p.getCoordinates(Geometry::Point::Coordinates::Reference));
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

  /**
  * @ingroup ShapeFunctionSpecializations
  * @brief H1 ShapeFunction
  */
  template <class Derived, class ... Ps, ShapeFunctionSpaceType Space>
  class ShapeFunction<Derived, H1<Math::Vector, Ps...>, Space>
    : public ShapeFunctionBase<ShapeFunction<Derived, H1<Math::Vector, Ps...>, Space>, H1<Math::Vector, Ps...>, Space>
  {
    public:
      using FES = H1<Math::Vector, Ps...>;
      using Parent = ShapeFunctionBase<ShapeFunction<Derived, FES, Space>, FES, Space>;

      ShapeFunction() = delete;

      constexpr
      ShapeFunction(const FES& fes)
        : Parent(fes)
      {}

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
      RangeShape getRangeShape() const
      {
        return { this->getFiniteElementSpace().getVectorDimension(), 1 };
      }

      inline
      constexpr
      size_t getDOFs(const Geometry::Polytope& element) const
      {
        return this->getFiniteElementSpace().getFiniteElement(element).getDOFs();
      }

      inline
      TensorBasis<Math::Vector> getTensorBasis(const Geometry::Point& p) const
      {
        const auto& fe = this->getFiniteElementSpace().getFiniteElement(p.getSimplex());
        const Math::Vector& coords = p.getCoordinates(Geometry::Point::Coordinates::Reference);
        return fe.getBasis(coords).transpose();
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
