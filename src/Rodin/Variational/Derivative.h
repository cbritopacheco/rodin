#ifndef RODIN_VARIATIONAL_DERIVATIVE_H
#define RODIN_VARIATIONAL_DERIVATIVE_H

#include <cassert>
#include <cstdlib>

#include "ForwardDecls.h"
#include "FiniteElementSpace.h"
#include "GridFunction.h"
#include "RealFunction.h"

namespace Rodin::FormLanguage
{
  template <class NestedDerived, class FES, Variational::ShapeFunctionSpaceType Space>
  struct Traits<Variational::Derivative<Variational::ShapeFunction<NestedDerived, FES, Space>>>
  {
    static constexpr Variational::ShapeFunctionSpaceType SpaceType = Space;
    using FESType = FES;
    using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;
    using OperandType = Variational::ShapeFunction<NestedDerived, FESType, Space>;
  };
}

namespace Rodin::Variational
{
  /**
   * @brief Base class for Grad classes.
   */
  template <class Operand, class Derived>
  class DerivativeBase;

  /**
   * @ingroup GradSpecializations
   */
  template <class FES, class Derived>
  class DerivativeBase<GridFunction<FES>, Derived>
    : public ScalarFunctionBase<
        typename FormLanguage::Traits<FES>::ScalarType, DerivativeBase<GridFunction<FES>, Derived>>
  {
    public:
      using FESType = FES;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using OperandType = GridFunction<FESType>;

      using Parent = ScalarFunctionBase<ScalarType, DerivativeBase<OperandType, Derived>>;

      /**
       * @brief Constructs the gradient of a @f$ \mathbb{P}_1 @f$ function @f$
       * u @f$.
       * @param[in] u P1 GridFunction
       */
      DerivativeBase(const OperandType& u)
        : m_u(u)
      {
        assert(u.getFiniteElementSpace().getVectorDimension() == 1);
      }

      /**
       * @brief Copy constructor
       */
      DerivativeBase(const DerivativeBase& other)
        : Parent(other),
          m_u(other.m_u)
      {}

      /**
       * @brief Move constructor
       */
      DerivativeBase(DerivativeBase&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u))
      {}

      constexpr
      size_t getDimension() const
      {
        return m_u.get().getFiniteElementSpace().getMesh().getSpaceDimension();
      }

      ScalarType getValue(const Geometry::Point& p) const
      {
        ScalarType out;
        getValue(out, p);
        return out;
      }

      void getValue(ScalarType& out, const Geometry::Point& p) const
      {
        const auto& polytope = p.getPolytope();
        const auto& polytopeMesh = polytope.getMesh();
        const auto& gf = getOperand();
        const auto& fes = gf.getFiniteElementSpace();
        const auto& fesMesh = fes.getMesh();
        if (polytopeMesh == fesMesh)
        {
          interpolate(out, p);
        }
        else if (const auto inclusion = fesMesh.inclusion(p))
        {
          interpolate(out, *inclusion);
        }
        else if (fesMesh.isSubMesh())
        {
          const auto& submesh = fesMesh.asSubMesh();
          const auto restriction = submesh.restriction(p);
          interpolate(out, *restriction);
        }
        else
        {
          assert(false);
        }
      }

      /**
       * @brief Interpolation function to be overriden in Derived type.
       */
      constexpr
      void interpolate(ScalarType& out, const Geometry::Point& p) const
      {
        static_cast<const Derived&>(*this).interpolate(out, p);
      }

      constexpr
      const OperandType& getOperand() const
      {
        return m_u.get();
      }

      /**
       * @brief Copy function to be overriden in Derived type.
       */
      DerivativeBase* copy() const noexcept override
      {
        return static_cast<const Derived&>(*this).copy();
      }

    private:
      std::reference_wrapper<const OperandType> m_u;
  };

  template <class NestedDerived, class FES, ShapeFunctionSpaceType SpaceType>
  class Derivative<ShapeFunction<NestedDerived, FES, SpaceType>> final
    : public ShapeFunctionBase<Derivative<ShapeFunction<NestedDerived, FES, SpaceType>>>
  {
    public:
      /// Finite element space type
      using FESType = FES;
      static constexpr ShapeFunctionSpaceType Space = SpaceType;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      /// Operand type
      using OperandType = ShapeFunction<NestedDerived, FESType, Space>;

      /// Parent class
      using Parent = ShapeFunctionBase<Derivative<OperandType>, FESType, Space>;

      Derivative(size_t i, const OperandType& u)
        : Parent(u.getFiniteElementSpace()),
          m_i(i),
          m_u(u)
      {}

      Derivative(const Derivative& other)
        : Parent(other),
          m_i(other.m_i),
          m_u(other.m_u)
      {}

      Derivative(Derivative&& other)
        : Parent(std::move(other)),
          m_i(other.m_i),
          m_u(std::move(other.m_u))
      {}

      constexpr
      const OperandType& getOperand() const
      {
        return m_u.get();
      }

      constexpr
      const auto& getLeaf() const
      {
        return getOperand().getLeaf();
      }

      constexpr
      RangeShape getRangeShape() const
      {
        return { 1, 1 };
      }

      constexpr
      size_t getDOFs(const Geometry::Polytope& element) const
      {
        return getOperand().getDOFs(element);
      }

      const Geometry::Point& getPoint() const
      {
        return m_p.value().get();
      }

      Derivative& setPoint(const Geometry::Point& p)
      {
        m_p = p;
        const auto& polytope = p.getPolytope();
        const size_t d = polytope.getDimension();
        const Index i = polytope.getIndex();
        const auto& fes = this->getFiniteElementSpace();
        const auto& fe = fes.getFiniteElement(d, i);
        const auto& rc = p.getReferenceCoordinates();
        const size_t dofs = this->getDOFs(p.getPolytope());
        m_gradients.resize(dofs);
        for (size_t local = 0; local < dofs; local++)
          m_gradients[local] = p.getJacobianInverse().transpose() * fe.getGradient(local)(rc);
        return *this;
      }

      auto getBasis(size_t local) const
      {
        return m_gradients[local].coeff(m_i);
      }

      Derivative* copy() const noexcept override
      {
        return new Derivative(*this);
      }

    private:
      size_t m_i;
      std::reference_wrapper<const OperandType> m_u;

      std::optional<std::reference_wrapper<const Geometry::Point>> m_p;

      std::vector<Math::SpatialVector<Real>> m_gradients;
  };

  template <class NestedDerived, class FES, ShapeFunctionSpaceType SpaceType>
  Derivative(size_t i, const ShapeFunction<NestedDerived, FES, SpaceType>& u)
    -> Derivative<ShapeFunction<NestedDerived, FES, SpaceType>>;

  /**
   * @brief %Utility function for computing @f$ \partial_x u @f$
   * @param[in] u GridFunction instance
   *
   * Given a scalar function @f$ u : \mathbb{R}^s \rightarrow \mathbb{R} @f$,
   * this function constructs the derivative in the @f$ x @f$ direction
   * @f$
   *   \dfrac{\partial u}{\partial x}
   * @f$
   */
  template <class Operand>
  auto Dx(const Operand& u)
  {
    return Derivative(0, u);
  }

  /**
   * @brief %Utility function for computing @f$ \partial_y u @f$
   * @param[in] u GridFunction instance
   *
   * Given a scalar function @f$ u : \mathbb{R}^s \rightarrow \mathbb{R} @f$,
   * this function constructs the derivative in the @f$ y @f$ direction
   * @f$
   *   \dfrac{\partial u}{\partial y}
   * @f$
   */
  template <class Operand>
  auto Dy(const Operand& u)
  {
    return Derivative(1, u);
  }

  /**
   * @brief %Utility function for computing @f$ \partial_z u @f$
   * @param[in] u GridFunction instance
   *
   * Given a scalar function @f$ u : \mathbb{R}^s \rightarrow \mathbb{R} @f$,
   * this function constructs the derivative in the @f$ y @f$ direction
   * @f$
   *   \dfrac{\partial u}{\partial z}
   * @f$
   */
  template <class Operand>
  auto Dz(const Operand& u)
  {
    return Derivative(2, u);
  }
}

#endif
