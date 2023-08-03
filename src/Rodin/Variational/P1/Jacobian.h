/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P1_JACOBIAN_H
#define RODIN_VARIATIONAL_P1_JACOBIAN_H

#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Variational/Jacobian.h"

#include "P1Element.h"
#include "Rodin/Geometry/IsoparametricTransformation.h"

#include "GridFunction.h"

namespace Rodin::Variational
{
  /**
   * @ingroup JacobianSpecializations
   * @brief Jacobian of an P1 GridFunction object.
   */
  template <class ... Ps>
  class Jacobian<GridFunction<P1<Math::Vector, Ps...>>> final
    : public MatrixFunctionBase<Jacobian<GridFunction<P1<Math::Vector, Ps...>>>>
  {
    public:
      using Operand = GridFunction<P1<Math::Vector, Ps...>>;
      using Parent = MatrixFunctionBase<Jacobian<GridFunction<P1<Math::Vector, Ps...>>>>;

      /**
       * @brief Constructs the Jacobian matrix of an @f$ H^1 (\Omega)^d @f$ function
       * @f$ u @f$.
       * @param[in] u Grid function to be differentiated
       */
      Jacobian(const Operand& u)
        :  m_u(u)
      {}

      Jacobian(const Jacobian& other)
        : Parent(other),
          m_u(other.m_u)
      {}

      Jacobian(Jacobian&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u))
      {}

      inline
      constexpr
      size_t getRows() const
      {
        return m_u.get().getFiniteElementSpace().getVectorDimension();
      }

      inline
      constexpr
      size_t getColumns() const
      {
        return m_u.get().getFiniteElementSpace().getMesh().getSpaceDimension();
      }

      Math::Matrix getValue(const Geometry::Point& p) const
      {
        const auto& polytope = p.getPolytope();
        const auto& d = polytope.getDimension();
        const auto& i = polytope.getIndex();
        const auto& gf = m_u.get();
        const auto& fes = gf.getFiniteElementSpace();
        const auto& vdim = fes.getVectorDimension();
        const auto& fe = fes.getFiniteElement(d, i);
        const auto& rc = p.getReferenceCoordinates();
        Math::Matrix jacobian(vdim, d);
        Math::Matrix res(vdim, d);
        res.setZero();
        for (size_t local = 0; local < fe.getCount(); local++)
        {
          fe.getJacobian(local)(jacobian, rc);
          res += gf.getValue(fes.getGlobalIndex({d, i}, local)).coeff(local % vdim) * jacobian;
        }
        return this->object(std::move(res)) * p.getJacobianInverse();
      }

      inline
      constexpr
      const Operand& getOperand() const
      {
        return m_u.get();
      }

      inline
      constexpr
      Jacobian& traceOf(Geometry::Attribute attr)
      {
        Parent::traceOf(attr);
        return *this;
      }

      inline
      constexpr
      Jacobian& traceOf(const FlatSet<Geometry::Attribute>& attrs)
      {
        Parent::traceOf(attrs);
        return *this;
      }

      inline Jacobian* copy() const noexcept override
      {
        return new Jacobian(*this);
      }

    private:
      std::reference_wrapper<const Operand> m_u;
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for Jacobian of a P1 GridFunction
   */
  template <class ... Ps>
  Jacobian(const GridFunction<P1<Math::Vector, Ps...>>&)
    -> Jacobian<GridFunction<P1<Math::Vector, Ps...>>>;

  /**
   * @ingroup JacobianSpecializations
   * @brief Jacobian of an P1 ShapeFunction object.
   */
  template <class ShapeFunctionDerived, ShapeFunctionSpaceType SpaceType, class ... Ts>
  class Jacobian<ShapeFunction<ShapeFunctionDerived, P1<Math::Vector, Ts...>, SpaceType>> final
    : public ShapeFunctionBase<Jacobian<ShapeFunction<ShapeFunctionDerived, P1<Math::Vector, Ts...>, SpaceType>>, P1<Math::Vector, Ts...>, SpaceType>
  {
    public:
      using FES = P1<Math::Vector, Ts...>;
      static constexpr ShapeFunctionSpaceType Space = SpaceType;

      using Operand = ShapeFunction<ShapeFunctionDerived, FES, Space>;
      using Parent = ShapeFunctionBase<Jacobian<Operand>, FES, Space>;

      Jacobian(const Operand& u)
        : Parent(u.getFiniteElementSpace()),
          m_u(u)
      {}

      Jacobian(const Jacobian& other)
        : Parent(other),
          m_u(other.m_u)
      {}

      Jacobian(Jacobian&& other)
        : Parent(std::move(other)),
          m_u(other.m_u)
      {}

      inline
      constexpr
      const Operand& getOperand() const
      {
        return m_u.get();
      }

      inline
      constexpr
      const FES& getFiniteElementSpace() const
      {
        return getOperand().getFiniteElementSpace();
      }

      inline
      constexpr
      const auto& getLeaf() const
      {
        return getOperand().getLeaf();
      }

      inline
      constexpr
      RangeShape getRangeShape() const
      {
        return { getOperand().getFiniteElementSpace().getMesh().getSpaceDimension(),
                 getOperand().getFiniteElementSpace().getVectorDimension() };
      }

      inline
      constexpr
      size_t getDOFs(const Geometry::Polytope& element) const
      {
        return getOperand().getDOFs(element);
      }

      inline
      constexpr
      auto getTensorBasis(const Geometry::Point& p) const
      {
        const size_t d = p.getPolytope().getDimension();
        const Index i = p.getPolytope().getIndex();
        const auto& fe = this->getFiniteElementSpace().getFiniteElement(d, i);
        const Math::Vector& rc = p.getCoordinates(Geometry::Point::Coordinates::Reference);
        return TensorBasis(fe.getCount(),
            [&](size_t local)
            { return this->object(fe.getJacobian(local)(rc)) * p.getJacobianInverse(); });
      }

      inline Jacobian* copy() const noexcept override
      {
        return new Jacobian(*this);
      }

    private:
      std::reference_wrapper<const Operand> m_u;
  };

  template <class ShapeFunctionDerived, class FES, ShapeFunctionSpaceType Space>
  Jacobian(const ShapeFunction<ShapeFunctionDerived, FES, Space>&)
    -> Jacobian<ShapeFunction<ShapeFunctionDerived, FES, Space>>;
}

#endif
