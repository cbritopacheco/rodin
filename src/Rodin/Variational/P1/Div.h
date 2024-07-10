/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P1_DIV_H
#define RODIN_VARIATIONAL_P1_DIV_H

#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Variational/Div.h"

#include "GridFunction.h"

namespace Rodin::FormLanguage
{
  template <class Scalar, class Mesh>
  struct Traits<Variational::Div<Variational::GridFunction<Variational::P1<Math::Vector<Scalar>, Mesh>>>>
  {
    using FESType = Variational::P1<Math::Vector<Scalar>, Mesh>;
    using ScalarType = Scalar;
    using OperandType = Variational::GridFunction<Variational::P1<Math::Vector<Scalar>, Mesh>>;
  };

  template <class NestedDerived, class Scalar, class Mesh, Variational::ShapeFunctionSpaceType Space>
  struct Traits<
    Variational::Div<
      Variational::ShapeFunction<NestedDerived, Variational::P1<Math::Vector<Scalar>, Mesh>, Space>>>
  {
    using FESType = Variational::P1<Math::Vector<Scalar>, Mesh>;
    static constexpr Variational::ShapeFunctionSpaceType SpaceType = Space;
    using ScalarType = Scalar;
    using OperandType =
      Variational::ShapeFunction<NestedDerived, Variational::P1<Math::Vector<Scalar>, Mesh>, Space>;
  };
}

namespace Rodin::Variational
{
  /**
   * @ingroup DivSpecializations
   * @brief Divient of a P1 GridFunction
   */
  template <class Scalar, class Mesh>
  class Div<GridFunction<P1<Math::Vector<Scalar>, Mesh>>> final
    : public DivBase<
        GridFunction<P1<Math::Vector<Scalar>, Mesh>>,
        Div<GridFunction<P1<Math::Vector<Scalar>, Mesh>>>>
  {
    public:
      using FESType = Variational::P1<Math::Vector<Scalar>, Mesh>;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using SpatialMatrixType = Math::SpatialMatrix<ScalarType>;

      /// Operand type
      using OperandType = GridFunction<FESType>;

      /// Parent class
      using Parent = DivBase<OperandType, Div<OperandType>>;

      /**
       * @brief Constructs the Divient of an @f$ \mathbb{P}^1 @f$ function
       * @f$ u @f$.
       * @param[in] u P1 GridFunction
       */
      Div(const OperandType& u)
        : Parent(u)
      {}

      /**
       * @brief Copy constructor
       */
      Div(const Div& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor
       */
      Div(Div&& other)
        : Parent(std::move(other))
      {}

      void interpolate(ScalarType& out, const Geometry::Point& p) const
      {
        const auto& polytope = p.getPolytope();
        const auto& d = polytope.getDimension();
        const auto& i = polytope.getIndex();
        const auto& mesh = polytope.getMesh();
        const size_t meshDim = mesh.getDimension();
        if (d == meshDim - 1) // Evaluating on a face
        {
          const auto& conn = mesh.getConnectivity();
          const auto& inc = conn.getIncidence({ meshDim - 1, meshDim }, i);
          const auto& pc = p.getPhysicalCoordinates();
          assert(inc.size() == 1 || inc.size() == 2);
          if (inc.size() == 1)
          {
            const auto& tracePolytope = mesh.getPolytope(meshDim, *inc.begin());
            const auto rc = tracePolytope->getTransformation().inverse(pc);
            const Geometry::Point np(*tracePolytope, std::cref(rc), pc);
            interpolate(out, np);
            return;
          }
          else
          {
            assert(inc.size() == 2);
            const auto& traceDomain = this->getTraceDomain();
            assert(traceDomain.size() > 0);
            if (traceDomain.size() == 0)
            {
              Alert::MemberFunctionException(*this, __func__)
                << "No trace domain provided: "
                << Alert::Notation::Predicate(true, "getTraceDomain().size() == 0")
                << ". Div at an interface with no trace domain is undefined."
                << Alert::Raise;
            }
            else
            {
              for (auto& idx : inc)
              {
                const auto& tracePolytope = mesh.getPolytope(meshDim, idx);
                if (traceDomain.count(tracePolytope->getAttribute()))
                {
                  const auto rc = tracePolytope->getTransformation().inverse(pc);
                  const Geometry::Point np(*tracePolytope, std::cref(rc), pc);
                  interpolate(out, np);
                  return;
                }
              }
              UndeterminedTraceDomainException(
                  *this, __func__, {d, i}, traceDomain.begin(), traceDomain.end()) << Alert::Raise;
            }
            return;
          }
        }
        else // Evaluating on a cell
        {
          assert(d == mesh.getDimension());
          const auto& gf = this->getOperand();
          const auto& fes = gf.getFiniteElementSpace();
          const auto& vdim = fes.getVectorDimension();
          const auto& fe = fes.getFiniteElement(d, i);
          const auto& rc = p.getReferenceCoordinates();
          SpatialMatrixType jacobian(vdim, d);
          out = 0;
          for (size_t local = 0; local < fe.getCount(); local++)
          {
            fe.getJacobian(local)(jacobian, rc);
            out += gf.getValue(fes.getGlobalIndex({d, i}, local)).coeff(local % vdim) * (jacobian * p.getJacobianInverse()).trace();
          }
        }
      }

      inline Div* copy() const noexcept override
      {
        return new Div(*this);
      }
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for Div of a P1 GridFunction
   */
  template <class ... Ts>
  Div(const GridFunction<P1<Ts...>>&) -> Div<GridFunction<P1<Ts...>>>;

  /**
   * @ingroup DivSpecializations
   */
  template <class NestedDerived, class Number, class Mesh, ShapeFunctionSpaceType Space>
  class Div<ShapeFunction<NestedDerived, P1<Math::Vector<Number>, Mesh>, Space>> final
    : public ShapeFunctionBase<Div<ShapeFunction<NestedDerived, P1<Math::Vector<Number>, Mesh>, Space>>>
  {
    public:
      using FESType = P1<Math::Vector<Number>, Mesh>;
      static constexpr ShapeFunctionSpaceType SpaceType = Space;

      using OperandType = ShapeFunction<NestedDerived, FESType, SpaceType>;

      using Parent = ShapeFunctionBase<Div<OperandType>, FESType, SpaceType>;

      /**
       * @brief Constructs Div object
       * @param[in] u ShapeFunction to be differentiated
       */
      Div(const OperandType& u)
        : Parent(u.getFiniteElementSpace()),
          m_u(u)
      {}

      Div(const Div& other)
        : Parent(other),
          m_u(other.m_u)
      {}

      Div(Div&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u))
      {}

      inline
      constexpr
      const OperandType& getOperand() const
      {
        return m_u.get();
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
        return { 1, 1 };
      }

      inline
      constexpr
      size_t getDOFs(const Geometry::Polytope& polytope) const
      {
        return getOperand().getDOFs(polytope);
      }

      inline
      const Geometry::Point& getPoint() const
      {
        return m_p.value().get();
      }

      Div& setPoint(const Geometry::Point& p)
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
        const auto& fe = this->getFiniteElementSpace().getFiniteElement(d, i);
        const auto& rc = p.getReferenceCoordinates();
        return (fe.getJacobian(local)(rc) * p.getJacobianInverse()).trace();
      }

      inline
      constexpr
      const auto& getFiniteElementSpace() const
      {
        return getOperand().getFiniteElementSpace();
      }

      inline Div* copy() const noexcept override
      {
        return new Div(*this);
      }

    private:
      std::reference_wrapper<const OperandType> m_u;

      std::optional<std::reference_wrapper<const Geometry::Point>> m_p;
  };

  template <class NestedDerived, class Number, class Mesh, ShapeFunctionSpaceType Space>
  Div(const ShapeFunction<NestedDerived, P1<Math::Vector<Number>, Mesh>, Space>&)
    -> Div<ShapeFunction<NestedDerived, P1<Math::Vector<Number>, Mesh>, Space>>;
}

#endif
