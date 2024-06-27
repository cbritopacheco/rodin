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
  template <class NestedDerived, class ... Ps, Variational::ShapeFunctionSpaceType SpaceType>
  struct Traits<
    Variational::Div<
      Variational::ShapeFunction<NestedDerived, Variational::P1<Math::Vector<Scalar>, Ps...>, SpaceType>>>
  {
    using FES = Variational::P1<Math::Vector<Scalar>, Ps...>;
    static constexpr Variational::ShapeFunctionSpaceType Space = SpaceType;
  };
}

namespace Rodin::Variational
{
  /**
   * @ingroup DivSpecializations
   * @brief Divient of a P1 GridFunction
   */
  template <class Context, class Mesh>
  class Div<GridFunction<P1<Math::Vector<Scalar>, Context, Mesh>>> final
    : public DivBase<Div<GridFunction<P1<Math::Vector<Scalar>, Context, Mesh>>>, GridFunction<P1<Math::Vector<Scalar>, Context, Mesh>>>
  {
    public:
      using FES = P1<Math::Vector<Scalar>, Context, Mesh>;

      /// Operand type
      using Operand = GridFunction<FES>;

      /// Parent class
      using Parent = DivBase<Div<Operand>, Operand>;

      /**
       * @brief Constructs the Divient of an @f$ \mathbb{P}^1 @f$ function
       * @f$ u @f$.
       * @param[in] u P1 GridFunction
       */
      Div(const Operand& u)
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

      void interpolate(Math::Vector<Scalar>& out, const Geometry::Point& p) const
      {
        Math::SpatialVector<Scalar> tmp;
        interpolate(tmp, p);
        out = std::move(tmp);
      }

      void interpolate(Scalar& out, const Geometry::Point& p) const
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
            const Math::SpatialVector<Scalar> rc = tracePolytope->getTransformation().inverse(pc);
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
                  const Math::SpatialVector<Scalar> rc = tracePolytope->getTransformation().inverse(pc);
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
          Math::SpatialMatrix<Scalar> jacobian(vdim, d);
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
  template <class NestedDerived, ShapeFunctionSpaceType SpaceType, class ... Ts>
  class Div<ShapeFunction<NestedDerived, P1<Math::Vector<Scalar>, Ts...>, SpaceType>> final
    : public ShapeFunctionBase<Div<ShapeFunction<NestedDerived, P1<Math::Vector<Scalar>, Ts...>, SpaceType>>>
  {
    public:
      using FES = P1<Math::Vector<Scalar>, Ts...>;
      static constexpr ShapeFunctionSpaceType Space = SpaceType;

      using Operand = ShapeFunction<NestedDerived, FES, Space>;
      using Parent = ShapeFunctionBase<Div<ShapeFunction<NestedDerived, FES, Space>>, FES, Space>;

      /**
       * @brief Constructs Div object
       * @param[in] u ShapeFunction to be differentiated
       */
      Div(const Operand& u)
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
      const Operand& getOperand() const
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
      constexpr
      auto getTensorBasis(const Geometry::Point& p) const
      {
        const size_t d = p.getPolytope().getDimension();
        const Index i = p.getPolytope().getIndex();
        const auto& fe = this->getFiniteElementSpace().getFiniteElement(d, i);
        const Math::Vector<Scalar>& rc = p.getCoordinates(Geometry::Point::Coordinates::Reference);
        return TensorBasis(fe.getCount(),
            [&](size_t local) -> Scalar
            { return (fe.getJacobian(local)(rc) * p.getJacobianInverse()).trace(); });
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
      std::reference_wrapper<const Operand> m_u;
  };

  template <class NestedDerived, ShapeFunctionSpaceType Space, class ... Ts>
  Div(const ShapeFunction<NestedDerived, P1<Math::Vector<Scalar>, Ts...>, Space>&)
    -> Div<ShapeFunction<NestedDerived, P1<Math::Vector<Scalar>, Ts...>, Space>>;
}

#endif
