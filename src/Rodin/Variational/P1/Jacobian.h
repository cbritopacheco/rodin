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
#include "Rodin/Variational/Exceptions/UndeterminedTraceDomainException.h"

#include "P1Element.h"
#include "Rodin/Geometry/IsoparametricTransformation.h"

#include "GridFunction.h"

namespace Rodin::FormLanguage
{
  template <class Number, class Mesh>
  struct Traits<Variational::Jacobian<Variational::GridFunction<Variational::P1<Math::Vector<Number>, Mesh>>>>
  {
    using FESType = Variational::P1<Math::Vector<Number>>;
    using OperandType = Variational::GridFunction<FESType>;
  };

  template <class NestedDerived, class Number, class Mesh, Variational::ShapeFunctionSpaceType Space>
  struct Traits<
    Variational::Jacobian<
      Variational::ShapeFunction<NestedDerived, Variational::P1<Math::Vector<Number>, Mesh>, Space>>>
  {
    using FESType = Variational::P1<Math::Vector<Number>>;
    static constexpr Variational::ShapeFunctionSpaceType SpaceType = Space;
    using OperandType = Variational::ShapeFunction<NestedDerived, FESType, Space>;
  };
}

namespace Rodin::Variational
{
  /**
   * @ingroup JacobianSpecializations
   * @brief Jacobian of an P1 GridFunction object.
   */
  template <class Number, class Mesh>
  class Jacobian<GridFunction<P1<Math::Vector<Number>, Mesh>>> final
    : public JacobianBase<
        Jacobian<GridFunction<P1<Math::Vector<Number>, Mesh>>>, GridFunction<P1<Math::Vector<Number>, Mesh>>>
  {
    public:
      using FESType = P1<Math::Vector<Number>, Mesh>;

      using OperandType = GridFunction<FESType>;

      using Parent = JacobianBase<Jacobian<OperandType>, OperandType>;

      /**
       * @brief Constructs the Jacobian matrix of an @f$ H^1 (\Omega)^d @f$ function
       * @f$ u @f$.
       * @param[in] u Grid function to be differentiated
       */
      Jacobian(const OperandType& u)
        : Parent(u)
      {}

      Jacobian(const Jacobian& other)
        : Parent(other)
      {}

      Jacobian(Jacobian&& other)
        : Parent(std::move(other))
      {}

      void interpolate(Math::SpatialMatrix<Scalar>& out, const Geometry::Point& p) const
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
                << ". Grad at an interface with no trace domain is undefined."
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
          Math::SpatialMatrix<Scalar> res(vdim, d);
          res.setZero();
          for (size_t local = 0; local < fe.getCount(); local++)
          {
            fe.getJacobian(local)(jacobian, rc);
            res += gf.getValue(fes.getGlobalIndex({d, i}, local)).coeff(local % vdim) * jacobian;
          }
          out = res * p.getJacobianInverse();
        }
      }

      inline Jacobian* copy() const noexcept override
      {
        return new Jacobian(*this);
      }
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for Jacobian of a P1 GridFunction
   */
  template <class Number, class Mesh>
  Jacobian(const GridFunction<P1<Math::Vector<Number>, Mesh>>&)
    -> Jacobian<GridFunction<P1<Math::Vector<Number>, Mesh>>>;

  /**
   * @ingroup JacobianSpecializations
   * @brief Jacobian of an P1 ShapeFunction object.
   */
  template <class ShapeFunctionDerived, class Number, class Mesh, ShapeFunctionSpaceType Space>
  class Jacobian<ShapeFunction<ShapeFunctionDerived, P1<Math::Vector<Number>, Mesh>, Space>> final
    : public ShapeFunctionBase<Jacobian<ShapeFunction<ShapeFunctionDerived, P1<Math::Vector<Number>, Mesh>, Space>>>
  {
    public:
      using FESType = P1<Math::Vector<Number>, Mesh>;
      static constexpr ShapeFunctionSpaceType SpaceType = Space;

      using OperandType = ShapeFunction<ShapeFunctionDerived, FESType, SpaceType>;

      using Parent = ShapeFunctionBase<Jacobian<OperandType>, FESType, SpaceType>;

      Jacobian(const OperandType& u)
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
      const OperandType& getOperand() const
      {
        return m_u.get();
      }

      inline
      constexpr
      const FESType& getFiniteElementSpace() const
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
      const Geometry::Point& getPoint() const
      {
        return m_p.value().get();
      }

      Jacobian& setPoint(const Geometry::Point& p)
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
        const Math::Vector<Scalar>& rc = p.getCoordinates(Geometry::Point::Coordinates::Reference);
        return this->object(fe.getJacobian(local)(rc)) * p.getJacobianInverse();
      }

      inline Jacobian* copy() const noexcept override
      {
        return new Jacobian(*this);
      }

    private:
      std::reference_wrapper<const OperandType> m_u;

      std::optional<std::reference_wrapper<const Geometry::Point>> m_p;
  };

  template <class ShapeFunctionDerived, class Number, class Mesh, ShapeFunctionSpaceType Space>
  Jacobian(const ShapeFunction<ShapeFunctionDerived, P1<Math::Vector<Number>, Mesh>, Space>&)
    -> Jacobian<ShapeFunction<ShapeFunctionDerived, P1<Math::Vector<Number>, Mesh>, Space>>;
}

#endif
