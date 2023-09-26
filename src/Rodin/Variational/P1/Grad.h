/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P1_GRADIENT_H
#define RODIN_VARIATIONAL_P1_GRADIENT_H

#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Variational/Grad.h"

#include "Rodin/Variational/Exceptions/UndeterminedTraceDomainException.h"

#include "GridFunction.h"

namespace Rodin::FormLanguage
{
  template <class ... Ps>
  struct Traits<
    Variational::Grad<Variational::GridFunction<Variational::P1<Scalar, Ps...>>>>
  {
    using FES = Variational::P1<Scalar, Ps...>;
    using Operand = Variational::GridFunction<FES>;
  };

  template <class NestedDerived, class ... Ps, Variational::ShapeFunctionSpaceType SpaceType>
  struct Traits<
    Variational::Grad<
      Variational::ShapeFunction<NestedDerived, Variational::P1<Scalar, Ps...>, SpaceType>>>
  {
    using FES = Variational::P1<Scalar, Ps...>;
    static constexpr Variational::ShapeFunctionSpaceType Space = SpaceType;
  };
}

namespace Rodin::Variational
{
  /**
   * @ingroup GradSpecializations
   * @brief Gradient of a P1 GridFunction
   */
  template <class Context, class Mesh>
  class Grad<GridFunction<P1<Scalar, Context, Mesh>>> final
    : public GradBase<Grad<GridFunction<P1<Scalar, Context, Mesh>>>, GridFunction<P1<Scalar, Context, Mesh>>>
  {
    public:
      using FES = P1<Scalar, Context, Mesh>;

      /// Operand type
      using Operand = GridFunction<FES>;

      /// Parent class
      using Parent = GradBase<Grad<Operand>, Operand>;

      /**
       * @brief Constructs the gradient of an @f$ \mathbb{P}^1 @f$ function
       * @f$ u @f$.
       * @param[in] u P1 GridFunction
       */
      Grad(const Operand& u)
        : Parent(u)
      {}

      /**
       * @brief Copy constructor
       */
      Grad(const Grad& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor
       */
      Grad(Grad&& other)
        : Parent(std::move(other))
      {}

      void interpolate(Math::Vector& out, const Geometry::Point& p) const
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
          const auto& rc = p.getPhysicalCoordinates();
          const auto& pc = p.getPhysicalCoordinates();
          assert(inc.size() == 1 || inc.size() == 2);
          if (inc.size() == 1)
          {
            const auto& tracePolytope = mesh.getPolytope(meshDim, *inc.begin());
            const Geometry::Point np(*tracePolytope, rc, pc);
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
                  const Geometry::Point np(*tracePolytope, rc, pc);
                  interpolate(out, np);
                  return;
                }
              }

              UndeterminedTraceDomainException(
                  *this, __func__, {d, i}, traceDomain.begin(), traceDomain.end()).raise();
            }
            return;
          }
        }
        else // Evaluating on a cell
        {
          assert(d == mesh.getDimension());
          const auto& gf = this->getOperand();
          const auto& fes = gf.getFiniteElementSpace();
          const auto& fe = fes.getFiniteElement(d, i);
          const auto& rc = p.getReferenceCoordinates();
          Math::SpatialVector grad(d);
          Math::SpatialVector res(d);
          res.setZero();
          for (size_t local = 0; local < fe.getCount(); local++)
          {
            fe.getGradient(local)(grad, rc);
            res += gf.getValue(fes.getGlobalIndex({d, i}, local)) * grad;
          }
          out = p.getJacobianInverse().transpose() * res;
        }
      }

      inline Grad* copy() const noexcept override
      {
        return new Grad(*this);
      }
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for Grad of a P1 GridFunction
   */
  template <class ... Ts>
  Grad(const GridFunction<P1<Ts...>>&) -> Grad<GridFunction<P1<Ts...>>>;

  /**
   * @ingroup GradSpecializations
   * @brief Gradient of a P1 ShapeFunction
   */
  template <class NestedDerived, class ... Ps, ShapeFunctionSpaceType SpaceType>
  class Grad<ShapeFunction<NestedDerived, P1<Scalar, Ps...>, SpaceType>> final
    : public ShapeFunctionBase<Grad<ShapeFunction<NestedDerived, P1<Scalar, Ps...>, SpaceType>>>
  {
    public:
      /// Finite element space type
      using FES = P1<Scalar, Ps...>;
      static constexpr ShapeFunctionSpaceType Space = SpaceType;

      /// Operand type
      using Operand = ShapeFunction<NestedDerived, P1<Scalar, Ps...>, Space>;

      /// Parent class
      using Parent = ShapeFunctionBase<Grad<Operand>, FES, Space>;

      Grad(const Operand& u)
        : Parent(u.getFiniteElementSpace()),
          m_u(u)
      {}

      Grad(const Grad& other)
        : Parent(other),
          m_u(other.m_u)
      {}

      Grad(Grad&& other)
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
        return { getOperand().getFiniteElementSpace().getMesh().getSpaceDimension(), 1 };
      }

      inline
      constexpr
      size_t getDOFs(const Geometry::Polytope& element) const
      {
        return getOperand().getDOFs(element);
      }

      inline
      auto getTensorBasis(const Geometry::Point& p) const
      {
        const size_t d = p.getPolytope().getDimension();
        const auto& fes = this->getFiniteElementSpace();
        const Index i = p.getPolytope().getIndex();
        const auto& fe = fes.getFiniteElement(d, i);
        const size_t dofs = fe.getCount();
        const auto& rc = p.getCoordinates(Geometry::Point::Coordinates::Reference);
        m_gradient.resize(dofs);
        for (size_t local = 0; local < dofs; local++)
          m_gradient[local] = fe.getGradient(local)(rc);
        return TensorBasis(dofs,
            [&](size_t local) { return p.getJacobianInverse().transpose() * m_gradient[local]; });
      }

      inline Grad* copy() const noexcept override
      {
        return new Grad(*this);
      }

    private:
      std::reference_wrapper<const Operand> m_u;
      mutable std::vector<Math::SpatialVector> m_gradient;
  };

  template <class NestedDerived, class ... Ps, ShapeFunctionSpaceType Space>
  Grad(const ShapeFunction<NestedDerived, P1<Scalar, Ps...>, Space>&)
    -> Grad<ShapeFunction<NestedDerived, P1<Scalar, Ps...>, Space>>;
}

#endif

