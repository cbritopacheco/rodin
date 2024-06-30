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
  template <class Range, class Mesh>
  struct Traits<Variational::Grad<Variational::GridFunction<Variational::P1<Range, Mesh>>>>
  {
    using FESType = Variational::P1<Range, Mesh>;
    using OperandType = Variational::GridFunction<FESType>;
  };

  template <class NestedDerived, class Range, class Mesh, Variational::ShapeFunctionSpaceType Space>
  struct Traits<
    Variational::Grad<
      Variational::ShapeFunction<NestedDerived, Variational::P1<Range, Mesh>, Space>>>
  {
    using FESType = Variational::P1<Range, Mesh>;
    static constexpr Variational::ShapeFunctionSpaceType SpaceType = Space;
    using OperandType = Variational::ShapeFunction<NestedDerived, FESType, SpaceType>;
  };
}

namespace Rodin::Variational
{
  /**
   * @ingroup GradSpecializations
   * @brief Gradient of a P1 GridFunction
   */
  template <class Mesh>
  class Grad<GridFunction<P1<Scalar, Mesh>>> final
    : public GradBase<Grad<GridFunction<P1<Scalar, Mesh>>>, GridFunction<P1<Scalar, Mesh>>>
  {
    public:
      using FESType = P1<Scalar, Mesh>;

      /// Operand type
      using OperandType = GridFunction<FESType>;

      /// Parent class
      using Parent = GradBase<Grad<OperandType>, OperandType>;

      /**
       * @brief Constructs the gradient of an @f$ \mathbb{P}^1 @f$ function
       * @f$ u @f$.
       * @param[in] u P1 GridFunction
       */
      Grad(const OperandType& u)
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

      void interpolate(Math::Vector<Scalar>& out, const Geometry::Point& p) const
      {
        Math::SpatialVector<Scalar> tmp;
        interpolate(tmp, p);
        out = std::move(tmp);
      }

      void interpolate(Math::SpatialVector<Scalar>& out, const Geometry::Point& p) const
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
          const auto& fe = fes.getFiniteElement(d, i);
          const auto& rc = p.getReferenceCoordinates();
          Math::SpatialVector<Scalar> grad(d);
          Math::SpatialVector<Scalar> res(d);
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
  template <class Number, class Mesh>
  Grad(const GridFunction<P1<Number, Mesh>>&) -> Grad<GridFunction<P1<Number, Mesh>>>;

  /**
   * @ingroup GradSpecializations
   * @brief Gradient of a P1 ShapeFunction
   */
  template <class NestedDerived, class Number, class Mesh, ShapeFunctionSpaceType SpaceType>
  class Grad<ShapeFunction<NestedDerived, P1<Number, Mesh>, SpaceType>> final
    : public ShapeFunctionBase<Grad<ShapeFunction<NestedDerived, P1<Number, Mesh>, SpaceType>>>
  {
    public:
      /// Finite element space type
      using FESType = P1<Number, Mesh>;
      static constexpr ShapeFunctionSpaceType Space = SpaceType;

      using NumberType = Number;

      /// Operand type
      using OperandType = ShapeFunction<NestedDerived, FESType, Space>;

      /// Parent class
      using Parent = ShapeFunctionBase<Grad<OperandType>, FESType, Space>;

      Grad(const OperandType& u)
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
        s_gradient.resize(dofs);
        for (size_t local = 0; local < dofs; local++)
          s_gradient[local] = fe.getGradient(local)(rc);
        return TensorBasis(dofs,
            [&](size_t local) { return p.getJacobianInverse().transpose() * s_gradient[local]; });
      }

      inline Grad* copy() const noexcept override
      {
        return new Grad(*this);
      }

    private:
      thread_local static std::vector<Math::SpatialVector<NumberType>> s_gradient;
      std::reference_wrapper<const OperandType> m_u;
  };

  template <class NestedDerived, class Number, class Mesh, ShapeFunctionSpaceType SpaceType>
  thread_local
  std::vector<Math::SpatialVector<Number>>
  Grad<ShapeFunction<NestedDerived, P1<Number, Mesh>, SpaceType>>::s_gradient;

  template <class NestedDerived, class Number, class Mesh, ShapeFunctionSpaceType Space>
  Grad(const ShapeFunction<NestedDerived, P1<Number, Mesh>, Space>&)
    -> Grad<ShapeFunction<NestedDerived, P1<Number, Mesh>, Space>>;
}

#endif

