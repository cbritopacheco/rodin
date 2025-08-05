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

namespace Rodin::FormLanguage
{
  template <class Range, class Data, class Mesh>
  struct Traits<Variational::Grad<Variational::GridFunction<Variational::P1<Range, Mesh>, Data>>>
  {
    using FESType = Variational::P1<Range, Mesh>;

    using OperandType = Variational::GridFunction<FESType, Data>;

    using RangeType = Range;
  };

  template <class NestedDerived, class Range, class Mesh, Variational::ShapeFunctionSpaceType Space>
  struct Traits<
    Variational::Grad<
      Variational::ShapeFunction<NestedDerived, Variational::P1<Range, Mesh>, Space>>>
  {
    using FESType = Variational::P1<Range, Mesh>;
    static constexpr Variational::ShapeFunctionSpaceType SpaceType = Space;

    using OperandType = Variational::ShapeFunction<NestedDerived, FESType, SpaceType>;

    using RangeType = Range;
  };
}

namespace Rodin::Variational
{
  /**
   * @ingroup GradSpecializations
   * @brief Gradient of a P1 GridFunction
   */
  template <class Range, class Data, class Mesh>
  class Grad<GridFunction<P1<Range, Mesh>, Data>> final
    : public GradBase<GridFunction<P1<Range, Mesh>, Data>, Grad<GridFunction<P1<Range, Mesh>, Data>>>
  {
    public:
      using FESType = P1<Range, Mesh>;

      using RangeType = Range;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using SpatialVectorType = Math::SpatialVector<ScalarType>;

      using OperandType = GridFunction<FESType, Data>;

      using Parent = GradBase<OperandType, Grad<OperandType>>;

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

      void interpolate(SpatialVectorType& out, const Geometry::Point& p) const
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
          const auto& fe = fes.getFiniteElement(d, i);
          const auto& rc = p.getReferenceCoordinates();
          SpatialVectorType grad(d);
          SpatialVectorType res(d);
          res.setZero();
          for (size_t local = 0; local < fe.getCount(); local++)
          {
            const auto basis = fe.getBasis(local);
            for (size_t i = 0; i < d; i++)
              grad(i) = basis.template getDerivative<1>(i)(rc);
            res += gf[fes.getGlobalIndex({d, i}, local)] * grad;
          }
          out = p.getJacobianInverse().transpose() * res;
        }
      }

      Grad* copy() const noexcept override
      {
        return new Grad(*this);
      }
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for Grad of a P1 GridFunction
   */
  template <class Range, class Data, class Mesh>
  Grad(const GridFunction<P1<Range, Mesh>, Data>&) -> Grad<GridFunction<P1<Range, Mesh>, Data>>;

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

      using ScalarType = Number;

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
          m_u(std::move(other.m_u)),
          m_gradient(std::move(other.m_gradient)),
          m_p(std::exchange(other.m_p, nullptr))
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
      size_t getDOFs(const Geometry::Polytope& element) const
      {
        return getOperand().getDOFs(element);
      }

      const Geometry::Point& getPoint() const
      {
        assert(m_p);
        return *m_p;
      }

      Grad& setPoint(const Geometry::Point& p)
      {
        if (m_p == &p)
          return *this;
        m_p = &p;
        const auto& polytope = p.getPolytope();
        const auto& rc = p.getReferenceCoordinates();
        const size_t d = polytope.getDimension();
        const Index i = polytope.getIndex();
        const auto& fes = this->getFiniteElementSpace();
        const auto& fe = fes.getFiniteElement(d, i);
        const size_t count = fe.getCount();
        m_gradient.resize(count);
        for (size_t local = 0; local < count; local++)
        {
          m_gradient[local].resize(d);
          const auto basis = fe.getBasis(local);
          for (size_t j = 0; j < d; j++)
            m_gradient[local](j) = basis.template getDerivative<1>(j)(rc);
        }
        return *this;
      }

      constexpr
      auto getBasis(size_t local) const
      {
        return getPoint().getJacobianInverse().transpose() * m_gradient[local];
      }

      Grad* copy() const noexcept override
      {
        return new Grad(*this);
      }

    private:
      std::reference_wrapper<const OperandType> m_u;

      std::vector<Math::SpatialVector<ScalarType>> m_gradient;

      const Geometry::Point* m_p;
  };

  template <class NestedDerived, class Number, class Mesh, ShapeFunctionSpaceType Space>
  Grad(const ShapeFunction<NestedDerived, P1<Number, Mesh>, Space>&)
    -> Grad<ShapeFunction<NestedDerived, P1<Number, Mesh>, Space>>;
}

#endif

