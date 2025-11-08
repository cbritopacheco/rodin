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

namespace Rodin::FormLanguage
{
  template <class Range, class Data, class Mesh>
  struct Traits<
    Variational::Jacobian<
      Variational::GridFunction<
        Variational::P1<Range, Mesh>, Data>>>
  {
    using FESType = Variational::P1<Range, Mesh>;
    using OperandType = Variational::GridFunction<FESType, Data>;
  };

  template <class NestedDerived, class Range, class Mesh, Variational::ShapeFunctionSpaceType Space>
  struct Traits<
    Variational::Jacobian<
      Variational::ShapeFunction<NestedDerived, Variational::P1<Range, Mesh>, Space>>>
  {
    using FESType = Variational::P1<Range, Mesh>;
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
  template <class Range, class Data, class Mesh>
  class Jacobian<GridFunction<P1<Range, Mesh>, Data>> final
    : public JacobianBase<
        GridFunction<P1<Range, Mesh>, Data>, Jacobian<GridFunction<P1<Range, Mesh>, Data>>>
  {
    public:
      using FESType = P1<Range, Mesh>;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using SpatialMatrixType = Math::SpatialMatrix<ScalarType>;

      using OperandType = GridFunction<FESType, Data>;

      using Parent = JacobianBase<OperandType, Jacobian<OperandType>>;

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

      void interpolate(SpatialMatrixType& out, const Geometry::Point& p) const
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
            Math::SpatialPoint rc;
            tracePolytope->getTransformation().inverse(rc, pc);
            const Geometry::Point np(*tracePolytope, std::cref(rc), pc);
            this->interpolate(out, np);
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
                << ". Jacobian at an interface with no trace domain is undefined."
                << Alert::Raise;
            }
            else
            {
              for (auto& idx : inc)
              {
                const auto& tracePolytope = mesh.getPolytope(meshDim, idx);
                if (traceDomain.count(tracePolytope->getAttribute()))
                {
                  Math::SpatialPoint rc;
                  tracePolytope->getTransformation().inverse(rc, pc);
                  const Geometry::Point np(*tracePolytope, std::cref(rc), pc);
                  this->interpolate(out, np);
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
          SpatialMatrixType jacobian(d, d);
          SpatialMatrixType res(vdim, d);
          res.setZero();
          for (size_t local = 0; local < fe.getCount(); local++)
          {
            const auto& basis = fe.getBasis(local);
            for (size_t j = 0; j < d; j++)
            {
              for (size_t k = 0; k < d; k++)
                jacobian(j, k) = basis.template getDerivative<1>(j, k)(rc);
            }
            res += gf[fes.getGlobalIndex({d, i}, local)] * jacobian;
          }
          out = res * p.getJacobianInverse();
        }
      }

      Jacobian* copy() const noexcept override
      {
        return new Jacobian(*this);
      }
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for Jacobian of a P1 GridFunction
   */
  template <class Range, class Data, class Mesh>
  Jacobian(const GridFunction<P1<Range, Mesh>, Data>&) -> Jacobian<GridFunction<P1<Range, Mesh>, Data>>;

  /**
   * @ingroup JacobianSpecializations
   * @brief Jacobian of an P1 ShapeFunction object.
   */
  template <class ShapeFunctionDerived, class Range, class Mesh, ShapeFunctionSpaceType Space>
  class Jacobian<ShapeFunction<ShapeFunctionDerived, P1<Range, Mesh>, Space>> final
    : public ShapeFunctionBase<Jacobian<ShapeFunction<ShapeFunctionDerived, P1<Range, Mesh>, Space>>>
  {
    public:
      using FESType = P1<Range, Mesh>;
      static constexpr ShapeFunctionSpaceType SpaceType = Space;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using SpatialMatrixType = Math::SpatialMatrix<ScalarType>;

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

      constexpr
      const OperandType& getOperand() const
      {
        return m_u.get();
      }

      constexpr
      const FESType& getFiniteElementSpace() const
      {
        return getOperand().getFiniteElementSpace();
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
        return m_p.value().get();
      }

      Jacobian& setPoint(const Geometry::Point& p)
      {
        m_p = p;
        const auto& polytope = p.getPolytope();
        const auto& rc = p.getReferenceCoordinates();
        const size_t d = polytope.getDimension();
        const Index i = polytope.getIndex();
        const auto& fes = this->getFiniteElementSpace();
        const auto& fe = fes.getFiniteElement(d, i);
        const size_t count = fe.getCount();
        const size_t vdim = fes.getVectorDimension();
        m_jacobian.resize(count);
        for (size_t local = 0; local < count; local++)
        {
          m_jacobian[local].resize(vdim, d);
          for (size_t i = 0; i < vdim; i++)
          {
            const auto& basis = fe.getBasis(local);
            for (size_t j = 0; j < d; j++)
              m_jacobian[local](i, j) = basis.template getDerivative<1>(i, j)(rc);
          }
        }
        return *this;
      }

      constexpr
      auto getBasis(size_t local) const
      {
        return m_jacobian[local] * getPoint().getJacobianInverse();
      }

      Jacobian* copy() const noexcept override
      {
        return new Jacobian(*this);
      }

    private:
      std::reference_wrapper<const OperandType> m_u;

      std::vector<SpatialMatrixType> m_jacobian;

      Optional<std::reference_wrapper<const Geometry::Point>> m_p;
  };

  template <class ShapeFunctionDerived, class Number, class Mesh, ShapeFunctionSpaceType Space>
  Jacobian(const ShapeFunction<ShapeFunctionDerived, P1<Math::Vector<Number>, Mesh>, Space>&)
    -> Jacobian<ShapeFunction<ShapeFunctionDerived, P1<Math::Vector<Number>, Mesh>, Space>>;
}

#endif
