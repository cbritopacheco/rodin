/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_H1_DIV_H
#define RODIN_VARIATIONAL_H1_DIV_H

/**
 * @file
 * @brief Divergence operator specialization for H1 vector functions.
 *
 * This file provides specialized implementations of the divergence operator
 * for H1 vector-valued GridFunctions and ShapeFunctions.
 *
 * For H1<K> elements, the divergence is computed as:
 * @f[
 *   \nabla \cdot \mathbf{u}|_K = \sum_{i=1}^{n_{dof}} \sum_{j=1}^d u_{i,j} \frac{\partial \phi_i}{\partial x_j}
 * @f]
 * where @f$ \phi_i @f$ are H1<K> basis functions and @f$ u_{i,j} @f$ are DOF values.
 *
 * For H1<K> functions, the divergence is polynomial of degree K-1 on each element.
 */

#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Variational/Div.h"

namespace Rodin::FormLanguage
{
  template <size_t K, class Scalar, class Data, class Mesh>
  struct Traits<Variational::Div<Variational::GridFunction<Variational::H1<K, Math::Vector<Scalar>, Mesh>, Data>>>
  {
    using FESType = Variational::H1<K, Math::Vector<Scalar>, Mesh>;
    using ScalarType = Scalar;
    using OperandType = Variational::GridFunction<Variational::H1<K, Math::Vector<Scalar>, Mesh>, Data>;
  };

  template <size_t K, class NestedDerived, class Scalar, class Mesh, Variational::ShapeFunctionSpaceType Space>
  struct Traits<
    Variational::Div<
      Variational::ShapeFunction<NestedDerived, Variational::H1<K, Math::Vector<Scalar>, Mesh>, Space>>>
  {
    using FESType = Variational::H1<K, Math::Vector<Scalar>, Mesh>;
    static constexpr Variational::ShapeFunctionSpaceType SpaceType = Space;
    using ScalarType = Scalar;
    using OperandType =
      Variational::ShapeFunction<NestedDerived, Variational::H1<K, Math::Vector<Scalar>, Mesh>, Space>;
  };
}

namespace Rodin::Variational
{
  /**
   * @ingroup DivSpecializations
   * @brief Divergence of an H1 vector GridFunction.
   *
   * Computes @f$ \nabla \cdot \mathbf{u} @f$ for vector-valued H1<K> functions.
   * The divergence is polynomial of degree K-1 on each element:
   * @f[
   *   (\nabla \cdot \mathbf{u})|_K = \sum_{i=1}^{n_{dof}} \mathbf{u}_i \cdot \nabla \phi_i
   * @f]
   *
   * Applications include:
   * - Incompressibility constraint: @f$ \nabla \cdot \mathbf{v} = 0 @f$
   * - Mass conservation in fluid mechanics
   * - Continuity equation
   * - Mixed formulations (Stokes, elasticity)
   */
  template <size_t K, class Scalar, class Data, class Mesh>
  class Div<GridFunction<H1<K, Math::Vector<Scalar>, Mesh>, Data>> final
    : public DivBase<
        GridFunction<H1<K, Math::Vector<Scalar>, Mesh>, Data>,
        Div<GridFunction<H1<K, Math::Vector<Scalar>, Mesh>, Data>>>
  {
    public:
      using FESType = Variational::H1<K, Math::Vector<Scalar>, Mesh>;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using SpatialMatrixType = Math::SpatialMatrix<ScalarType>;

      /// Operand type
      using OperandType = GridFunction<FESType, Data>;

      /// Parent class
      using Parent = DivBase<OperandType, Div<OperandType>>;

      /**
       * @brief Constructs the divergence of an H1 vector function.
       * @param[in] u H1 vector GridFunction to differentiate
       *
       * Creates divergence operator @f$ \nabla \cdot \mathbf{u} @f$ where
       * @f$ \mathbf{u} \in [H^1]^d @f$.
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
            Math::SpatialPoint rc;
            tracePolytope->getTransformation().inverse(rc, pc);
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
                  Math::SpatialPoint rc;
                  tracePolytope->getTransformation().inverse(rc, pc);
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
          SpatialMatrixType jacobian(d, d);
          out = 0;
          for (size_t local = 0; local < fe.getCount(); local++)
          {
            for (size_t j = 0; j < d; j++)
            {
              const auto& basis = fe.getBasis(local);
              for (size_t k = 0; k < d; k++)
                jacobian(j, k) = basis.template getDerivative<1>(j, k)(rc);
            }
            out += gf[fes.getGlobalIndex({d, i}, local)] * (jacobian * p.getJacobianInverse()).trace();
          }
        }
      }

      Div* copy() const noexcept override
      {
        return new Div(*this);
      }
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for Div of an H1 GridFunction
   */
  template <size_t K, class Scalar, class Data, class Mesh>
  Div(const GridFunction<H1<K, Math::Vector<Scalar>, Mesh>, Data>&)
    -> Div<GridFunction<H1<K, Math::Vector<Scalar>, Mesh>, Data>>;

  /**
   * @ingroup DivSpecializations
   * @brief Divergence of an H1 vector ShapeFunction.
   *
   * Represents @f$ \nabla \cdot \mathbf{v} @f$ for H1 test/trial functions.
   * Used in weak formulations like:
   * @f[
   *   \int_\Omega p (\nabla \cdot \mathbf{v}) \, dx
   * @f]
   * for pressure-velocity coupling in Stokes/Navier-Stokes equations.
   */
  template <size_t K, class NestedDerived, class Number, class Mesh, ShapeFunctionSpaceType Space>
  class Div<ShapeFunction<NestedDerived, H1<K, Math::Vector<Number>, Mesh>, Space>> final
    : public ShapeFunctionBase<Div<ShapeFunction<NestedDerived, H1<K, Math::Vector<Number>, Mesh>, Space>>>
  {
    public:
      using FESType = H1<K, Math::Vector<Number>, Mesh>;
      static constexpr ShapeFunctionSpaceType SpaceType = Space;

      using OperandType = ShapeFunction<NestedDerived, FESType, SpaceType>;

      using Parent = ShapeFunctionBase<Div<OperandType>, FESType, SpaceType>;

      using ScalarType = Number;

      /**
       * @brief Constructs divergence of H1 ShapeFunction.
       * @param[in] u ShapeFunction to differentiate
       *
       * Creates @f$ \nabla \cdot \mathbf{u} @f$ for use in bilinear forms.
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
          m_u(std::move(other.m_u)),
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
      size_t getDOFs(const Geometry::Polytope& polytope) const
      {
        return getOperand().getDOFs(polytope);
      }

      const Geometry::Point& getPoint() const
      {
        assert(m_p);
        return *m_p;
      }

      Div& setPoint(const Geometry::Point& p)
      {
        m_p = &p;
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
          const auto& basis = fe.getBasis(local);
          for (size_t i = 0; i < vdim; i++)
          {
            for (size_t j = 0; j < d; j++)
              m_jacobian[local](i, j) = basis.template getDerivative<1>(i, j)(rc);
          }
        }
        return *this;
      }

      constexpr
      ScalarType getBasis(size_t local) const
      {
        const auto& p = this->getPoint();
        return (m_jacobian[local] * p.getJacobianInverse()).trace();
      }

      constexpr
      const auto& getFiniteElementSpace() const
      {
        return getOperand().getFiniteElementSpace();
      }

      Div* copy() const noexcept override
      {
        return new Div(*this);
      }

    private:
      std::reference_wrapper<const OperandType> m_u;

      std::vector<Math::SpatialMatrix<ScalarType>> m_jacobian;

      const Geometry::Point* m_p;
  };

  template <size_t K, class NestedDerived, class Number, class Mesh, ShapeFunctionSpaceType Space>
  Div(const ShapeFunction<NestedDerived, H1<K, Math::Vector<Number>, Mesh>, Space>&)
    -> Div<ShapeFunction<NestedDerived, H1<K, Math::Vector<Number>, Mesh>, Space>>;
}

#endif
