/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P1_DIV_H
#define RODIN_VARIATIONAL_P1_DIV_H

/**
 * @file
 * @brief Divergence operator specialization for P1 vector functions.
 *
 * This file provides specialized implementations of the divergence operator
 * for P1 vector-valued GridFunctions and ShapeFunctions.
 *
 * For P1 elements, the divergence is computed as:
 * @f[
 *   \nabla \cdot \mathbf{u}|_K = \sum_{i=1}^{n_v} \sum_{j=1}^d u_{i,j} \frac{\partial \phi_i}{\partial x_j}
 * @f]
 * where @f$ \phi_i @f$ are P1 basis functions and @f$ u_{i,j} @f$ are DOF values.
 *
 * Since P1 basis gradients are constant on each element, the divergence is
 * also piecewise constant.
 */

#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Variational/Div.h"
#include "Rodin/Variational/IntegrationPoint.h"

namespace Rodin::FormLanguage
{
  /// @brief Type traits for @c Div over a grid function: exposes the finite element
  /// space, the scalar type and the operand type.
  template <class Scalar, class Data, class Mesh>
  struct Traits<Variational::Div<Variational::GridFunction<Variational::P1<Math::SpatialVector<Scalar>, Mesh>, Data>>>
  {
      /// @brief Finite element space type.
      using FESType = Variational::P1<Math::SpatialVector<Scalar>, Mesh>;
      /// @brief Scalar value type.
      using ScalarType = Scalar;
      /// @brief Operand type.
      using OperandType =
        Variational::GridFunction<Variational::P1<Math::SpatialVector<Scalar>, Mesh>,
          Data>;
  };

  /// @brief Type traits for @c Div over a shape function: exposes the finite element
  /// space, the shape function space, the scalar type and the operand type.
  template <class NestedDerived, class Scalar, class Mesh, Variational::ShapeFunctionSpaceType Space>
  struct Traits<
    Variational::Div<
      Variational::ShapeFunction<NestedDerived, Variational::P1<Math::SpatialVector<Scalar>, Mesh>, Space>>>
  {
      /// @brief Finite element space type.
      using FESType = Variational::P1<Math::SpatialVector<Scalar>, Mesh>;
      /// @brief Shape function space the expression belongs to, trial or test.
      static constexpr Variational::ShapeFunctionSpaceType SpaceType = Space;
      /// @brief Scalar value type.
      using ScalarType = Scalar;
      /// @brief Operand type.
      using OperandType = Variational::ShapeFunction<NestedDerived,
        Variational::P1<Math::SpatialVector<Scalar>, Mesh>, Space>;
  };
}

namespace Rodin::Variational
{
  /**
   * @ingroup DivSpecializations
   * @brief Divergence of a P1 vector GridFunction.
   *
   * Computes @f$ \nabla \cdot \mathbf{u} @f$ for vector-valued P1 functions.
   * The divergence is piecewise constant on each element:
   * @f[
   *   (\nabla \cdot \mathbf{u})|_K = \sum_{i=1}^{n_v} \mathbf{u}_i \cdot \nabla \phi_i
   * @f]
   *
   * Applications include:
   * - Incompressibility constraint: @f$ \nabla \cdot \mathbf{v} = 0 @f$
   * - Mass conservation in fluid mechanics
   * - Continuity equation
   */
  template <class Scalar, class Data, class Mesh>
  class Div<GridFunction<P1<Math::SpatialVector<Scalar>, Mesh>, Data>> final
    : public DivBase<
        GridFunction<P1<Math::SpatialVector<Scalar>, Mesh>, Data>,
        Div<GridFunction<P1<Math::SpatialVector<Scalar>, Mesh>, Data>>>
  {
    public:
      /// @brief Finite element space type.
      using FESType = P1<Math::SpatialVector<Scalar>, Mesh>;
      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;
      /// @brief Operand type.
      using OperandType = GridFunction<FESType, Data>;
      /// @brief Parent class type.
      using Parent = DivBase<OperandType, Div<OperandType>>;

      /// @brief Constructs the expression from its operand.
      Div(const OperandType& u)
        : Parent(u)
      {}

      /// @brief Copy constructor.
      Div(const Div& other)
        : Parent(other)
      {}

      /// @brief Move constructor.
      Div(Div&& other)
        : Parent(std::move(other))
      {}

      /// @brief Interpolates at an integration point.
      void interpolate(ScalarType& out, const IntegrationPoint& ip) const
      {
        interpolate(out, ip.getPoint());
      }

      /// @brief Interpolates at a geometric point.
      void interpolate(ScalarType& out, const Geometry::Point& p) const
      {
        const auto& polytope = p.getPolytope();
        const auto d = polytope.getDimension();
        const auto i = polytope.getIndex();
        const auto& mesh = polytope.getMesh();
        const size_t meshDim = mesh.getDimension();

        if (d == meshDim - 1) // face eval
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
            const auto& traceDomain = this->getTraceDomain();
            if (traceDomain.size() == 0)
            {
              Alert::MemberFunctionException(*this, __func__)
                << "No trace domain provided: "
                << Alert::Notation::Predicate(true, "getTraceDomain().size() == 0")
                << ". Div at an interface with no trace domain is undefined."
                << Alert::Raise;
            }

            for (auto& idx : inc)
            {
              const auto& tracePolytope = mesh.getPolytope(meshDim, idx);
              const auto a = tracePolytope->getAttribute();
              if (a && traceDomain.count(*a))
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
            return;
          }
        }

        // cell eval
        assert(d == mesh.getDimension());

        const auto& gf  = this->getOperand();
        const auto& fes = gf.getFiniteElementSpace();
        const size_t vdim = fes.getVectorDimension();
        assert(vdim == d); // divergence is defined as trace of Jacobian: requires vdim == spatial dim.

        const auto geom = polytope.getGeometry();
        const P1Element<ScalarType> feScalar(geom);
        const size_t nv = feScalar.getCount();

        const auto& rc = p.getReferenceCoordinates();
        const auto& Jinv = p.getJacobianInverse();

        // out = Σ_v u(v) \cdot (J^{-T} \nabla_hat φ_v)
        // but compute as: out = Σ_v Σ_j u_j(v) * (\nabla_hat φ_v)^T * (Jinv)_{:,j}
        out = ScalarType(0);

        // Precompute column sums: s_j = Σ_k ghat(k) * Jinv(k,j)  (this is (J^{-T} ghat)_j)
        for (size_t v = 0; v < nv; ++v)
        {
          // \nabla_hat φ_v
          Math::SpatialVector<ScalarType> ghat(d);
          for (size_t k = 0; k < d; ++k)
            ghat(k) = feScalar.getBasis(v).template getDerivative<1>(k)(rc);

          // phys grad components: gphys(j) = Σ_k Jinv(k,j) * ghat(k)
          // and add u(v)\cdotgphys
          for (size_t j = 0; j < d; ++j)
          {
            ScalarType gphysJ = ScalarType(0);
            for (size_t k = 0; k < d; ++k)
              gphysJ += Jinv(k, j) * ghat(k);

            const size_t local = v * vdim + j; // component j at vertex v
            out += gf[fes.getGlobalIndex({d, i}, local)] * gphysJ;
          }
        }
      }

      constexpr
      /// @brief Returns the polynomial order used on a mesh entity.
      Optional<size_t> getOrder(const Geometry::Polytope& geom) const noexcept
      {
        const size_t k = P1Element<ScalarType>(geom.getGeometry()).getOrder();
        return (k == 0) ? 0 : (k - 1);
      }

      /// @brief Creates a polymorphic copy.
      Div* copy() const noexcept override
      {
        return new Div(*this);
      }
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for Div of a P1 GridFunction
   */
  template <class Scalar, class Data, class Mesh>
  Div(const GridFunction<P1<Math::SpatialVector<Scalar>, Mesh>, Data>&)
    -> Div<GridFunction<P1<Math::SpatialVector<Scalar>, Mesh>, Data>>;

  /**
   * @ingroup DivSpecializations
   * @brief Divergence of a P1 vector ShapeFunction.
   *
   * Represents @f$ \nabla \cdot \mathbf{v} @f$ for P1 test/trial functions.
   * Used in weak formulations like:
   * @f[
   *   \int_\Omega p (\nabla \cdot \mathbf{v}) \, dx
   * @f]
   * for pressure-velocity coupling in Stokes/Navier-Stokes equations.
   */
  template <class NestedDerived, class Scalar, class Mesh, ShapeFunctionSpaceType Space>
  class Div<ShapeFunction<NestedDerived, P1<Math::SpatialVector<Scalar>, Mesh>, Space>> final
    : public ShapeFunctionBase<Div<ShapeFunction<NestedDerived, P1<Math::SpatialVector<Scalar>, Mesh>, Space>>>
  {
    public:
      /// @brief Finite element space type.
      using FESType = P1<Math::SpatialVector<Scalar>, Mesh>;
      /// @brief Shape function space the expression belongs to, trial or test.
      static constexpr ShapeFunctionSpaceType SpaceType = Space;

      /// @brief Scalar value type.
      using ScalarType = Scalar;

      /// @brief Operand type.
      using OperandType = ShapeFunction<NestedDerived, FESType, SpaceType>;

      /// @brief Parent class type.
      using Parent = ShapeFunctionBase<Div<OperandType>, FESType, SpaceType>;

      /// @brief Per-cell tabulation cache.
      struct Cache
      {
        /// @brief Key identifying the cell a tabulation was computed for.
        struct CellKey
        {
          /// @brief Mesh the cached tabulation belongs to.
          const void* mesh = nullptr;
          /// @brief Topological dimension of the cached polytope.
          size_t d = 0;
          /// @brief Index of the cached polytope.
          Index i = 0;
          /// @brief Geometry of the cached polytope.
          Geometry::Polytope::Type geom = Geometry::Polytope::Type::Point;
          /// @brief Order of the geometric transformation.
          int transOrder = 1;
          /// @brief Vector dimension of the finite element space.
          size_t vdim = 0;
          /// @brief Whether the key holds a cached entry.
          bool valid = false;

          /// @brief Tests whether the key holds a cached entry.
          explicit operator bool() const noexcept { return valid; }

          /// @brief Equality comparison.
          bool operator==(const CellKey& o) const noexcept
          {
            if (!valid || !o.valid) return false;
            return mesh == o.mesh && d == o.d && i == o.i
                && geom == o.geom && transOrder == o.transOrder && vdim == o.vdim;
          }

          /// @brief Resets the key, invalidating the cached entry.
          void operator=(std::initializer_list<int>) noexcept
          {
            valid = false;
            mesh = nullptr;
            d = 0;
            i = 0;
            geom = Geometry::Polytope::Type::Point;
            transOrder = 1;
            vdim = 0;
          }
        };

        /// @brief Key identifying the quadrature point a tabulation was computed for.
        struct QpKey
        {
          /// @brief Quadrature formula the cached tabulation belongs to.
          const QF::QuadratureFormulaBase* qf = nullptr;
          /// @brief Index of the quadrature point.
          size_t qp = 0;
          /// @brief Whether the key holds a cached entry.
          bool valid = false;

          /// @brief Tests whether the key holds a cached entry.
          explicit operator bool() const noexcept { return valid; }

          /// @brief Equality comparison.
          bool operator==(const QpKey& o) const noexcept
          {
            if (!valid || !o.valid) return false;
            return qf == o.qf && qp == o.qp;
          }

          /// @brief Resets the key, invalidating the cached entry.
          void operator=(std::initializer_list<int>) noexcept
          {
            valid = false;
            qf = nullptr;
            qp = 0;
          }
        };

        // Cached divergence of each vector basis (scalar), size = vdim * nvertices
        std::vector<ScalarType> div;

        /// @brief Key of the cached cell tabulation.
        CellKey cellKey;
        /// @brief Key of the cached quadrature-point tabulation.
        QpKey qpKey;
      };

      /// @brief Constructs the expression from its operand.
      Div(const OperandType& u)
        : Parent(u.getFiniteElementSpace()),
          m_u(u),
          m_ip(nullptr)
      {}

      /// @brief Copy constructor.
      Div(const Div& other)
        : Parent(other),
          m_u(other.m_u),
          m_ip(nullptr),
          m_cache(other.m_cache)
      {}

      /// @brief Move constructor.
      Div(Div&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)),
          m_ip(std::exchange(other.m_ip, nullptr)),
          m_cache(std::move(other.m_cache))
      {}

      constexpr
      /// @brief Gets the operand function.
      const OperandType& getOperand() const
      {
        return m_u.get();
      }

      constexpr
      /// @brief Gets the operand in the shape function expression.
      const auto& getLeaf() const
      {
        return getOperand().getLeaf();
      }

      constexpr
      /// @brief Gets the global DOF indices for a polytope.
      size_t getDOFs(const Geometry::Polytope& polytope) const
      {
        return getOperand().getDOFs(polytope);
      }

      constexpr
      /// @brief Gets the integration point the expression is evaluated at.
      const IntegrationPoint& getIntegrationPoint() const
      {
        assert(m_ip);
        return *m_ip;
      }

      /// @brief Tests whether a polytope geometry is a tensor-product one.
      static constexpr bool isTensorGeom(Geometry::Polytope::Type g) noexcept
      {
        return g == Geometry::Polytope::Type::Quadrilateral
            || g == Geometry::Polytope::Type::Wedge
            || g == Geometry::Polytope::Type::Hexahedron;
      }

      /// @brief Sets the integration point the expression is evaluated at.
      Div& setIntegrationPoint(const IntegrationPoint& ip)
      {
        m_ip = &ip;

        const auto& pt   = ip.getPoint();
        const auto& poly = pt.getPolytope();
        const auto& mesh = poly.getMesh();

        const size_t d    = poly.getDimension();
        const Index  i    = poly.getIndex();
        const auto   geom = poly.getGeometry();

        const int transOrder = poly.getTransformation().getOrder();
        const auto* qf = ip.getQuadratureFormula();

        const auto& fes = this->getFiniteElementSpace();
        const size_t vdim = fes.getVectorDimension();
        assert(vdim == d);

        // cell key
        typename Cache::CellKey ckey;
        ckey.mesh = static_cast<const void*>(&mesh);
        ckey.d = d;
        ckey.i = i;
        ckey.geom = geom;
        ckey.transOrder = transOrder;
        ckey.vdim = vdim;
        ckey.valid = true;

        const bool cellChanged = !(m_cache.cellKey == ckey);
        if (cellChanged)
        {
          m_cache.cellKey = ckey;
          m_cache.qpKey = {};

          const size_t nv = Geometry::Polytope::Traits(geom).getVertexCount();
          m_cache.div.resize(vdim * nv);
        }

        const bool needsQp = (transOrder > 1) || isTensorGeom(geom);

        typename Cache::QpKey qkey;
        if (needsQp)
        {
          qkey.qf = qf;
          qkey.qp = qf ? ip.getIndex() : 0;
          qkey.valid = true;
        }
        else
        {
          qkey.qf = nullptr;
          qkey.qp = 0;
          qkey.valid = true;
        }

        const bool qpChanged = !qf || !(m_cache.qpKey == qkey);
        if (cellChanged || qpChanged)
        {
          m_cache.qpKey = qkey;

          const P1Element<ScalarType> feScalar(geom);
          const size_t nv = feScalar.getCount();

          const auto& rc =
            qf ? qf->getPoint(ip.getIndex()) : pt.getReferenceCoordinates();

          const auto& Jinv = pt.getJacobianInverse();

          // For basis (v,c): div(φ_{v,c}) = (J^{-T} \nabla_hat φ_v)_c
          // because only component c is nonzero.
          for (size_t v = 0; v < nv; ++v)
          {
            Math::SpatialVector<ScalarType> ghat;
            ghat.resize(d);
            for (size_t k = 0; k < d; ++k)
              ghat(k) = feScalar.getBasis(v).template getDerivative<1>(k)(rc);

            for (size_t c = 0; c < d; ++c)
            {
              ScalarType gphysC = ScalarType(0);
              for (size_t k = 0; k < d; ++k)
                gphysC += Jinv(k, c) * ghat(k);

              m_cache.div[v * vdim + c] = gphysC;
            }
          }
        }

        return *this;
      }

      /// @brief Gets the basis function of a local degree of freedom.
      ScalarType getBasis(size_t local) const
      {
        assert(m_cache.cellKey);
        assert(local < m_cache.div.size());
        return m_cache.div[local];
      }

      constexpr
      /// @brief Returns the polynomial order used on a mesh entity.
      Optional<size_t> getOrder(const Geometry::Polytope& geom) const noexcept
      {
        const auto k = getOperand().getOrder(geom);
        if (!k.has_value())
          return std::nullopt;
        return (*k == 0) ? 0 : (*k - 1);
      }

      Div* copy() const noexcept override
      {
        return new Div(*this);
      }

    private:
      std::reference_wrapper<const OperandType> m_u;
      const IntegrationPoint* m_ip;
      Cache m_cache;
  };

  /// @brief Deduction guide for @c Div.
  template <class NestedDerived, class Number, class Mesh, ShapeFunctionSpaceType Space>
  Div(const ShapeFunction<NestedDerived, P1<Math::SpatialVector<Number>, Mesh>, Space>&)
    -> Div<ShapeFunction<NestedDerived, P1<Math::SpatialVector<Number>, Mesh>, Space>>;
}

#endif
