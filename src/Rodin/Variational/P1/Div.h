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
  /**
   * @brief Type traits for @c Div over a grid function: exposes the finite element
   * space, the scalar type and the operand type.
   */
  template <class Scalar, class Data, class Mesh>
  struct Traits<Variational::Div<
    Variational::GridFunction<Variational::P1<Math::SpatialVector<Scalar>, Mesh>, Data>>>
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

  /**
   * @brief Type traits for @c Div over a shape function: exposes the finite element
   * space, the shape function space, the scalar type and the operand type.
   */
  template <class NestedDerived, class Scalar, class Mesh,
    Variational::ShapeFunctionSpaceType Space>
  struct Traits<Variational::Div<Variational::ShapeFunction<NestedDerived,
    Variational::P1<Math::SpatialVector<Scalar>, Mesh>, Space>>>
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
    : public DivBase<GridFunction<P1<Math::SpatialVector<Scalar>, Mesh>, Data>,
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

      /**
       * @brief Constructs the expression from its operand.
       * @param u Operand expression.
       */
      Div(const OperandType& u)
        : Parent(u)
      {}

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      Div(const Div& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      Div(Div&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Interpolates at an integration point.
       * @param out Storage for the computed result.
       * @param ip Integration point at which the expression is evaluated.
       */
      void interpolate(ScalarType& out, const IntegrationPoint& ip) const
      {
        const auto& p = ip.getPoint();
        const auto& polytope = p.getPolytope();
        const auto* qf = ip.getReferenceQuadratureFormula();
        if (!qf || polytope.getDimension() != polytope.getMesh().getDimension())
        {
          interpolate(out, p);
          return;
        }
        const size_t d = polytope.getDimension();
        const Index i = polytope.getIndex();

        const auto& gf = this->getOperand();
        const auto& fes = gf.getFiniteElementSpace();
        const size_t vdim = fes.getVectorDimension();
        const size_t spaceDimension = polytope.getMesh().getSpaceDimension();
        if (vdim != spaceDimension)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Divergence requires one vector component per ambient coordinate."
            << Alert::Raise;
        }

        const P1Element<ScalarType> feS(polytope.getGeometry());
        const size_t nscalar = feS.getCount();
        const auto JinvT = p.getJacobianInverse().transpose();

        Math::SpatialVector<ScalarType> ref(static_cast<std::uint8_t>(d));
        Math::SpatialVector<ScalarType> phys(spaceDimension);

        out = ScalarType(0);
        for (size_t alpha = 0; alpha < nscalar; ++alpha)
        {
          const auto gref = feS.getTabulation(*qf).getGradient(ip.getIndex(), alpha);
          for (size_t j = 0; j < d; ++j)
            ref(static_cast<std::uint8_t>(j)) = gref[j];

          phys = JinvT * ref;

          for (size_t comp = 0; comp < vdim; ++comp)
          {
            const size_t local = alpha * vdim + comp;
            out += gf[fes.getGlobalIndex({d, i}, local)] *
              phys(static_cast<std::uint8_t>(comp));
          }
        }
      }

      /**
       * @brief Interpolates at a geometric point.
       * @param out Storage for the computed result.
       * @param p Point at which the operation is evaluated.
       */
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
          const auto& inc = conn.getIncidence({meshDim - 1, meshDim}, i);
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
              *this, __func__, {d, i}, traceDomain.begin(), traceDomain.end())
              << Alert::Raise;
            return;
          }
        }

        // cell eval
        assert(d == mesh.getDimension());

        const auto& gf = this->getOperand();
        const auto& fes = gf.getFiniteElementSpace();
        const size_t vdim = fes.getVectorDimension();
        if (vdim != mesh.getSpaceDimension())
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Divergence requires one vector component per ambient coordinate."
            << Alert::Raise;
        }

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
          for (size_t j = 0; j < vdim; ++j)
          {
            ScalarType gphysJ = ScalarType(0);
            for (size_t k = 0; k < d; ++k)
              gphysJ += Jinv(k, j) * ghat(k);

            const size_t local = v * vdim + j; // component j at vertex v
            out += gf[fes.getGlobalIndex({d, i}, local)] * gphysJ;
          }
        }
      }

      /**
       * @brief Returns the polynomial order used on a mesh entity.
       * @param geom Reference geometry.
       * @returns Polynomial order on the entity, or an empty optional when no order is available.
       */
      constexpr Optional<size_t> getOrder(const Geometry::Polytope& geom) const noexcept
      {
        const size_t k = P1Element<ScalarType>(geom.getGeometry()).getOrder();
        return (k == 0) ? 0 : (k - 1);
      }

      /**
       * @brief Creates a polymorphic copy.
       * @returns Pointer to a newly allocated copy; the caller owns the returned object.
       */
      Div* copy() const noexcept override
      {
        return new Div(*this);
      }
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for Div of a P1 GridFunction
   * @param u Operand expression.
   */
  template <class Scalar, class Data, class Mesh>
  Div(const GridFunction<P1<Math::SpatialVector<Scalar>, Mesh>, Data>& u)
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
  class Div<ShapeFunction<NestedDerived, P1<Math::SpatialVector<Scalar>, Mesh>, Space>>
    final
    : public ShapeFunctionBase<
        Div<ShapeFunction<NestedDerived, P1<Math::SpatialVector<Scalar>, Mesh>, Space>>>
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
          const Geometry::MeshBase* mesh;

          /// @brief Topological dimension of the cached polytope.
          size_t dimension;

          /// @brief Index of the cached polytope.
          Index index;

          /// @brief Geometry of the cached polytope.
          Geometry::Polytope::Type geometry;

          /// @brief Order of the geometric transformation.
          size_t transformationOrder;

          /// @brief Vector dimension of the finite element space.
          size_t vectorDimension;

          /// @brief Whether the key holds a cached entry.
          bool valid;

          /// @brief Constructs an invalid cell key.
          CellKey()
            : mesh(nullptr),
              dimension(0),
              index(0),
              geometry(Geometry::Polytope::Type::Point),
              transformationOrder(1),
              vectorDimension(0),
              valid(false)
          {}

          /**
           * @brief Tests whether the key holds a cached entry.
           * @returns True if the key identifies a cached entry; false otherwise.
           */
          explicit operator bool() const noexcept
          {
            return valid;
          }

          /**
           * @brief Equality comparison.
           * @returns Whether the operands compare equal.
           * @param o Key to compare with this key.
           */
          bool operator==(const CellKey& o) const noexcept
          {
            if (!valid || !o.valid)
              return false;
            return mesh == o.mesh && dimension == o.dimension && index == o.index &&
              geometry == o.geometry &&
              transformationOrder == o.transformationOrder &&
              vectorDimension == o.vectorDimension;
          }

          /// @brief Resets the key, invalidating the cached entry.
          void reset() noexcept
          {
            valid = false;
            mesh = nullptr;
            dimension = 0;
            index = 0;
            geometry = Geometry::Polytope::Type::Point;
            transformationOrder = 1;
            vectorDimension = 0;
          }
        };

        /// @brief Key identifying the quadrature point a tabulation was computed for.
        struct QuadraturePointKey
        {
          /// @brief Quadrature formula the cached tabulation belongs to.
          const QF::QuadratureFormula* formula;

          /// @brief Index of the quadrature point.
          size_t index;

          /// @brief Whether the key holds a cached entry.
          bool valid;

          /// @brief Constructs an invalid quadrature cache key.
          QuadraturePointKey()
            : formula(nullptr),

              index(0),
              valid(false)
          {}

          /**
           * @brief Tests whether the key holds a cached entry.
           * @returns True if the key identifies a cached entry; false otherwise.
           */
          explicit operator bool() const noexcept
          {
            return valid;
          }

          /**
           * @brief Equality comparison.
           * @returns Whether the operands compare equal.
           * @param o Key to compare with this key.
           */
          bool operator==(const QuadraturePointKey& o) const noexcept
          {
            if (!valid || !o.valid)
              return false;
            return formula == o.formula && index == o.index;
          }

          /// @brief Resets the key, invalidating the cached entry.
          void reset() noexcept
          {
            valid = false;
            formula = nullptr;

            index = 0;
          }
        };

        // Cached divergence of each vector basis (scalar), size = vdim * nvertices
        /// @brief Cached divergence values.
        std::vector<ScalarType> div;

        /// @brief Key of the cached cell tabulation.
        CellKey cellKey;

        /// @brief Key of the cached quadrature-point tabulation.
        QuadraturePointKey quadraturePointKey;
      };

      /**
       * @brief Constructs the expression from its operand.
       * @param u Operand expression.
       */
      Div(const OperandType& u)
        : Parent(u.getFiniteElementSpace()),
          m_u(u),
          m_ip(nullptr)
      {}

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      Div(const Div& other)
        : Parent(other),
          m_u(other.m_u),
          m_ip(nullptr)
      {}

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      Div(Div&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)),
          m_ip(nullptr)
      {}

      /**
       * @brief Gets the operand function.
       * @returns The operand function.
       */
      constexpr const OperandType& getOperand() const
      {
        return m_u.get();
      }

      /**
       * @brief Gets the operand in the shape function expression.
       * @returns The operand in the shape function expression.
       */
      constexpr const auto& getLeaf() const
      {
        return getOperand().getLeaf();
      }

      /**
       * @brief Returns the number of local basis functions for a polytope.
       * @param polytope Mesh entity used by this operation.
       * @returns Number of local basis functions on the selected entity.
       */
      constexpr size_t getDOFs(const Geometry::Polytope& polytope) const
      {
        return getOperand().getDOFs(polytope);
      }

      /**
       * @brief Gets the integration point the expression is evaluated at.
       * @returns The integration point the expression is evaluated at.
       */
      constexpr const IntegrationPoint& getIntegrationPoint() const
      {
        assert(m_ip);
        return *m_ip;
      }

      /**
       * @brief Tests whether a polytope geometry is a tensor-product one.
       * @param g Mesh entity used by this operation.
       * @returns Whether the geometry is a tensor-product reference cell.
       */
      static constexpr bool isTensorGeom(Geometry::Polytope::Type g) noexcept
      {
        return g == Geometry::Polytope::Type::Quadrilateral ||
          g == Geometry::Polytope::Type::Wedge ||
          g == Geometry::Polytope::Type::Hexahedron;
      }

      /**
       * @brief Starts a fresh cell binding for physical basis derivatives.
       *
       * Affine simplex derivatives are prepared at the first point and reused
       * throughout this binding. Curved or varying reference derivatives are
       * prepared at each new quadrature sample. Rebinding the same entity
       * discards the preceding physical results while retaining storage.
       * @param polytope Entity whose geometry remains fixed during this binding.
       * @returns This expression with no point bound.
       */
      Div& setPolytope(const Geometry::Polytope& polytope)
      {
        m_ip = nullptr;
        m_binding.reset();
        m_cache.cellKey.reset();
        m_cache.quadraturePointKey.reset();
        m_binding.mesh = &polytope.getMesh();
        m_binding.dimension = polytope.getDimension();
        m_binding.index = polytope.getIndex();
        m_binding.geometry = polytope.getGeometry();
        m_binding.transformationOrder = polytope.getTransformation().getOrder();
        m_binding.vectorDimension = this->getFiniteElementSpace().getVectorDimension();
        m_binding.valid = true;
        return *this;
      }

      /**
       * @brief Sets the integration point the expression is evaluated at.
       *
       * Pyramidal scalar basis derivatives are rational and sample-dependent.
       * Even with a constant inverse map, their contractions defining basis
       * divergence must be recomputed at each reference sample. Tensor-product
       * and pyramidal reference derivatives retain sample-indexed reuse.
       *
       * @param ip Integration point at which the expression is evaluated.
       * @returns Reference to this object after the operation.
       */
      Div& setIntegrationPoint(const IntegrationPoint& ip)
      {
        m_ip = &ip;

        const auto& pt = ip.getPoint();
        const auto& poly = pt.getPolytope();
        const auto& mesh = poly.getMesh();

        const size_t d = poly.getDimension();
        const Index i = poly.getIndex();
        const auto geom = poly.getGeometry();

        const auto& transformation = poly.getTransformation();
        const size_t transformationOrder = transformation.getOrder();
        const auto* qf = ip.getReferenceQuadratureFormula();

        const auto& fes = this->getFiniteElementSpace();
        const size_t vdim = fes.getVectorDimension();
        if (vdim != mesh.getSpaceDimension())
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Divergence requires one vector component per ambient coordinate."
            << Alert::Raise;
        }

        // cell key
        typename Cache::CellKey cellKey;
        cellKey.mesh = &mesh;
        cellKey.dimension = d;
        cellKey.index = i;
        cellKey.geometry = geom;
        cellKey.transformationOrder = transformationOrder;
        cellKey.vectorDimension = vdim;
        cellKey.valid = true;

        const bool cellChanged = !(m_binding == cellKey) || !(m_cache.cellKey == cellKey);
        if (cellChanged)
        {
          m_cache.cellKey = cellKey;
          m_cache.quadraturePointKey.reset();

          const size_t nv = Geometry::Polytope::Traits(geom).getVertexCount();
          m_cache.div.resize(vdim * nv);
        }

        const bool needsQp = (transformationOrder > 1) || isTensorGeom(geom) ||
          geom == Geometry::Polytope::Type::Pyramid;

        typename Cache::QuadraturePointKey quadraturePointKey;
        if (needsQp)
        {
          quadraturePointKey.formula = qf;
          quadraturePointKey.index = qf ? ip.getIndex() : 0;
          quadraturePointKey.valid = true;
        }
        else
        {
          quadraturePointKey.formula = nullptr;
          quadraturePointKey.index = 0;
          quadraturePointKey.valid = true;
        }

        const bool quadraturePointChanged =
          (needsQp && !qf) || !(m_cache.quadraturePointKey == quadraturePointKey);
        if (cellChanged || quadraturePointChanged)
        {
          m_cache.quadraturePointKey = quadraturePointKey;

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

            for (size_t c = 0; c < vdim; ++c)
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

      /**
       * @brief Gets the basis function of a local degree of freedom.
       * @param local Index in the local numbering.
       * @returns Value of the selected local basis function at the evaluation point.
       */
      ScalarType getBasis(size_t local) const
      {
        assert(m_cache.cellKey);
        assert(local < m_cache.div.size());
        return m_cache.div[local];
      }

      /**
       * @brief Returns the polynomial order used on a mesh entity.
       * @param geom Reference geometry.
       * @returns Polynomial order on the entity, or an empty optional when no order is available.
       */
      constexpr Optional<size_t> getOrder(const Geometry::Polytope& geom) const noexcept
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
      typename Cache::CellKey m_binding;
  };

  /**
   * @brief Deduction guide for @c Div.
   * @param u Operand expression.
   */
  template <class NestedDerived, class Number, class Mesh, ShapeFunctionSpaceType Space>
  Div(const ShapeFunction<NestedDerived, P1<Math::SpatialVector<Number>, Mesh>, Space>& u)
    -> Div<ShapeFunction<NestedDerived, P1<Math::SpatialVector<Number>, Mesh>, Space>>;
}

#endif
