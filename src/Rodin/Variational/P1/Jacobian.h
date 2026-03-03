/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Jacobian.h
 * @brief Jacobian operator specialization for P1 vector functions.
 *
 * This file provides specialized implementations of the Jacobian matrix
 * operator for P1 vector-valued GridFunctions and ShapeFunctions.
 *
 * For a vector function @f$ \mathbf{u} : \Omega \to \mathbb{R}^d @f$, the
 * Jacobian is the @f$ d \times d @f$ matrix:
 * @f[
 *   \mathbf{J}(\mathbf{u}) = \begin{pmatrix}
 *     \frac{\partial u_1}{\partial x_1} & \cdots & \frac{\partial u_1}{\partial x_d} \\
 *     \vdots & \ddots & \vdots \\
 *     \frac{\partial u_d}{\partial x_1} & \cdots & \frac{\partial u_d}{\partial x_d}
 *   \end{pmatrix}
 * @f]
 *
 * For P1 elements, the Jacobian is piecewise constant on each element since
 * P1 basis function gradients are constant per element.
 *
 * @see Jacobian, P1, Grad, Div
 */
#ifndef RODIN_VARIATIONAL_P1_JACOBIAN_H
#define RODIN_VARIATIONAL_P1_JACOBIAN_H

#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Variational/Jacobian.h"
#include "Rodin/Variational/Exceptions/UndeterminedTraceDomainException.h"
#include "Rodin/Variational/Mult.h"

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
   * @brief Jacobian of a P1 vector GridFunction.
   *
   * Computes the Jacobian matrix of a P1 vector-valued grid function:
   * @f[
   *   \mathbf{J}(\mathbf{u})|_K = \sum_{i=1}^{n_v \cdot d} u_i \nabla \phi_i \otimes \mathbf{e}_j
   * @f]
   * where @f$ j = i \mod d @f$ and @f$ n_v @f$ is the number of vertices.
   *
   * The Jacobian is constant on each element for P1 functions.
   *
   * ## Applications
   * - Strain tensor in elasticity: @f$ \boldsymbol{\varepsilon} = \frac{1}{2}(\mathbf{J} + \mathbf{J}^T) @f$
   * - Velocity gradient in fluid mechanics
   * - Deformation gradient in nonlinear mechanics
   *
   * @tparam Range Value range type (typically Math::Vector<Scalar>)
   * @tparam Data Data storage type
   * @tparam Mesh Mesh type
   */
  template <class Data, class Mesh, class Scalar>
  class Jacobian<GridFunction<P1<Math::Vector<Scalar>, Mesh>, Data>> final
    : public JacobianBase<
        GridFunction<P1<Math::Vector<Scalar>, Mesh>, Data>,
        Jacobian<GridFunction<P1<Math::Vector<Scalar>, Mesh>, Data>>>
  {
    public:
      using RangeType = Math::Matrix<Scalar>;

      using FESType = P1<Math::Vector<Scalar>, Mesh>;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using SpatialMatrixType = Math::SpatialMatrix<ScalarType>;

      using OperandType = GridFunction<FESType, Data>;

      using Parent = JacobianBase<OperandType, Jacobian<OperandType>>;

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

        if (d == meshDim - 1) // face evaluation: same logic as your original
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

            if (traceDomain.size() == 0)
            {
              Alert::MemberFunctionException(*this, __func__)
                << "No trace domain provided: "
                << Alert::Notation::Predicate(true, "getTraceDomain().size() == 0")
                << ". Jacobian at an interface with no trace domain is undefined."
                << Alert::Raise;
            }

            for (auto idx : inc)
            {
              const auto& tracePolytope = mesh.getPolytope(meshDim, idx);

              const auto a = tracePolytope->getAttribute(); // Optional<Attribute>
              if (!a)
                continue; // or: treat as error (see below)

              if (traceDomain.contains(*a)) // or count(*a)
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
            return;
          }
        }

        // cell evaluation
        assert(d == mesh.getDimension());

        const auto& gf = this->getOperand();
        const auto& fes = gf.getFiniteElementSpace();
        const size_t vdim = fes.getVectorDimension();

        const auto geom = polytope.getGeometry();
        const P1Element<ScalarType> fe_scalar(geom);
        const size_t nv = fe_scalar.getCount();

        // Reference coordinates
        const auto& rc = p.getReferenceCoordinates();

        // Build the reference jacobian accumulator:
        // G = Σ_v u(v) ⊗ ∇_hat φ_v   (vdim x d)
        SpatialMatrixType G(vdim, d);
        G.setZero();

        // For each vertex basis (scalar), accumulate:
        // for each component c: G(c, :) += u_c(v) * ∇_hat φ_v^T
        for (size_t v = 0; v < nv; ++v)
        {
          // ∇_hat φ_v (size d)
          // Avoid GradientFunction() to prevent extra vector construction.
          Math::SpatialVector<ScalarType> ghat(d);
          for (size_t k = 0; k < d; ++k)
            ghat(k) = fe_scalar.getBasis(v).template getDerivative<1>(k)(rc);

          // vertex value u(v) is stored in gf as vdim dofs at that vertex
          // local index = v*vdim + c
          for (size_t c = 0; c < vdim; ++c)
          {
            const size_t local = v * vdim + c;
            const ScalarType uc = gf[fes.getGlobalIndex({d, i}, local)];

            // Row update: G(c, col) += uc * ghat(col)
            for (size_t col = 0; col < d; ++col)
              G(c, col) += uc * ghat(col);
          }
        }

        // Physical jacobian: Jx = G * J^{-1}
        out = G * p.getJacobianInverse();
      }

      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope& polytope) const noexcept
      {
        const size_t k = P1Element<ScalarType>(polytope.getGeometry()).getOrder();
        return (k == 0) ? 0 : (k - 1);
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
    static_assert(std::is_same_v<Range, Math::Vector<typename FormLanguage::Traits<Range>::ScalarType>>,
                  "Jacobian<P1> specialization is intended for vector-valued P1.");

    public:
      using FESType = P1<Range, Mesh>;
      static constexpr ShapeFunctionSpaceType SpaceType = Space;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using RangeType = Math::Matrix<ScalarType>;

      using SpatialMatrixType = Math::SpatialMatrix<ScalarType>;

      using OperandType = ShapeFunction<ShapeFunctionDerived, FESType, SpaceType>;

      using Parent = ShapeFunctionBase<Jacobian<OperandType>, FESType, SpaceType>;

      struct Cache
      {
        struct CellKey
        {
          const void* mesh = nullptr;
          size_t d = 0;
          Index i = 0;
          Geometry::Polytope::Type geom = Geometry::Polytope::Type::Point;
          int transOrder = 1;
          size_t vdim = 0;
          bool valid = false;

          explicit operator bool() const noexcept { return valid; }

          bool operator==(const CellKey& o) const noexcept
          {
            if (!valid || !o.valid) return false;
            return mesh == o.mesh && d == o.d && i == o.i
                && geom == o.geom && transOrder == o.transOrder && vdim == o.vdim;
          }

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

        struct QpKey
        {
          const QF::QuadratureFormulaBase* qf = nullptr;
          size_t qp = 0;
          bool valid = false;

          explicit operator bool() const noexcept { return valid; }

          bool operator==(const QpKey& o) const noexcept
          {
            if (!valid || !o.valid) return false;
            return qf == o.qf && qp == o.qp;
          }

          void operator=(std::initializer_list<int>) noexcept
          {
            valid = false;
            qf = nullptr;
            qp = 0;
          }
        };

        // Cached physical Jacobians, one per vector DOF: size = vdim * nvertices
        std::vector<SpatialMatrixType> jac;

        CellKey cellKey;
        QpKey qpKey;
      };

      Jacobian(const OperandType& u)
        : Parent(u.getFiniteElementSpace()),
          m_u(u),
          m_ip(nullptr)
      {}

      Jacobian(const Jacobian& other)
        : Parent(other),
          m_u(other.m_u),
          m_ip(nullptr),
          m_cache(other.m_cache)
      {}

      Jacobian(Jacobian&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)),
          m_ip(std::exchange(other.m_ip, nullptr)),
          m_cache(std::move(other.m_cache))
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

      constexpr
      const IntegrationPoint& getIntegrationPoint() const
      {
        assert(m_ip);
        return *m_ip;
      }

      static constexpr bool isTensorGeom(Geometry::Polytope::Type g) noexcept
      {
        return g == Geometry::Polytope::Type::Quadrilateral
            || g == Geometry::Polytope::Type::Wedge
            || g == Geometry::Polytope::Type::Hexahedron;
      }

      Jacobian& setIntegrationPoint(const IntegrationPoint& ip)
      {
        m_ip = &ip;

        const auto& pt   = ip.getPoint();
        const auto& poly = pt.getPolytope();
        const auto& mesh = poly.getMesh();

        const size_t d    = poly.getDimension();
        const Index  i    = poly.getIndex();
        const auto   geom = poly.getGeometry();

        const int transOrder = poly.getTransformation().getOrder();

        const auto& fes = this->getFiniteElementSpace();
        const size_t vdim = fes.getVectorDimension();

        // ---- cell key
        typename Cache::CellKey ckey;
        ckey.mesh = static_cast<const void*>(&mesh);
        ckey.d = d;
        ckey.i = i;
        ckey.geom = geom;
        ckey.transOrder = transOrder;
        ckey.vdim = vdim;
        ckey.valid = true;

        const bool cell_changed = !(m_cache.cellKey == ckey);
        if (cell_changed)
        {
          m_cache.cellKey = ckey;
          m_cache.qpKey = {}; // invalidate qp cache

          const size_t nv = Geometry::Polytope::Traits(geom).getVertexCount();
          const size_t count = vdim * nv;

          m_cache.jac.resize(count);
          for (auto& J : m_cache.jac)
          {
            // Jacobian is (vdim x d)
            J.resize(vdim, d);
            J.setZero();
          }
        }

        // ---- decide if qp needed
        const bool needs_qp = (transOrder > 1) || isTensorGeom(geom);

        typename Cache::QpKey qkey;
        if (needs_qp)
        {
          qkey.qf = &ip.getQuadratureFormula();
          qkey.qp = ip.getIndex();
          qkey.valid = true;
        }
        else
        {
          qkey.qf = nullptr;
          qkey.qp = 0;
          qkey.valid = true;
        }

        const bool qp_changed = !(m_cache.qpKey == qkey);
        if (cell_changed || qp_changed)
        {
          m_cache.qpKey = qkey;

          // Scalar P1 element gives the scalar basis derivatives
          const P1Element<ScalarType> fe_scalar(geom);
          const size_t nv = fe_scalar.getCount();

          // rc at this qp
          const auto& qf = ip.getQuadratureFormula();
          const size_t qp = ip.getIndex();
          const auto& rc = qf.getPoint(qp);

          // J^{-1} at this integration point
          const auto Jinv = pt.getJacobianInverse();

          // Precompute scalar reference gradients ghat_v (size d) for each vertex basis
          // then fill structured Jacobians.
          std::vector<Math::SpatialVector<ScalarType>> ghat;
          ghat.resize(nv);
          for (size_t v = 0; v < nv; ++v)
          {
            ghat[v].resize(d);
            for (size_t k = 0; k < d; ++k)
              ghat[v](k) = fe_scalar.getBasis(v).template getDerivative<1>(k)(rc);
          }

          // Fill each vector basis Jacobian:
          // local = v*vdim + c
          // Jhat has only row c = ghat[v]^T
          // Jphys = Jhat * Jinv
          for (size_t v = 0; v < nv; ++v)
          {
            for (size_t c = 0; c < vdim; ++c)
            {
              const size_t local = v * vdim + c;
              auto& Jphys = m_cache.jac[local];
              Jphys.setZero();

              // row c of Jhat is ghat[v]^T, so
              // row c of Jphys = ghat[v]^T * Jinv
              // i.e. (1 x d) * (d x d) -> (1 x d)
              for (size_t col = 0; col < d; ++col)
              {
                ScalarType s = ScalarType(0);
                for (size_t k = 0; k < d; ++k)
                  s += ghat[v](k) * Jinv(k, col);
                Jphys(c, col) = s;
              }
            }
          }
        }

        return *this;
      }

      auto getBasis(size_t local) const
      {
        assert(m_cache.cellKey);
        assert(local < m_cache.jac.size());
        return m_cache.jac[local].getData().topLeftCorner(
          m_cache.cellKey.vdim,
          m_cache.cellKey.d);
      }

      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope& geom) const noexcept
      {
        const auto k = getOperand().getOrder(geom);
        if (!k.has_value())
          return std::nullopt;
        return (*k == 0) ? 0 : (*k - 1);
      }

      Jacobian* copy() const noexcept override
      {
        return new Jacobian(*this);
      }

    private:
      std::reference_wrapper<const OperandType> m_u;
      const IntegrationPoint* m_ip;
      Cache m_cache;
  };

  template <class ShapeFunctionDerived, class Number, class Mesh, ShapeFunctionSpaceType Space>
  Jacobian(const ShapeFunction<ShapeFunctionDerived, P1<Math::Vector<Number>, Mesh>, Space>&)
    -> Jacobian<ShapeFunction<ShapeFunctionDerived, P1<Math::Vector<Number>, Mesh>, Space>>;
}

#endif
