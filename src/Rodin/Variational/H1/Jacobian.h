/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Jacobian.h
 * @brief Jacobian operator specialization for H1 vector functions.
 *
 * This file provides specialized implementations of the Jacobian matrix
 * operator for H1 vector-valued GridFunctions and ShapeFunctions.
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
 * For H1<K> elements, the Jacobian is polynomial of degree K-1 on each element.
 *
 * @see Jacobian, H1, Grad, Div
 */
#ifndef RODIN_VARIATIONAL_H1_JACOBIAN_H
#define RODIN_VARIATIONAL_H1_JACOBIAN_H

#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Variational/Jacobian.h"
#include "Rodin/Variational/Exceptions/UndeterminedTraceDomainException.h"

namespace Rodin::FormLanguage
{
  template <size_t K, class Range, class Data, class Mesh>
  struct Traits<
    Variational::Jacobian<
      Variational::GridFunction<
        Variational::H1<K, Range, Mesh>, Data>>>
  {
    using FESType = Variational::H1<K, Range, Mesh>;
    using OperandType = Variational::GridFunction<FESType, Data>;
  };

  template <size_t K, class NestedDerived, class Range, class Mesh, Variational::ShapeFunctionSpaceType Space>
  struct Traits<
    Variational::Jacobian<
      Variational::ShapeFunction<NestedDerived, Variational::H1<K, Range, Mesh>, Space>>>
  {
    using FESType = Variational::H1<K, Range, Mesh>;
    static constexpr Variational::ShapeFunctionSpaceType SpaceType = Space;
    using OperandType = Variational::ShapeFunction<NestedDerived, FESType, Space>;
  };
}

namespace Rodin::Variational
{
  /**
   * @ingroup JacobianSpecializations
   * @brief Jacobian of an H1 vector GridFunction.
   *
   * Computes the Jacobian matrix of an H1<K> vector-valued grid function:
   * @f[
   *   \mathbf{J}(\mathbf{u})|_K = \sum_{i=1}^{n_{dof} \cdot d} u_i \nabla \phi_i \otimes \mathbf{e}_j
   * @f]
   * where @f$ j = i \mod d @f$ and @f$ n_{dof} @f$ is the number of DOFs per element.
   *
   * The Jacobian is polynomial of degree K-1 on each element for H1<K> functions.
   *
   * ## Applications
   * - Strain tensor in elasticity: @f$ \boldsymbol{\varepsilon} = \frac{1}{2}(\mathbf{J} + \mathbf{J}^T) @f$
   * - Velocity gradient in fluid mechanics
   * - Deformation gradient in nonlinear mechanics
   *
   * @tparam K Polynomial degree
   * @tparam Range Value range type (typically Math::Vector<Scalar>)
   * @tparam Data Data storage type
   * @tparam Mesh Mesh type
   */
  template <size_t K, class Scalar, class Data, class Mesh>
  class Jacobian<GridFunction<H1<K, Math::Vector<Scalar>, Mesh>, Data>> final
    : public JacobianBase<
        GridFunction<H1<K, Math::Vector<Scalar>, Mesh>, Data>,
        Jacobian<GridFunction<H1<K, Math::Vector<Scalar>, Mesh>, Data>>>
  {
    public:
      using FESType = H1<K, Math::Vector<Scalar>, Mesh>;

      using OperandType = GridFunction<FESType, Data>;

      using Parent = JacobianBase<OperandType, Jacobian<OperandType>>;

      using RangeType = Math::Matrix<Scalar>;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using SpatialMatrixType = Math::SpatialMatrix<ScalarType>;

      using SpatialVectorType = Math::SpatialVector<ScalarType>;

      Jacobian(const OperandType& u) : Parent(u) {}
      Jacobian(const Jacobian& other) : Parent(other) {}
      Jacobian(Jacobian&& other) : Parent(std::move(other)) {}

      void interpolate(SpatialMatrixType& out, const Geometry::Point& p) const
      {
        const auto& polytope = p.getPolytope();
        const auto& d = polytope.getDimension();
        const auto& i = polytope.getIndex();
        const auto& mesh = polytope.getMesh();
        const size_t meshDim = mesh.getDimension();

        if (d == meshDim - 1) // face
        {
          const auto& conn = mesh.getConnectivity();
          const auto& inc  = conn.getIncidence({ meshDim - 1, meshDim }, i);
          const auto& pc   = p.getPhysicalCoordinates();

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
            const auto& traceDomain = this->getTraceDomain();
            if (traceDomain.size() == 0)
            {
              Alert::MemberFunctionException(*this, __func__)
                << "No trace domain provided: "
                << Alert::Notation::Predicate(true, "getTraceDomain().size() == 0")
                << ". Jacobian at an interface with no trace domain is undefined."
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
                this->interpolate(out, np);
                return;
              }
            }

            UndeterminedTraceDomainException(
                *this, __func__, {d, i}, traceDomain.begin(), traceDomain.end()) << Alert::Raise;
            return;
          }
        }

        // cell
        assert(d == mesh.getDimension());

        const auto& gf  = this->getOperand();
        const auto& fes = gf.getFiniteElementSpace();
        const size_t vdim = fes.getVectorDimension();

        const auto geom = polytope.getGeometry();
        const auto& rc  = p.getReferenceCoordinates();

        // scalar element gives scalar basis gradients in reference coords
        const auto feS = H1Element<K, ScalarType>(geom);
        const size_t nscalar = feS.getCount();

        // phys gradient mapping
        const auto JinvT = p.getJacobianInverse().transpose();

        static thread_local SpatialVectorType s_ref;
        static thread_local SpatialVectorType s_phys;
        s_ref.resize(d);
        s_phys.resize(d);

        SpatialMatrixType res(vdim, d);
        res.setZero();

        // Local ordering must match your vector H1 element:
        // local = alpha * vdim + comp
        for (size_t alpha = 0; alpha < nscalar; ++alpha)
        {
          const auto gref = feS.getBasis(alpha).getGradient()(rc); // length d

          for (size_t j = 0; j < d; ++j)
            s_ref(j) = gref(j);

          s_phys = JinvT * s_ref; // ∇_x φ_alpha

          for (size_t comp = 0; comp < vdim; ++comp)
          {
            const size_t local = alpha * vdim + comp;

            // J(u) row 'comp' gets contribution u_comp(alpha) * (∇φ_alpha)^T
            const auto uval = gf[fes.getGlobalIndex({d, i}, local)];
            if (comp < vdim)
            {
              for (size_t j = 0; j < d; ++j)
                res(comp, j) += uval * s_phys(j);
            }
          }
        }

        out = std::move(res);
      }

      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope& polytope) const noexcept
      {
        const size_t k = H1Element<K, ScalarType>(polytope.getGeometry()).getOrder();
        return (k == 0) ? 0 : (k - 1);
      }

      Jacobian* copy() const noexcept override { return new Jacobian(*this); }
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for Jacobian of an H1 GridFunction
   */
  template <size_t K, class Range, class Data, class Mesh>
  Jacobian(const GridFunction<H1<K, Range, Mesh>, Data>&) -> Jacobian<GridFunction<H1<K, Range, Mesh>, Data>>;

  /**
   * @ingroup JacobianSpecializations
   * @brief Jacobian of an H1 ShapeFunction object.
   */
  template <size_t K, class NestedDerived, class Scalar, class Mesh, ShapeFunctionSpaceType Space>
  class Jacobian<ShapeFunction<NestedDerived, H1<K, Math::Vector<Scalar>, Mesh>, Space>> final
    : public ShapeFunctionBase<
        Jacobian<ShapeFunction<NestedDerived, H1<K, Math::Vector<Scalar>, Mesh>, Space>>,
        H1<K, Math::Vector<Scalar>, Mesh>,
        Space>
  {
    public:
      using FESType = H1<K, Math::Vector<Scalar>, Mesh>;
      static constexpr ShapeFunctionSpaceType SpaceType = Space;

      using OperandType = ShapeFunction<NestedDerived, FESType, SpaceType>;

      using Parent      =
        ShapeFunctionBase<
          Jacobian<ShapeFunction<NestedDerived, FESType, SpaceType>>,
          FESType,
          SpaceType>;

      using ScalarType        = typename FormLanguage::Traits<FESType>::ScalarType;

      using RangeType         = Math::Matrix<ScalarType>;

      using SpatialMatrixType = Math::SpatialMatrix<ScalarType>;

      using SpatialVectorType = Math::SpatialVector<ScalarType>;

      struct Cache
      {
        struct Key
        {
          Geometry::Polytope::Type geom = Geometry::Polytope::Type::Point;
          size_t dim = 0;
          Index cell = 0;

          const QF::QuadratureFormulaBase* qf = nullptr;
          size_t qp = 0;

          bool valid = false;

          explicit operator bool() const noexcept { return valid; }

          bool operator==(const Key& o) const noexcept
          {
            if (!valid || !o.valid)
              return false;
            return geom == o.geom
                && dim  == o.dim
                && cell == o.cell
                && qf   == o.qf
                && qp   == o.qp;
          }

          void invalidate() noexcept
          {
            valid = false;
            geom = Geometry::Polytope::Type::Point;
            dim = 0;
            cell = 0;
            qf = nullptr;
            qp = 0;
          }
        };

        // minimal cache: physical gradients for scalar basis indices alpha
        std::vector<SpatialVectorType> grad_phys; // size = nscalar, each size = dim
        Key key;
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
      const auto& getFiniteElementSpace() const
      {
        return getOperand().getFiniteElementSpace();
      }

      constexpr
      const IntegrationPoint& getIntegrationPoint() const
      {
        assert(m_ip);
        return *m_ip;
      }

      Jacobian& setIntegrationPoint(const IntegrationPoint& ip)
      {
        m_ip = &ip;

        const auto& p  = ip.getPoint();
        const auto& qf = ip.getQuadratureFormula();
        const size_t qp = ip.getIndex();

        const auto& poly = p.getPolytope();
        const auto  geom = poly.getGeometry();
        const size_t d   = poly.getDimension();
        const Index  cell = poly.getIndex();

        typename Cache::Key key;
        key.geom  = geom;
        key.dim   = d;
        key.cell  = cell;
        key.qf    = &qf;
        key.qp    = qp;
        key.valid = true;

        if (m_cache.key == key)
          return *this;

        m_cache.key = key;

        // scalar element for geometry-only tabulation
        const H1Element<K, ScalarType> feS(geom);
        const size_t nscalar = feS.getCount();

        if (m_cache.grad_phys.size() != nscalar)
          m_cache.grad_phys.resize(nscalar);

        for (auto& g : m_cache.grad_phys)
          if (g.size() != d) g.resize(d);

        const auto& tab   = feS.getTabulation(qf);
        const auto  JinvT = p.getJacobianInverse().transpose();

        static thread_local SpatialVectorType s_ref;
        s_ref.resize(d);

        for (size_t alpha = 0; alpha < nscalar; ++alpha)
        {
          const auto gref = tab.getGradient(qp, alpha); // span size d
          for (size_t j = 0; j < d; ++j)
            s_ref(j) = gref[j];

          m_cache.grad_phys[alpha] = JinvT * s_ref;
        }

        return *this;
      }

      const RangeType& getBasis(size_t local) const
      {
        assert(m_cache.key);

        const auto& fes  = this->getFiniteElementSpace();
        const size_t vdim = fes.getVectorDimension();

        const auto& p = this->getIntegrationPoint().getPoint();
        const size_t d = p.getPolytope().getDimension();

        const size_t alpha = local / vdim;
        const size_t comp  = local % vdim;

        assert(alpha < m_cache.grad_phys.size());

        static thread_local RangeType s_J;
        if (static_cast<size_t>(s_J.rows()) != vdim || static_cast<size_t>(s_J.cols()) != d)
          s_J.resize(vdim, d);

        s_J.setZero();

        // Only row comp is non-zero
        for (size_t j = 0; j < d; ++j)
          s_J(comp, j) = m_cache.grad_phys[alpha](j);

        return s_J;
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

  template <size_t K, class ShapeFunctionDerived, class Number, class Mesh, ShapeFunctionSpaceType Space>
  Jacobian(const ShapeFunction<ShapeFunctionDerived, H1<K, Math::Vector<Number>, Mesh>, Space>&)
    -> Jacobian<ShapeFunction<ShapeFunctionDerived, H1<K, Math::Vector<Number>, Mesh>, Space>>;
}

#endif
