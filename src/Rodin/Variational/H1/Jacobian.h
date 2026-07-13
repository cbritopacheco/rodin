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

#include <vector>

#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Variational/IntegrationPoint.h"
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
    /// @brief Finite element space type.
      using FESType = Variational::H1<K, Range, Mesh>;
    /// @brief Operand type.
      using OperandType = Variational::GridFunction<FESType, Data>;
  };

  template <size_t K, class NestedDerived, class Range, class Mesh, Variational::ShapeFunctionSpaceType Space>
  struct Traits<
    Variational::Jacobian<
      Variational::ShapeFunction<NestedDerived, Variational::H1<K, Range, Mesh>, Space>>>
  {
    /// @brief Finite element space type.
      using FESType = Variational::H1<K, Range, Mesh>;
      static constexpr Variational::ShapeFunctionSpaceType SpaceType = Space;
    /// @brief Operand type.
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
   * @tparam Range Value range type (typically Math::SpatialVector<Scalar>)
   * @tparam Data Data storage type
   * @tparam Mesh Mesh type
   */
  template <size_t K, class Scalar, class Data, class Mesh>
  class Jacobian<GridFunction<H1<K, Math::SpatialVector<Scalar>, Mesh>, Data>> final
    : public JacobianBase<
        GridFunction<H1<K, Math::SpatialVector<Scalar>, Mesh>, Data>,
        Jacobian<GridFunction<H1<K, Math::SpatialVector<Scalar>, Mesh>, Data>>>
  {
    public:
      /// @brief Finite element space type.
      using FESType = H1<K, Math::SpatialVector<Scalar>, Mesh>;

      /// @brief Operand type.
      using OperandType = GridFunction<FESType, Data>;

      /// @brief Parent class type.
      using Parent = JacobianBase<OperandType, Jacobian<OperandType>>;

      /// @brief Range (evaluation value) type.
      using RangeType = Math::SpatialMatrix<Scalar>;

      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using SpatialMatrixType = Math::SpatialMatrix<ScalarType>;

      using SpatialVectorType = Math::SpatialVector<ScalarType>;

      Jacobian(const OperandType& u) : Parent(u) {}

      Jacobian(const Jacobian& other)
        : Parent(other)
      {}

      Jacobian(Jacobian&& other)
        : Parent(std::move(other))
      {}

      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope& polytope) const noexcept
      {
        const size_t k = H1Element<K, ScalarType>(polytope.getGeometry()).getOrder();
        return (k == 0) ? 0 : (k - 1);
      }

      Jacobian* copy() const noexcept override { return new Jacobian(*this); }

      void interpolate(SpatialMatrixType& out, const IntegrationPoint& ip) const
      {
        const auto& p = ip.getPoint();
        const auto& polytope = p.getPolytope();
        const size_t d = polytope.getDimension();
        const Index i = polytope.getIndex();

        const auto& gf  = this->getOperand();
        const auto& fes = gf.getFiniteElementSpace();
        const size_t vdim = fes.getVectorDimension();

        const auto feS = H1Element<K, ScalarType>(polytope.getGeometry());
        const size_t nscalar = feS.getCount();
        const auto* qf = ip.getQuadratureFormula();
        assert(qf);
        const auto& tab = feS.getTabulation(*qf);
        const auto JinvT = p.getJacobianInverse().transpose();

        SpatialVectorType ref(d);
        SpatialVectorType phys(d);

        out.resize(vdim, d);
        out.setZero();

        for (size_t alpha = 0; alpha < nscalar; ++alpha)
        {
          const auto gref = tab.getGradient(ip.getIndex(), alpha);
          for (size_t j = 0; j < d; ++j)
            ref(j) = gref[j];

          phys = JinvT * ref;

          for (size_t comp = 0; comp < vdim; ++comp)
          {
            const size_t local = alpha * vdim + comp;
            const auto uval = gf[fes.getGlobalIndex({d, i}, local)];
            for (size_t j = 0; j < d; ++j)
              out(comp, j) += uval * phys(j);
          }
        }
      }

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
                interpolate(out, np);
                return;
              }
            }

            UndeterminedTraceDomainException(
                *this, __func__, {d, i}, traceDomain.begin(), traceDomain.end()) << Alert::Raise;
            return;
          }
        }

        assert(d == mesh.getDimension());

        const auto& gf  = this->getOperand();
        const auto& fes = gf.getFiniteElementSpace();
        const size_t vdim = fes.getVectorDimension();

        const auto geom = polytope.getGeometry();
        const auto feS = H1Element<K, ScalarType>(geom);
        const size_t nscalar = feS.getCount();
        const auto& rc = p.getReferenceCoordinates();
        const auto JinvT = p.getJacobianInverse().transpose();

        SpatialVectorType ref(d);
        SpatialVectorType phys(d);

        out.resize(vdim, d);
        out.setZero();

        for (size_t alpha = 0; alpha < nscalar; ++alpha)
        {
          const auto gref = feS.getBasis(alpha).getGradient()(rc);
          for (size_t j = 0; j < d; ++j)
            ref(j) = gref(j);

          phys = JinvT * ref;

          for (size_t comp = 0; comp < vdim; ++comp)
          {
            const size_t local = alpha * vdim + comp;
            const auto uval = gf[fes.getGlobalIndex({d, i}, local)];
            for (size_t j = 0; j < d; ++j)
              out(comp, j) += uval * phys(j);
          }
        }
      }
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
  class Jacobian<ShapeFunction<NestedDerived, H1<K, Math::SpatialVector<Scalar>, Mesh>, Space>> final
    : public ShapeFunctionBase<
        Jacobian<ShapeFunction<NestedDerived, H1<K, Math::SpatialVector<Scalar>, Mesh>, Space>>,
        H1<K, Math::SpatialVector<Scalar>, Mesh>,
        Space>
  {
    public:
      /// @brief Finite element space type.
      using FESType = H1<K, Math::SpatialVector<Scalar>, Mesh>;
      static constexpr ShapeFunctionSpaceType SpaceType = Space;

      /// @brief Operand type.
      using OperandType = ShapeFunction<NestedDerived, FESType, SpaceType>;

      using Parent      =
        ShapeFunctionBase<
          Jacobian<ShapeFunction<NestedDerived, FESType, SpaceType>>,
          FESType,
          SpaceType>;

      using ScalarType        = typename FormLanguage::Traits<FESType>::ScalarType;

      using RangeType         = Math::SpatialMatrix<ScalarType>;

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
        /// @brief Cached physical gradients per scalar DOF (size = nscalar).
        std::vector<SpatialVectorType> gradPhys;
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
        const auto* qf = ip.getQuadratureFormula();
        const size_t qp = qf ? ip.getIndex() : 0;

        const auto& poly = p.getPolytope();
        const auto  geom = poly.getGeometry();
        const size_t d   = poly.getDimension();
        const Index  cell = poly.getIndex();

        typename Cache::Key key;
        key.geom  = geom;
        key.dim   = d;
        key.cell  = cell;
        key.qf    = qf;
        key.qp    = qp;
        key.valid = true;

        if (qf && m_cache.key == key)
          return *this;

        m_cache.key = key;

        // scalar element for geometry-only tabulation
        const H1Element<K, ScalarType> feS(geom);
        const size_t nscalar = feS.getCount();

        if (m_cache.gradPhys.size() != nscalar)
          m_cache.gradPhys.resize(nscalar);

        for (auto& g : m_cache.gradPhys)
          if (g.size() != d) g.resize(d);

        const auto* tab = qf ? &feS.getTabulation(*qf) : nullptr;
        const auto& rc = p.getReferenceCoordinates();
        const auto  JinvT = p.getJacobianInverse().transpose();

        SpatialVectorType ref(d);

        for (size_t alpha = 0; alpha < nscalar; ++alpha)
        {
          for (size_t j = 0; j < d; ++j)
            ref(j) =
              qf
                ? tab->getGradient(qp, alpha)[j]
                : feS.getBasis(alpha).template getDerivative<1>(j)(rc);

          m_cache.gradPhys[alpha] = JinvT * ref;
        }

        return *this;
      }

      RangeType getBasis(size_t local) const
      {
        assert(m_cache.key);

        const auto& fes = this->getFiniteElementSpace();
        const std::uint8_t vdim = fes.getVectorDimension();

        const auto& p = this->getIntegrationPoint().getPoint();
        const std::uint8_t d = p.getPolytope().getDimension();

        const size_t alpha = local / vdim;
        const size_t comp  = local % vdim;

        assert(alpha < m_cache.gradPhys.size());
        assert(vdim <= RangeType::MaxSize);
        assert(d <= RangeType::MaxSize);

        RangeType J(vdim, d);
        J.setZero();

        for (std::uint8_t j = 0; j < d; ++j)
          J(comp, j) = m_cache.gradPhys[alpha](j);

        return J;
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
  Jacobian(const ShapeFunction<ShapeFunctionDerived, H1<K, Math::SpatialVector<Number>, Mesh>, Space>&)
    -> Jacobian<ShapeFunction<ShapeFunctionDerived, H1<K, Math::SpatialVector<Number>, Mesh>, Space>>;
}

#endif
