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
 * On a cell of reference dimension @f$d@f$ embedded in @f$\mathbb{R}^s@f$,
 * a field with @f$m@f$ components has an @f$m\times s@f$ tangential Jacobian
 * @f$J(u)=\widehat J(u)J_K^{\dagger}@f$. The geometric inverse has shape
 * @f$d\times s@f$; ambient derivative columns are retained independently of
 * the reference and field dimensions. When @f$m=s=d@f$, this is the matrix.
 * @f[
 *   \mathbf{J}(\mathbf{u}) = \begin{pmatrix}
 *     \frac{\partial u_1}{\partial x_1} & \cdots & \frac{\partial u_1}{\partial x_d} \\
 *     \vdots & \ddots & \vdots \\
 *     \frac{\partial u_d}{\partial x_1} & \cdots & \frac{\partial u_d}{\partial x_d}
 *   \end{pmatrix}
 * @f]
 *
 * The physical derivative includes the current inverse geometric Jacobian;
 * a curved mapping does not generally preserve reference polynomial degree.
 *
 * @see <a href="_jacobian_8h.html">Jacobian</a>
 * @see <a href="_variational_2_h1_8h.html">H1</a>
 * @see <a href="_grad_8h.html">Grad</a>
 * @see <a href="_div_8h.html">Div</a>
 */
#ifndef RODIN_VARIATIONAL_H1_JACOBIAN_H
#define RODIN_VARIATIONAL_H1_JACOBIAN_H

#include <array>
#include <vector>

#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Variational/IntegrationPoint.h"
#include "Rodin/Variational/Jacobian.h"
#include "Rodin/Variational/Exceptions/UndeterminedTraceDomainException.h"

namespace Rodin::FormLanguage
{
  /**
   * @brief Type traits for @c Jacobian over a grid function, exposing the finite element
   * space and the operand type.
   */
  template <size_t K, class Range, class Data, class Mesh>
  struct Traits<Variational::Jacobian<
    Variational::GridFunction<Variational::H1<K, Range, Mesh>, Data>>>
  {
    /// @brief Finite element space type.
    using FESType = Variational::H1<K, Range, Mesh>;

    /// @brief Operand type.
    using OperandType = Variational::GridFunction<FESType, Data>;
  };

  /**
   * @brief Type traits for @c Jacobian over a shape function, exposing the finite element
   * space, the shape function space and the operand type.
   */
  template <size_t K, class NestedDerived, class Range, class Mesh,
    Variational::ShapeFunctionSpaceType Space>
  struct Traits<Variational::Jacobian<
    Variational::ShapeFunction<NestedDerived, Variational::H1<K, Range, Mesh>, Space>>>
  {
    /// @brief Finite element space type.
    using FESType = Variational::H1<K, Range, Mesh>;

    /// @brief Shape function space the expression belongs to, trial or test.
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
   * Computes the Jacobian matrix of an H1<K> vector-valued grid function.
   * @f[
   *   \mathbf{J}(\mathbf{u})|_K
   *   = \sum_{\alpha=1}^{N_s}\sum_{c=0}^{m-1} u_{\alpha c}\,
   *     \mathbf{e}_c\otimes
   *     \bigl((J_K^{\dagger})^T\widehat\nabla\widehat\phi_\alpha\bigr)
   * @f]
   * where @f$N_s@f$ is the scalar basis count and @f$m@f$ the field component
   * count. Reference gradients have @f$d@f$ entries; physical gradients retain
   * all @f$s@f$ ambient entries.
   *
   * The inverse geometric Jacobian is evaluated at the current sample. Curved
   * mappings need not produce polynomial physical derivatives.
   *
   * # Applications
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
    : public JacobianBase<GridFunction<H1<K, Math::SpatialVector<Scalar>, Mesh>, Data>,
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

      /// @brief Small spatial matrix value type.
      using SpatialMatrixType = Math::SpatialMatrix<ScalarType>;

      /// @brief Small spatial vector value type.
      using SpatialVectorType = Math::SpatialVector<ScalarType>;

      /**
       * @brief Constructs the expression from its operand.
       * @param u Operand expression.
       */
      Jacobian(const OperandType& u)
        : Parent(u)
      {}

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      Jacobian(const Jacobian& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      Jacobian(Jacobian&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Returns the polynomial order used on a mesh entity.
       * @param polytope Mesh entity used by this operation.
       * @returns Polynomial order on the entity, or an empty optional when no order is available.
       */
      constexpr Optional<size_t> getOrder(
        const Geometry::Polytope& polytope) const noexcept
      {
        const size_t k =
          H1Element<K, ScalarType>(polytope.getGeometry()).getDegree().value_or(2 * K);
        // A reference-linear basis can have a varying physical derivative
        // when its coordinate map is non-affine or its degree is unknown.
        if (k == 1)
        {
          const auto degree = polytope.getTransformation().getDegree();
          if (!degree || *degree > 1)
            return 1;
        }
        return (k == 0) ? 0 : (k - 1);
      }

      /**
       * @brief Creates a polymorphic copy.
       * @returns Pointer to a newly allocated copy; the caller owns the returned object.
       */
      Jacobian* copy() const noexcept override
      {
        return new Jacobian(*this);
      }

      /**
       * @brief Interpolates at an integration point.
       * @param out Storage for the computed result.
       * @param ip Integration point at which the expression is evaluated.
       */
      void interpolate(SpatialMatrixType& out, const IntegrationPoint& ip) const
      {
        const auto& p = ip.getPoint();
        if (!ip.getReferenceQuadratureFormula() ||
          p.getPolytope().getDimension() != p.getPolytope().getMesh().getDimension())
        {
          interpolate(out, p);
          return;
        }
        const auto& polytope = p.getPolytope();
        const size_t d = polytope.getDimension();
        const size_t spaceDimension = polytope.getMesh().getSpaceDimension();
        const Index i = polytope.getIndex();

        const auto& gf = this->getOperand();
        const auto& fes = gf.getFiniteElementSpace();
        const size_t vdim = fes.getVectorDimension();

        const auto feS = H1Element<K, ScalarType>(polytope.getGeometry());
        const size_t nscalar = feS.getCount();
        const auto* qf = ip.getReferenceQuadratureFormula();
        assert(qf);
        const auto& tab = feS.getTabulation(*qf);
        const auto JinvT = p.getJacobianInverse().transpose();

        SpatialVectorType ref(d);
        SpatialVectorType phys(spaceDimension);

        out.resize(vdim, spaceDimension);
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
            for (size_t j = 0; j < spaceDimension; ++j)
              out(comp, j) += uval * phys(j);
          }
        }
      }

      /**
       * @brief Interpolates at a geometric point.
       * @param out Storage for the computed result.
       * @param p Point at which the operation is evaluated.
       */
      void interpolate(SpatialMatrixType& out, const Geometry::Point& p) const
      {
        const auto& polytope = p.getPolytope();
        const auto& d = polytope.getDimension();
        const auto& i = polytope.getIndex();
        const auto& mesh = polytope.getMesh();
        const size_t meshDim = mesh.getDimension();
        const size_t spaceDimension = mesh.getSpaceDimension();

        if (d == meshDim - 1) // face
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
              *this, __func__, {d, i}, traceDomain.begin(), traceDomain.end())
              << Alert::Raise;
            return;
          }
        }

        assert(d == mesh.getDimension());

        const auto& gf = this->getOperand();
        const auto& fes = gf.getFiniteElementSpace();
        const size_t vdim = fes.getVectorDimension();

        const auto geom = polytope.getGeometry();
        const auto feS = H1Element<K, ScalarType>(geom);
        const auto& rc = p.getReferenceCoordinates();
        const auto JinvT = p.getJacobianInverse().transpose();

        SpatialVectorType ref(d);
        SpatialVectorType phys(spaceDimension);

        out.resize(vdim, spaceDimension);
        out.setZero();

        if (static_cast<size_t>(rc.size()) < d)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Reference gradient evaluation requires " << d
            << " coordinate components, got " << rc.size() << "." << Alert::Raise;
        }
        const auto accumulate = [&](size_t alpha, const auto& gradient) {
          for (size_t j = 0; j < d; ++j)
            ref(j) = gradient(j);
          phys = JinvT * ref;
          for (size_t comp = 0; comp < vdim; ++comp)
          {
            const size_t local = alpha * vdim + comp;
            const auto uval = gf[fes.getGlobalIndex({d, i}, local)];
            for (size_t j = 0; j < spaceDimension; ++j)
              out(comp, j) += uval * phys(j);
          }
        };
        if (geom == Geometry::Polytope::Type::Tetrahedron ||
          geom == Geometry::Polytope::Type::Pyramid)
        {
          const auto tabulate = [&]<size_t Count>() {
            std::array<ScalarType, 3 * Count> storage;
            feS.tabulate(typename H1Element<K, ScalarType>::Gradients{storage, 3}, rc);
            for (size_t alpha = 0; alpha < Count; ++alpha)
            {
              for (size_t j = 0; j < 3; ++j)
                ref(j) = storage[3 * alpha + j];
              accumulate(alpha, ref);
            }
          };
          if (geom == Geometry::Polytope::Type::Tetrahedron)
            tabulate.template operator()<FeketeTetrahedron<K>::Count>();
          else
            tabulate.template operator()<PyramidIndex<K>::Count>();
        }
        else
        {
          for (size_t alpha = 0; alpha < feS.getCount(); ++alpha)
          {
            const auto basis = feS.getBasis(alpha);
            for (size_t j = 0; j < d; ++j)
              ref(j) = basis.template getDerivative<1>(j)(rc);
            accumulate(alpha, ref);
          }
        }
      }
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for Jacobian of an H1 GridFunction
   * @param u Operand expression.
   */
  template <size_t K, class Range, class Data, class Mesh>
  Jacobian(const GridFunction<H1<K, Range, Mesh>, Data>& u)
    -> Jacobian<GridFunction<H1<K, Range, Mesh>, Data>>;

  /**
   * @ingroup JacobianSpecializations
   * @brief Jacobian of an H1 ShapeFunction object.
   *
   * Each basis Jacobian has one row per field component and one column per
   * ambient coordinate. Reference gradients have @f$d@f$ entries, while their
   * physical images under @f$(J_K^{\dagger})^T@f$ have @f$s@f$ entries, including
   * on embedded cells and zero-dimensional points.
   */
  template <size_t K, class NestedDerived, class Scalar, class Mesh,
    ShapeFunctionSpaceType Space>
  class Jacobian<
    ShapeFunction<NestedDerived, H1<K, Math::SpatialVector<Scalar>, Mesh>, Space>>
    final : public ShapeFunctionBase<Jacobian<ShapeFunction<NestedDerived,
                                       H1<K, Math::SpatialVector<Scalar>, Mesh>, Space>>,
              H1<K, Math::SpatialVector<Scalar>, Mesh>, Space>
  {
    public:
      /// @brief Finite element space type.
      using FESType = H1<K, Math::SpatialVector<Scalar>, Mesh>;

      /// @brief Shape function space the expression belongs to, trial or test.
      static constexpr ShapeFunctionSpaceType SpaceType = Space;

      /// @brief Operand type.
      using OperandType = ShapeFunction<NestedDerived, FESType, SpaceType>;

      /// @brief Parent class type.
      using Parent =
        ShapeFunctionBase<Jacobian<ShapeFunction<NestedDerived, FESType, SpaceType>>,
          FESType, SpaceType>;

      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      /// @brief Range (evaluation value) type.
      using RangeType = Math::SpatialMatrix<ScalarType>;

      /// @brief Small spatial matrix value type.
      using SpatialMatrixType = Math::SpatialMatrix<ScalarType>;

      /// @brief Small spatial vector value type.
      using SpatialVectorType = Math::SpatialVector<ScalarType>;

      /// @brief Physical derivatives prepared by the current integration-point bind.
      struct Cache
      {
        /// @brief Prepared physical derivatives in local basis order.
        std::vector<SpatialVectorType> gradPhys;
      };

      /**
       * @brief Constructs the expression from its operand.
       * @param u Operand expression.
       */
      Jacobian(const OperandType& u)
        : Parent(u.getFiniteElementSpace()),
          m_u(u),
          m_ip(nullptr)
      {}

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      Jacobian(const Jacobian& other)
        : Parent(other),
          m_u(other.m_u),
          m_ip(nullptr),
          m_cache(other.m_cache)
      {}

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      Jacobian(Jacobian&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)),
          m_ip(std::exchange(other.m_ip, nullptr)),
          m_cache(std::move(other.m_cache))
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
       * @param element Finite element used by the operation.
       * @returns Number of local basis functions on the selected entity.
       */
      constexpr size_t getDOFs(const Geometry::Polytope& element) const
      {
        return getOperand().getDOFs(element);
      }

      /**
       * @brief Gets the finite element space.
       * @returns The finite element space.
       */
      constexpr const auto& getFiniteElementSpace() const
      {
        return getOperand().getFiniteElementSpace();
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
       * @brief Starts a fresh cell binding and releases the previous point binding.
       * @param polytope Entity whose basis will be evaluated.
       * @returns This expression after rebinding its evaluation state.
       */
      Jacobian& setPolytope([[maybe_unused]] const Geometry::Polytope& polytope)
      {
        m_ip = nullptr;
        return *this;
      }

      /**
       * @brief Sets the integration point the expression is evaluated at.
       *
       * Every explicit bind prepares physical derivatives from the supplied
       * point. Cell and formula identities do not identify geometry state.
       * Deformation may preserve both while changing the inverse Jacobian.
       * Only reference tabulation and working storage are reused across binds.
       *
       * @param ip Integration point at which the expression is evaluated.
       * @returns Reference to this object after the operation.
       */
      Jacobian& setIntegrationPoint(const IntegrationPoint& ip)
      {
        m_ip = &ip;

        const auto& p = ip.getPoint();
        const auto* qf = ip.getReferenceQuadratureFormula();
        const size_t qp = qf ? ip.getIndex() : 0;

        const auto& poly = p.getPolytope();
        const auto geom = poly.getGeometry();
        const size_t d = poly.getDimension();

        // scalar element for geometry-only tabulation
        const H1Element<K, ScalarType> feS(geom);
        const size_t nscalar = feS.getCount();

        if (m_cache.gradPhys.size() != nscalar)
          m_cache.gradPhys.resize(nscalar);

        for (auto& g : m_cache.gradPhys)
        {
          if (g.size() != d)
            g.resize(d);
        }

        const auto* tab = qf ? &feS.getTabulation(*qf) : nullptr;
        const auto& rc = p.getReferenceCoordinates();
        const auto JinvT = p.getJacobianInverse().transpose();

        SpatialVectorType ref(d);

        for (size_t alpha = 0; alpha < nscalar; ++alpha)
        {
          for (size_t j = 0; j < d; ++j)
          {
            ref(j) = qf ? tab->getGradient(qp, alpha)[j]
                        : feS.getBasis(alpha).template getDerivative<1>(j)(rc);
          }

          m_cache.gradPhys[alpha] = JinvT * ref;
        }

        return *this;
      }

      /**
       * @brief Gets the basis function of a local degree of freedom.
       * @param local Index in the local numbering.
       * @returns Value of the selected local basis function at the evaluation point.
       */
      RangeType getBasis(size_t local) const
      {
        const auto& fes = this->getFiniteElementSpace();
        const std::uint8_t vdim = fes.getVectorDimension();

        const auto& p = this->getIntegrationPoint().getPoint();
        const std::uint8_t spaceDimension = p.getPolytope().getMesh().getSpaceDimension();

        const size_t alpha = local / vdim;
        const size_t comp = local % vdim;

        assert(alpha < m_cache.gradPhys.size());
        assert(vdim <= RangeType::MaxSize);
        assert(spaceDimension <= RangeType::MaxSize);

        RangeType J(vdim, spaceDimension);
        J.setZero();

        for (std::uint8_t j = 0; j < spaceDimension; ++j)
          J(comp, j) = m_cache.gradPhys[alpha](j);

        return J;
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
        // A reference-linear basis can have a varying physical derivative
        // when its coordinate map is non-affine or its degree is unknown.
        if (*k == 1)
        {
          const auto degree = geom.getTransformation().getDegree();
          if (!degree || *degree > 1)
            return 1;
        }
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

  /**
   * @brief Deduction guide for @c Jacobian.
   * @param u Operand expression.
   */
  template <size_t K, class ShapeFunctionDerived, class Number, class Mesh,
    ShapeFunctionSpaceType Space>
  Jacobian(
    const ShapeFunction<ShapeFunctionDerived, H1<K, Math::SpatialVector<Number>, Mesh>,
      Space>& u) -> Jacobian<ShapeFunction<ShapeFunctionDerived,
                   H1<K, Math::SpatialVector<Number>, Mesh>, Space>>;
}

#endif
