/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_H1_GRAD_H
#define RODIN_VARIATIONAL_H1_GRAD_H

/**
 * @file
 * @brief Gradient operator specialization for H1 (higher-order Lagrange) functions.
 *
 * For H1<K> functions, the gradient is polynomial of degree K-1 on each element:
 * @f[
 *   \nabla u|_K = \sum_{i=1}^{n} u_i \nabla \phi_i
 * @f]
 * where @f$ \phi_i @f$ are the H1<K> basis functions.
 */

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/Point.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Variational/Grad.h"
#include "Rodin/Variational/IntegrationPoint.h"
#include "Rodin/Variational/ShapeFunction.h"
#include "Rodin/Variational/Exceptions/UndeterminedTraceDomainException.h"

#include "ForwardDecls.h"

namespace Rodin::Variational
{
  /**
   * @addtogroup GradSpecializations
   */

  /**
   * @ingroup GradSpecializations
   * @brief Gradient of a GridFunction on H1<K> space
   */
  template <size_t K, class Scalar, class Mesh, class Data>
    requires(!FormLanguage::IsMatrixRange<Scalar>::Value)
  class Grad<GridFunction<H1<K, Scalar, Mesh>, Data>> final
    : public GradBase<GridFunction<H1<K, Scalar, Mesh>, Data>,
        Grad<GridFunction<H1<K, Scalar, Mesh>, Data>>>
  {
    public:
      /// @brief Finite element space type.
      using FESType = H1<K, Scalar, Mesh>;
      /// @brief Range (evaluation value) type.
      using RangeType = typename FormLanguage::Traits<FESType>::RangeType;
      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<RangeType>::ScalarType;
      /// @brief Small spatial vector value type.
      using SpatialVectorType = Math::SpatialVector<ScalarType>;
      /// @brief Operand type.
      using OperandType = GridFunction<FESType, Data>;
      /// @brief Parent class type.
      using Parent = GradBase<OperandType, Grad<OperandType>>;

      /// @brief Constructs the expression from its operand.
      Grad(const OperandType& u) : Parent(u) {}
      /// @brief Copy constructor.
      Grad(const Grad& other) : Parent(other) {}
      /// @brief Move constructor.
      Grad(Grad&& other) : Parent(std::move(other)) {}

      /// @brief Interpolates at an integration point.
      void interpolate(SpatialVectorType& out, const IntegrationPoint& ip) const
      {
        const auto& p = ip.getPoint();
        const auto& polytope = p.getPolytope();
        const size_t d = polytope.getDimension();
        const Index  i = polytope.getIndex();

        const auto& gf  = this->getOperand();
        const auto& fes = gf.getFiniteElementSpace();
        const auto& fe  = fes.getFiniteElement(d, i);
        const auto* qf = ip.getQuadratureFormula();
        assert(qf);
        const auto& tab = fe.getTabulation(*qf);
        const auto JinvT = p.getJacobianInverse().transpose();

        SpatialVectorType ref(static_cast<std::uint8_t>(d));
        ref.setZero();

        for (size_t local = 0; local < fe.getCount(); ++local)
        {
          const auto gref = tab.getGradient(ip.getIndex(), local);
          const auto uval = gf[fes.getGlobalIndex({d, i}, local)];
          for (size_t j = 0; j < d; ++j)
            ref(static_cast<std::uint8_t>(j)) += uval * gref[j];
        }

        out = JinvT * ref;
      }

      /// @brief Interpolates at a geometric point.
      void interpolate(SpatialVectorType& out, const Geometry::Point& p) const
      {
        const auto& polytope = p.getPolytope();
        const auto& d = polytope.getDimension();
        const auto& i = polytope.getIndex();
        const auto& mesh = polytope.getMesh();
        const size_t meshDim = mesh.getDimension();

        if (d == meshDim - 1) // face
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

          assert(inc.size() == 2);
          const auto& traceDomain = this->getTraceDomain();
          if (traceDomain.size() == 0)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "No trace domain provided: "
              << Alert::Notation::Predicate(true, "getTraceDomain().size() == 0")
              << ". Grad at an interface with no trace domain is undefined."
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
              *this, __func__, {d, i}, traceDomain.begin(), traceDomain.end())
            << Alert::Raise;
          return;
        }
        else // cell
        {
          out.resize(d);
          out.setZero();

          assert(d == mesh.getDimension());

          const auto& gf  = this->getOperand();
          const auto& fes = gf.getFiniteElementSpace();
          const auto& fe  = fes.getFiniteElement(d, i);
          const auto& rc  = p.getReferenceCoordinates();

          for (size_t local = 0; local < fe.getCount(); ++local)
          {
            const auto& basis = fe.getBasis(local);
            out += basis.getGradient()(rc) * gf[fes.getGlobalIndex({d, i}, local)];
          }

          out = p.getJacobianInverse().transpose() * out;
        }
      }

      /// @brief Returns the polynomial order used on a mesh entity.
      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope& geom) const noexcept
      {
        const size_t k = H1Element<K, ScalarType>(geom.getGeometry()).getOrder();
        return (k == 0) ? 0 : (k - 1);
      }

      /// @brief Creates a polymorphic copy.
      Grad* copy() const noexcept override
      {
        return new Grad(*this);
      }
  };

  /**
   * @ingroup GradSpecializations
   * @brief Gradient of a ShapeFunction on H1<K> space
   */
  template <size_t K, class NestedDerived, class Scalar, class Mesh,
    ShapeFunctionSpaceType SpaceType>
    requires(!FormLanguage::IsMatrixRange<Scalar>::Value)
  class Grad<ShapeFunction<NestedDerived, H1<K, Scalar, Mesh>, SpaceType>> final
    : public ShapeFunctionBase<
        Grad<ShapeFunction<NestedDerived, H1<K, Scalar, Mesh>, SpaceType>>>
  {
    public:
      /// @brief Finite element space type.
      using FESType = H1<K, Scalar, Mesh>;
      /// @brief Shape function space the expression belongs to, trial or test.
      static constexpr ShapeFunctionSpaceType Space = SpaceType;

      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;
      /// @brief Range (evaluation value) type.
      using RangeType  = Math::SpatialVector<ScalarType>;
      /// @brief Operand type.
      using OperandType = ShapeFunction<NestedDerived, FESType, Space>;
      /// @brief Parent class type.
      using Parent = ShapeFunctionBase<Grad<OperandType>, FESType, Space>;

      /// @brief Small spatial vector value type.
      using SpatialVectorType = Math::SpatialVector<ScalarType>;

      /// @brief Per-cell tabulation cache.
      struct Cache
      {
        /// @brief Key identifying a cached tabulation.
          struct Key
          {
          /// @brief Geometry of the cached polytope.
              Geometry::Polytope::Type geom = Geometry::Polytope::Type::Point;
          /// @brief Spatial dimension.
              size_t dim = 0;
          /// @brief Cached cell tabulation.
              Index cell = 0;

          /// @brief Quadrature formula the cached tabulation belongs to.
              const QF::QuadratureFormulaBase* qf = nullptr;
          /// @brief Index of the quadrature point.
              size_t qp = 0;

          /// @brief Whether the key holds a cached entry.
              bool valid = false;

          /// @brief Tests whether the key holds a cached entry.
              explicit operator bool() const noexcept
              {
                return valid;
              }

          /// @brief Equality comparison.
              bool operator==(const Key& o) const noexcept
              {
                if (!valid || !o.valid)
                  return false;
                return geom == o.geom && dim == o.dim && cell == o.cell && qf == o.qf &&
                  qp == o.qp;
              }

          /// @brief Resets the key, invalidating the cached entry.
              void operator=(std::initializer_list<int>) noexcept
              {
                valid = false;
                geom = Geometry::Polytope::Type::Point;
                dim = 0;
                cell = 0;
                qf = nullptr;
                qp = 0;
              }
        };

        /// @brief Cached physical gradients per DOF (size = ndof).
        std::vector<SpatialVectorType> gradPhys;
        /// @brief Key identifying the cached entry.
        Key key;
      };

      /// @brief Constructs the expression from its operand.
      Grad(const OperandType& u)
        : Parent(u.getFiniteElementSpace()),
          m_u(u),
          m_ip(nullptr)
      {}

      /// @brief Copy constructor.
      Grad(const Grad& other)
        : Parent(other),
          m_u(other.m_u),
          m_ip(nullptr),
          m_cache(other.m_cache)
      {}

      /// @brief Move constructor.
      Grad(Grad&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)),
          m_ip(std::exchange(other.m_ip, nullptr)),
          m_cache(std::move(other.m_cache))
      {}

      /// @brief Gets the operand function.
      constexpr
      const OperandType& getOperand() const
      {
        return m_u.get();
      }

      /// @brief Gets the operand in the shape function expression.
      constexpr
      const auto& getLeaf() const
      {
        return getOperand().getLeaf();
      }

      /// @brief Gets the global DOF indices for a polytope.
      constexpr
      size_t getDOFs(const Geometry::Polytope& element) const
      {
        return getOperand().getDOFs(element);
      }

      /// @brief Gets the integration point the expression is evaluated at.
      constexpr
      const IntegrationPoint& getIntegrationPoint() const
      {
        assert(m_ip);
        return *m_ip;
      }

      /// @brief Sets the integration point the expression is evaluated at.
      Grad& setIntegrationPoint(const IntegrationPoint& ip)
      {
        m_ip = &ip;

        const auto& p  = ip.getPoint();
        const auto* qf = ip.getQuadratureFormula();
        const size_t qp = qf ? ip.getIndex() : 0;

        const auto& poly = p.getPolytope();
        const size_t d   = poly.getDimension();
        const Index  cell = poly.getIndex();
        const auto geom = poly.getGeometry();

        typename Cache::Key key;
        key.geom  = geom;
        key.dim   = d;
        key.cell  = cell;
        key.qf = qf;
        key.qp    = qp;
        key.valid = true;

        const bool recompute = !qf || !(m_cache.key == key);
        if (!recompute)
          return *this;

        m_cache.key = key;

        const auto& fes = this->getFiniteElementSpace();
        const auto& fe  = fes.getFiniteElement(d, cell);
        const size_t ndof = fe.getCount();

        // Ensure storage sized once.
        if (m_cache.gradPhys.size() != ndof)
          m_cache.gradPhys.resize(ndof);

        for (auto& g : m_cache.gradPhys)
        {
          if (g.size() != d)
            g.resize(d);
        }

        // Reference gradients from tabulation when integrating, otherwise
        // directly from the basis at the supplied point.
        const auto* tab = qf ? &fe.getTabulation(*qf) : nullptr;
        const auto& rc = p.getReferenceCoordinates();
        const auto JinvT = p.getJacobianInverse().transpose();

        SpatialVectorType ref(d);

        for (size_t a = 0; a < ndof; ++a)
        {
          for (size_t ii = 0; ii < d; ++ii)
            ref(ii) = qf ? tab->getGradient(qp, a)[ii]
                         : fe.getBasis(a).template getDerivative<1>(ii)(rc);

          m_cache.gradPhys[a] = JinvT * ref;
        }

        return *this;
      }

      /// @brief Gets the basis function of a local degree of freedom.
      RangeType getBasis(size_t local) const
      {
        assert(m_cache.key);
        assert(local < m_cache.gradPhys.size());
        return m_cache.gradPhys[local];
      }

      /// @brief Returns the polynomial order used on a mesh entity.
      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope& geom) const noexcept
      {
        const auto k = getOperand().getOrder(geom);
        if (!k.has_value())
          return std::nullopt;
        return (*k == 0) ? 0 : (*k - 1);
      }

      Grad* copy() const noexcept override
      {
        return new Grad(*this);
      }

    private:
      std::reference_wrapper<const OperandType> m_u;

      const IntegrationPoint* m_ip;
      Cache m_cache;
  };
}

namespace Rodin::FormLanguage
{
  /// @brief Rank-three physical gradient range for the matrix H1 solution.
  template <size_t K, class Scalar, class Mesh, class Data>
  struct Traits<Variational::Grad<Variational::GridFunction<
    Variational::H1<K, Math::SpatialMatrix<Scalar>, Mesh>, Data>>>
  {
      /// @brief Finite element space of this family.
      using FES = Variational::H1<K, Math::SpatialMatrix<Scalar>, Mesh>;
      /// @brief Finite element space type.
      using FESType = FES;

      /// @brief Operand type.
      using OperandType = Variational::GridFunction<FESType, Data>;

      /// @brief Range (evaluation value) type.
      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;
      /// @brief Rank-three derivative tensor.
      using RangeType = Math::SpatialTensor<ScalarType>;
  };

  /// @brief Rank-three physical gradient range for matrix H1 shape functions.
  template <size_t K, class Scalar, class Mesh, class Derived,
    Variational::ShapeFunctionSpaceType Space>
  struct Traits<Variational::Grad<Variational::ShapeFunction<Derived,
    Variational::H1<K, Math::SpatialMatrix<Scalar>, Mesh>, Space>>>
  {
      /// @brief Finite element space of this family.
      using FES = Variational::H1<K, Math::SpatialMatrix<Scalar>, Mesh>;
      /// @brief Finite element space type.
      using FESType = FES;
      /// @brief Scalar type of matrix or tensor entries.
      using ScalarType = typename Traits<FES>::ScalarType;
      /// @brief Evaluated matrix, tensor, or scalar range type.
      using RangeType = Math::SpatialTensor<ScalarType>;
      /// @brief Cloned or referenced expression operand type.
      using OperandType = Variational::ShapeFunction<Derived, FES, Space>;
      /// @brief Trial or test shape-function space.
      static constexpr auto SpaceType = Space;
  };
}

namespace Rodin::Variational
{
  /**
   * @ingroup GradSpecializations
   * @brief Physical gradient of a matrix H1 solution, with derivative axis last.
   * Uses this family's scalar basis and component ordering.
   */
  template <size_t K, class Scalar, class Mesh, class Data>
  class Grad<GridFunction<H1<K, Math::SpatialMatrix<Scalar>, Mesh>, Data>> final
    : public FunctionBase<
        Grad<GridFunction<H1<K, Math::SpatialMatrix<Scalar>, Mesh>, Data>>>
  {
    public:
      /// @brief Finite element space of this family.
      using FES = H1<K, Math::SpatialMatrix<Scalar>, Mesh>;
      /// @brief CRTP or finite element base class.
      using Parent = FunctionBase<Grad>;
      /// @brief Cloned or referenced expression operand type.
      using OperandType = GridFunction<FES, Data>;
      /// @brief Scalar type of matrix or tensor entries.
      using ScalarType = typename FormLanguage::Traits<FES>::ScalarType;
      /// @brief Evaluated matrix, tensor, or scalar range type.
      using RangeType = Math::SpatialTensor<ScalarType>;
      /// @brief Constructs the matrix gradient with the derivative axis last.
      explicit Grad(const OperandType& operand)
        : m_operand(operand)
      {}
      /// @brief Constructs the matrix gradient with the derivative axis last.
      Grad(const Grad& other)
        : Parent(other),
          m_operand(other.m_operand)
      {}
      /// @brief Constructs the matrix gradient with the derivative axis last.
      Grad(Grad&& other)
        : Parent(std::move(other)),
          m_operand(other.m_operand)
      {}
      /// @brief Returns the differentiated or indexed operand.
      const OperandType& getOperand() const
      {
        return m_operand.get();
      }
      /// @brief Evaluates the expression at the supplied physical or integration point.
      RangeType getValue(const Geometry::Point& point) const
      {
        const auto& fes = getOperand().getFiniteElementSpace();
        const auto p = getDerivativePoint(point);
        const auto& poly = p.getPolytope();
        const auto& dofs = fes.getDOFs(poly.getDimension(), poly.getIndex());
        RangeType value(
          fes.getRows(), fes.getColumns(), fes.getMesh().getSpaceDimension());
        value.setZero();
        const auto& fe = fes.getFiniteElement(poly.getDimension(), poly.getIndex());
        const auto& inverse = p.getJacobianInverse();
        const size_t components = fes.getRows() * fes.getColumns();
        for (size_t a = 0; a < static_cast<size_t>(dofs.size()); ++a)
        {
          const auto basis = fe.getScalarElement().getBasis(a / components);
          const size_t row = (a % components) / fes.getColumns();
          const size_t column = a % fes.getColumns();
          for (size_t l = 0; l < poly.getDimension(); ++l)
          {
            const auto derivative =
              basis.template getDerivative<1>(l)(p.getReferenceCoordinates());
            for (size_t k = 0; k < fes.getMesh().getSpaceDimension(); ++k)
              value(row, column, k) += getOperand()[dofs[a]] * derivative * inverse(l, k);
          }
        }
        return value;
      }
      /// @brief Evaluates the expression at the supplied physical or integration point.
      RangeType getValue(const IntegrationPoint& point) const
      {
        return getValue(point.getPoint());
      }
      /// @brief Returns the polynomial order when it is known.
      Optional<size_t> getOrder(const Geometry::Polytope& poly) const noexcept
      {
        const auto order = getOperand().getOrder(poly);
        return order ? Optional<size_t>(*order ? *order - 1 : 0) : std::nullopt;
      }
      Grad* copy() const noexcept override
      {
        return new Grad(*this);
      }

    private:
      /// @brief Resolves the volume-side point for the selected trace domain.
      Geometry::Point getDerivativePoint(const Geometry::Point& point) const
      {
        const auto& mesh = getOperand().getFiniteElementSpace().getMesh();
        if (!mesh.isLocalPoint(point))
        {
          if (auto included = mesh.inclusion(point))
            return getDerivativePoint(*included);
          if (mesh.isSubMesh())
            if (auto restricted = mesh.asSubMesh().restriction(point))
              return getDerivativePoint(*restricted);
          Alert::Exception() << "Derivative point is outside the finite element mesh."
                             << Alert::Raise;
        }
        const auto& poly = point.getPolytope();
        const size_t dimension = mesh.getDimension();
        if (poly.getDimension() == dimension)
          return point;
        const auto& adjacent = mesh.getConnectivity().getIncidence(
          {poly.getDimension(), dimension}, poly.getIndex());
        for (auto index : adjacent)
        {
          const auto cell = mesh.getPolytope(dimension, index);
          const auto attribute = cell->getAttribute();
          if (adjacent.size() != 1 &&
            (!attribute || !this->getTraceDomain().contains(*attribute)))
            continue;
          Math::SpatialPoint reference;
          cell->getTransformation().inverse(reference, point.getPhysicalCoordinates());
          return Geometry::Point(*cell, reference);
        }
        Alert::Exception() << "Matrix derivative requires a determined trace domain."
                           << Alert::Raise;
        return point;
      }

      std::reference_wrapper<const OperandType> m_operand;
  };

  /**
   * @ingroup GradSpecializations
   * @brief Physical gradients of matrix H1 trial or test bases.
   */
  template <size_t K, class Scalar, class Mesh, class Derived,
    ShapeFunctionSpaceType Space>
  class Grad<ShapeFunction<Derived, H1<K, Math::SpatialMatrix<Scalar>, Mesh>, Space>>
    final
    : public ShapeFunctionBase<
        Grad<ShapeFunction<Derived, H1<K, Math::SpatialMatrix<Scalar>, Mesh>, Space>>,
        H1<K, Math::SpatialMatrix<Scalar>, Mesh>, Space>
  {
    public:
      /// @brief Finite element space of this family.
      using FES = H1<K, Math::SpatialMatrix<Scalar>, Mesh>;
      /// @brief CRTP or finite element base class.
      using Parent = ShapeFunctionBase<Grad, FES, Space>;
      /// @brief Cloned or referenced expression operand type.
      using OperandType = ShapeFunction<Derived, FES, Space>;
      /// @brief Scalar type of matrix or tensor entries.
      using ScalarType = typename FormLanguage::Traits<FES>::ScalarType;
      /// @brief Evaluated matrix, tensor, or scalar range type.
      using RangeType = Math::SpatialTensor<ScalarType>;
      /// @brief Constructs the matrix gradient with the derivative axis last.
      explicit Grad(const OperandType& operand)
        : Parent(operand.getFiniteElementSpace()),
          m_operand(operand.copy())
      {}
      /// @brief Constructs the matrix gradient with the derivative axis last.
      Grad(const Grad& other)
        : Parent(other),
          m_operand(other.m_operand->copy())
      {}
      /// @brief Constructs the matrix gradient with the derivative axis last.
      Grad(Grad&& other)
        : Parent(std::move(other)),
          m_operand(std::move(other.m_operand))
      {}
      /// @brief Returns the differentiated or indexed operand.
      const OperandType& getOperand() const
      {
        return *m_operand;
      }
      /// @brief Returns the leaf shape function used for assembly.
      const auto& getLeaf() const
      {
        return getOperand().getLeaf();
      }
      /// @brief Returns the local basis count for the selected polytope.
      size_t getDOFs(const Geometry::Polytope& poly) const
      {
        return getOperand().getDOFs(poly);
      }
      /// @brief Returns the currently bound integration point.
      const IntegrationPoint& getIntegrationPoint() const
      {
        assert(m_point);
        return *m_point;
      }
      /// @brief Binds the integration point and prepares local basis values.
      Grad& setIntegrationPoint(const IntegrationPoint& point)
      {
        m_point = &point;
        const auto& fes = this->getFiniteElementSpace();
        const auto& p = point.getPoint();
        const auto& poly = p.getPolytope();
        const size_t count =
          fes.getFiniteElement(poly.getDimension(), poly.getIndex()).getCount();
        m_basis.resize(count);
        const auto& fe = fes.getFiniteElement(poly.getDimension(), poly.getIndex());
        const auto& inverse = p.getJacobianInverse();
        const size_t components = fes.getRows() * fes.getColumns();
        for (size_t a = 0; a < count; ++a)
        {
          auto& gradient = m_basis[a];
          gradient =
            RangeType(fes.getRows(), fes.getColumns(), fes.getMesh().getSpaceDimension());
          gradient.setZero();
          const auto basis = fe.getScalarElement().getBasis(a / components);
          const size_t row = (a % components) / fes.getColumns();
          const size_t column = a % fes.getColumns();
          for (size_t l = 0; l < poly.getDimension(); ++l)
          {
            const auto derivative =
              basis.template getDerivative<1>(l)(p.getReferenceCoordinates());
            for (size_t k = 0; k < fes.getMesh().getSpaceDimension(); ++k)
              gradient(row, column, k) += derivative * inverse(l, k);
          }
        }
        return *this;
      }
      /// @brief Returns a basis value at the bound integration point.
      const RangeType& getBasis(size_t local) const
      {
        return m_basis.at(local);
      }
      /// @brief Returns the polynomial order when it is known.
      Optional<size_t> getOrder(const Geometry::Polytope& poly) const noexcept
      {
        const auto order = getOperand().getOrder(poly);
        return order ? Optional<size_t>(*order ? *order - 1 : 0) : std::nullopt;
      }
      Grad* copy() const noexcept override
      {
        return new Grad(*this);
      }

    private:
      std::unique_ptr<OperandType> m_operand;
      const IntegrationPoint* m_point = nullptr;
      std::vector<RangeType> m_basis;
  };
}

#endif
