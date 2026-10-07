/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P1_GRAD_H
#define RODIN_VARIATIONAL_P1_GRAD_H

/**
 * @file
 * @brief Gradient operator specialization for P1 (piecewise linear) functions.
 *
 * For P1 functions, the gradient is constant on each element:
 * @f[
 *   \nabla u|_K = \sum_{i=1}^{n+1} u_i \nabla \phi_i
 * @f]
 * where @f$ \phi_i @f$ are the P1 basis functions and @f$ n @f$ is the spatial dimension.
 */

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Variational/Grad.h"
#include "Rodin/Variational/IntegrationPoint.h"
#include "Rodin/Variational/ShapeFunction.h"
#include "Rodin/Variational/Exceptions/UndeterminedTraceDomainException.h"

#include "ForwardDecls.h"

namespace Rodin::Variational
{
  /// @addtogroup GradSpecializations

  /// @brief Base class for Grad classes.
  template <class Operand, class Derived>
  class GradBase;

  /// @brief Gradient of a P1 grid function.
  template <class Scalar, class Mesh, class Data>
    requires(!FormLanguage::IsMatrixRange<Scalar>::Value)
  class Grad<GridFunction<P1<Scalar, Mesh>, Data>> final
    : public GradBase<GridFunction<P1<Scalar, Mesh>, Data>,
        Grad<GridFunction<P1<Scalar, Mesh>, Data>>>
  {
    public:
      /// @brief Finite element space type.
      using FESType = P1<Scalar, Mesh>;

      /// @brief Range (evaluation value) type.
      using RangeType = Math::SpatialVector<Scalar>;

      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<RangeType>::ScalarType;

      /// @brief Small spatial vector value type.
      using SpatialVectorType = Math::SpatialVector<ScalarType>;

      /// @brief Operand type.
      using OperandType = GridFunction<FESType, Data>;

      /// @brief Parent class type.
      using Parent = GradBase<OperandType, Grad<OperandType>>;

      /**
       * @brief Constructs the gradient of a P1 function @f$ u @f$.
       * @param[in] u P1 GridFunction (piecewise linear)
       *
       * @note The gradient is constant on each element for P1 functions.
       */
      Grad(const OperandType& u)
        : Parent(u)
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Grad object to copy
       */
      Grad(const Grad& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Grad object to move from
       */
      Grad(Grad&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Interpolates at an integration point.
       * @param out Storage for the computed result.
       * @param ip Integration point at which the expression is evaluated.
       */
      void interpolate(SpatialVectorType& out, const IntegrationPoint& ip) const
      {
        interpolate(out, ip.getPoint());
      }

      /**
       * @brief Interpolates the gradient at a given point.
       * @param[out] out Output spatial vector for gradient
       * @param[in] p Point at which to evaluate gradient
       *
       * Computes @f$ \nabla u(p) @f$ using P1 basis function gradients.
       * Handles evaluation on faces by projecting to adjacent cells.
       */
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
                << ". Grad at an interface with no trace domain is undefined."
                << Alert::Raise;
            }
            else
            {
              for (auto& idx : inc)
              {
                const auto& tracePolytope = mesh.getPolytope(meshDim, idx);
                const Optional<Geometry::Attribute> a = tracePolytope->getAttribute();
                if (!a || !traceDomain.contains(*a)) // or traceDomain.count(*a)
                  continue;
                Math::SpatialPoint rc;
                tracePolytope->getTransformation().inverse(rc, pc);
                const Geometry::Point np(*tracePolytope, std::cref(rc), pc);
                this->interpolate(out, np);
                return;
              }
              UndeterminedTraceDomainException(
                  *this, __func__, {d, i}, traceDomain.begin(), traceDomain.end()) << Alert::Raise;
            }
            return;
          }
        }
        else // Evaluating on a cell
        {
          SpatialVectorType res(d);
          res.setZero();

          assert(d == mesh.getDimension());
          const auto& gf = this->getOperand();
          const auto& fes = gf.getFiniteElementSpace();
          const auto& fe = fes.getFiniteElement(d, i);
          const auto& rc = p.getReferenceCoordinates();
          for (size_t local = 0; local < fe.getCount(); local++)
          {
            const auto& basis = fe.getBasis(local);
            basis.getGradient()(rc);
            res += basis.getGradient()(rc) * gf[fes.getGlobalIndex({d, i}, local)];
          }
          out = p.getJacobianInverse().transpose() * res;
        }
      }

      /**
       * @brief Returns the polynomial order used on a mesh entity.
       * @param geom Reference geometry.
       * @returns Polynomial order on the entity, or an empty optional when no order is available.
       */
      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope& geom) const noexcept
      {
        const size_t k = P1Element<ScalarType>(geom.getGeometry()).getOrder();
        return (k == 0) ? 0 : (k - 1);
      }

      /**
       * @brief Creates a polymorphic copy.
       * @returns Pointer to a newly allocated copy; the caller owns the returned object.
       */
      Grad* copy() const noexcept override
      {
        return new Grad(*this);
      }
  };

  /**
   * @ingroup GradSpecializations
   * @brief Gradient of a ShapeFunction
   */
  template <class NestedDerived, class Scalar, class Mesh,
    ShapeFunctionSpaceType SpaceType>
    requires(!FormLanguage::IsMatrixRange<Scalar>::Value)
  class Grad<ShapeFunction<NestedDerived, P1<Scalar, Mesh>, SpaceType>> final
    : public ShapeFunctionBase<
        Grad<ShapeFunction<NestedDerived, P1<Scalar, Mesh>, SpaceType>>>
  {
    public:
      /// @brief Finite element space type.
      using FESType = P1<Scalar, Mesh>;
      /// @brief Shape function space the expression belongs to, trial or test.
      static constexpr ShapeFunctionSpaceType Space = SpaceType;

      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      /// @brief Small spatial vector value type.
      using SpatialVectorType = Math::SpatialVector<ScalarType>;

      /// @brief Range (evaluation value) type.
      using RangeType = Math::SpatialVector<ScalarType>;

      /// @brief Operand type.
      using OperandType = ShapeFunction<NestedDerived, FESType, Space>;

      /// @brief Parent class type.
      using Parent = ShapeFunctionBase<Grad<OperandType>, FESType, Space>;

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
          /// @brief Whether the key holds a cached entry.
              bool valid = false;

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
                return mesh == o.mesh && d == o.d && i == o.i && geom == o.geom &&
                  transOrder == o.transOrder;
              }

              /**
               * @brief Resets the key, invalidating the cached entry.
               * @param other Object to copy from.
               */
              void operator=(std::initializer_list<int> other) noexcept
              {
                valid = false;
                mesh = nullptr;
                d = 0;
                i = 0;
                geom = Geometry::Polytope::Type::Point;
                transOrder = 1;
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
            bool operator==(const QpKey& o) const noexcept
            {
              if (!valid || !o.valid)
                return false;
              return qf == o.qf && qp == o.qp;
            }

            /**
             * @brief Resets the key, invalidating the cached entry.
             * @param reset Initializer-list tag; its contents are ignored when invalidating the key.
             */
            void operator=([[maybe_unused]] std::initializer_list<int> reset) noexcept
            {
              valid = false;
              qf = nullptr;
              qp = 0;
            }
        };

        // Cached physical gradients \nabla_x φ_a (one per scalar basis function)
        /// @brief Cached gradient values.
        std::vector<SpatialVectorType> grad;

        /// @brief Key of the cached cell tabulation.
        CellKey cellKey;
        /// @brief Key of the cached quadrature-point tabulation.
        QpKey qpKey;
      };

      /**
       * @brief Constructs the expression from its operand.
       * @param u Operand expression.
       */
      Grad(const OperandType& u)
        : Parent(u.getFiniteElementSpace()),
          m_u(u),
          m_ip(nullptr)
      {}

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      Grad(const Grad& other)
        : Parent(other),
          m_u(other.m_u),
          m_ip(nullptr),
          m_cache(other.m_cache)
      {}

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      Grad(Grad&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)),
          m_ip(std::exchange(other.m_ip, nullptr)),
          m_cache(std::move(other.m_cache))
      {}

      /**
       * @brief Gets the operand function.
       * @returns The operand function.
       */
      constexpr
      const OperandType& getOperand() const
      {
        return m_u.get();
      }

      /**
       * @brief Gets the operand in the shape function expression.
       * @returns The operand in the shape function expression.
       */
      constexpr
      const auto& getLeaf() const
      {
        return getOperand().getLeaf();
      }

      /**
       * @brief Returns the number of local basis functions for a polytope.
       * @param element Finite element used by the operation.
       * @returns Number of local basis functions on the selected entity.
       */
      constexpr
      size_t getDOFs(const Geometry::Polytope& element) const
      {
        // Gradient has same number of scalar DOFs as the operand basis count.
        return getOperand().getDOFs(element);
      }

      /**
       * @brief Gets the integration point the expression is evaluated at.
       * @returns The integration point the expression is evaluated at.
       */
      constexpr
      const IntegrationPoint& getIntegrationPoint() const
      {
        assert(m_ip);
        return *m_ip;
      }

      /**
       * @brief Sets the integration point the expression is evaluated at.
       * @param ip Integration point at which the expression is evaluated.
       * @returns Reference to this object after the operation.
       */
      Grad& setIntegrationPoint(const IntegrationPoint& ip)
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

        // ---- cell key: allocate/size once per cell
        typename Cache::CellKey ckey;
        ckey.mesh = static_cast<const void*>(&mesh);
        ckey.d = d;
        ckey.i = i;
        ckey.geom = geom;
        ckey.transOrder = transOrder;
        ckey.valid = true;

        const bool cellChanged = !(m_cache.cellKey == ckey);
        if (cellChanged)
        {
          m_cache.cellKey = ckey;
          m_cache.qpKey = {}; // invalidate qp cache

          const size_t nv = Geometry::Polytope::Traits(geom).getVertexCount();
          m_cache.grad.resize(nv);
          for (auto& gvec : m_cache.grad)
          {
            gvec.resize(d);
            gvec.setZero();
          }
        }

        // ---- decide if gradients depend on quadrature point
        const bool tensorRef = (geom == Geometry::Polytope::Type::Quadrilateral) ||
          (geom == Geometry::Polytope::Type::Wedge) ||
          (geom == Geometry::Polytope::Type::Hexahedron);

        const bool needsQp = (transOrder > 1) || tensorRef;

        typename Cache::QpKey qkey;
        if (needsQp)
        {
          qkey.qf = qf;
          qkey.qp = qf ? ip.getIndex() : 0;
          qkey.valid = true;
        }
        else
        {
          // one state per cell
          qkey.qf = nullptr;
          qkey.qp = 0;
          qkey.valid = true;
        }

        const bool qpChanged = !qf || !(m_cache.qpKey == qkey);
        if (cellChanged || qpChanged)
        {
          m_cache.qpKey = qkey;

          const P1Element<ScalarType> fe(geom);
          const size_t nv = fe.getCount();

          const auto& rc =
            qf ? qf->getPoint(ip.getIndex()) : pt.getReferenceCoordinates();

          // J^{-T} at this integration point (constant for affine maps)
          const auto JinvT = pt.getJacobianInverse().transpose();

          // Compute physical gradients: \nabla_x φ_a = J^{-T} \nabla_hat φ_a
          for (size_t a = 0; a < nv; ++a)
          {
            // Reference gradient (size d). Build without using GradientFunction()
            // to avoid constructing thread_local vectors repeatedly.
            Math::SpatialVector<ScalarType> ghat(d);
            for (size_t k = 0; k < d; ++k)
              ghat(k) = fe.getBasis(a).template getDerivative<1>(k)(rc);

            m_cache.grad[a] = JinvT * ghat;
          }
        }

        return *this;
      }

      /**
       * @brief Gets the basis function of a local degree of freedom.
       * @param local Index in the local numbering.
       * @returns Value of the selected local basis function at the evaluation point.
       */
      constexpr const SpatialVectorType& getBasis(size_t local) const
      {
        assert(m_cache.cellKey);
        assert(local < m_cache.grad.size());
        return m_cache.grad[local];
      }

      /**
       * @brief Returns the polynomial order used on a mesh entity.
       * @param geom Reference geometry.
       * @returns Polynomial order on the entity, or an empty optional when no order is available.
       */
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
  /// @brief Rank-three physical gradient range for the matrix P1 solution.
  template <class Scalar, class Mesh, class Data>
  struct Traits<Variational::Grad<
    Variational::GridFunction<Variational::P1<Math::SpatialMatrix<Scalar>, Mesh>, Data>>>
  {
      /// @brief Finite element space of this family.
      using FES = Variational::P1<Math::SpatialMatrix<Scalar>, Mesh>;
      /// @brief Finite element space type.
      using FESType = FES;

      /// @brief Operand type.
      using OperandType = Variational::GridFunction<FESType, Data>;

      /// @brief Range (evaluation value) type.
      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;
      /// @brief Rank-three derivative tensor.
      using RangeType = Math::SpatialTensor<ScalarType>;
  };

  /// @brief Rank-three physical gradient range for matrix P1 shape functions.
  template <class Scalar, class Mesh, class Derived,
    Variational::ShapeFunctionSpaceType Space>
  struct Traits<Variational::Grad<Variational::ShapeFunction<Derived,
    Variational::P1<Math::SpatialMatrix<Scalar>, Mesh>, Space>>>
  {
      /// @brief Finite element space of this family.
      using FES = Variational::P1<Math::SpatialMatrix<Scalar>, Mesh>;
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
   * @brief Physical gradient of a matrix P1 solution, with derivative axis last.
   * Uses this family's scalar basis and component ordering.
   */
  template <class Scalar, class Mesh, class Data>
  class Grad<GridFunction<P1<Math::SpatialMatrix<Scalar>, Mesh>, Data>> final
    : public FunctionBase<Grad<GridFunction<P1<Math::SpatialMatrix<Scalar>, Mesh>, Data>>>
  {
    public:
      /// @brief Finite element space of this family.
      using FES = P1<Math::SpatialMatrix<Scalar>, Mesh>;
      /// @brief CRTP or finite element base class.
      using Parent = FunctionBase<Grad>;
      /// @brief Cloned or referenced expression operand type.
      using OperandType = GridFunction<FES, Data>;
      /// @brief Scalar type of matrix or tensor entries.
      using ScalarType = typename FormLanguage::Traits<FES>::ScalarType;
      /// @brief Evaluated matrix, tensor, or scalar range type.
      using RangeType = Math::SpatialTensor<ScalarType>;
      /**
       * @brief Constructs the matrix gradient with the derivative axis last.
       * @param operand Operand expression.
       */
      explicit Grad(const OperandType& operand)
        : m_operand(operand)
      {}
      /**
       * @brief Constructs the matrix gradient with the derivative axis last.
       * @param other Object to copy from.
       */
      Grad(const Grad& other)
        : Parent(other),
          m_operand(other.m_operand)
      {}
      /**
       * @brief Constructs the matrix gradient with the derivative axis last.
       * @param other Object to move from.
       */
      Grad(Grad&& other)
        : Parent(std::move(other)),
          m_operand(other.m_operand)
      {}
      /**
       * @brief Returns the differentiated or indexed operand.
       * @returns The differentiated or indexed operand.
       */
      const OperandType& getOperand() const
      {
        return m_operand.get();
      }
      /**
       * @brief Evaluates the expression at the supplied physical or integration point.
       * @param point Point at which the operation is evaluated.
       * @returns Value of the expression at the supplied evaluation point.
       */
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
      /**
       * @brief Evaluates the expression at the supplied physical or integration point.
       * @param point Point at which the operation is evaluated.
       * @returns Value of the expression at the supplied evaluation point.
       */
      RangeType getValue(const IntegrationPoint& point) const
      {
        return getValue(point.getPoint());
      }
      /**
       * @brief Returns the polynomial order when it is known.
       * @param poly Mesh entity used by this operation.
       * @returns Polynomial order on the entity, or an empty optional when no order is available.
       */
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
      /**
       * @brief Resolves the volume-side point for the selected trace domain.
       * @param point Evaluation point whose volume-side trace is resolved.
       * @returns Volume-side evaluation point selected by the trace domain.
       */
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
   * @brief Physical gradients of matrix P1 trial or test bases.
   */
  template <class Scalar, class Mesh, class Derived, ShapeFunctionSpaceType Space>
  class Grad<ShapeFunction<Derived, P1<Math::SpatialMatrix<Scalar>, Mesh>, Space>> final
    : public ShapeFunctionBase<
        Grad<ShapeFunction<Derived, P1<Math::SpatialMatrix<Scalar>, Mesh>, Space>>,
        P1<Math::SpatialMatrix<Scalar>, Mesh>, Space>
  {
    public:
      /// @brief Finite element space of this family.
      using FES = P1<Math::SpatialMatrix<Scalar>, Mesh>;
      /// @brief CRTP or finite element base class.
      using Parent = ShapeFunctionBase<Grad, FES, Space>;
      /// @brief Cloned or referenced expression operand type.
      using OperandType = ShapeFunction<Derived, FES, Space>;
      /// @brief Scalar type of matrix or tensor entries.
      using ScalarType = typename FormLanguage::Traits<FES>::ScalarType;
      /// @brief Evaluated matrix, tensor, or scalar range type.
      using RangeType = Math::SpatialTensor<ScalarType>;
      /**
       * @brief Constructs the matrix gradient with the derivative axis last.
       * @param operand Operand expression.
       */
      explicit Grad(const OperandType& operand)
        : Parent(operand.getFiniteElementSpace()),
          m_operand(operand.copy())
      {}
      /**
       * @brief Constructs the matrix gradient with the derivative axis last.
       * @param other Object to copy from.
       */
      Grad(const Grad& other)
        : Parent(other),
          m_operand(other.m_operand->copy())
      {}
      /**
       * @brief Constructs the matrix gradient with the derivative axis last.
       * @param other Object to move from.
       */
      Grad(Grad&& other)
        : Parent(std::move(other)),
          m_operand(std::move(other.m_operand))
      {}
      /**
       * @brief Returns the differentiated or indexed operand.
       * @returns The differentiated or indexed operand.
       */
      const OperandType& getOperand() const
      {
        return *m_operand;
      }
      /**
       * @brief Returns the leaf shape function used for assembly.
       * @returns The leaf shape function used for assembly.
       */
      const auto& getLeaf() const
      {
        return getOperand().getLeaf();
      }
      /**
       * @brief Returns the local basis count for the selected polytope.
       * @param poly Mesh entity used by this operation.
       * @returns Number of local basis functions on the selected entity.
       */
      size_t getDOFs(const Geometry::Polytope& poly) const
      {
        return getOperand().getDOFs(poly);
      }
      /**
       * @brief Returns the currently bound integration point.
       * @returns The currently bound integration point.
       */
      const IntegrationPoint& getIntegrationPoint() const
      {
        assert(m_point);
        return *m_point;
      }
      /**
       * @brief Binds the integration point and prepares local basis values.
       * @param point Point at which the operation is evaluated.
       * @returns Reference to this object after the operation.
       */
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
      /**
       * @brief Returns a basis value at the bound integration point.
       * @param local Index in the local numbering.
       * @returns Value of the selected local basis function at the evaluation point.
       */
      const RangeType& getBasis(size_t local) const
      {
        return m_basis.at(local);
      }
      /**
       * @brief Returns the polynomial order when it is known.
       * @param poly Mesh entity used by this operation.
       * @returns Polynomial order on the entity, or an empty optional when no order is available.
       */
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
