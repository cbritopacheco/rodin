/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P0_GRAD_H
#define RODIN_VARIATIONAL_P0_GRAD_H

/**
 * @file
 * @brief Gradient operator specialization for P0 (piecewise constant) functions.
 *
 * For P0 functions, the gradient is zero within each element since
 * @f$ u|_K = \text{const} @f$ implies @f$ \nabla u|_K = 0 @f$.
 */

#include "Rodin/Math/Vector.h"
#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Variational/Grad.h"
#include "Rodin/Variational/IntegrationPoint.h"
#include "Rodin/Variational/ShapeFunction.h"

#include "Rodin/Variational/Exceptions/UndeterminedTraceDomainException.h"

namespace Rodin::FormLanguage
{
  /// @brief Type traits for @c Grad over a grid function: exposes the finite element
  /// space and the operand type.
  template <class Range, class Data, class Mesh>
    requires(!FormLanguage::IsMatrixRange<Range>::Value)
  struct Traits<
    Variational::Grad<Variational::GridFunction<Variational::P0<Range, Mesh>, Data>>>
  {
      /// @brief Finite element space type.
      using FESType = Variational::P0<Range, Mesh>;
      /// @brief Operand type.
      using OperandType = Variational::GridFunction<FESType, Data>;
  };

  /// @brief Type traits for @c Grad over a shape function: exposes the finite element
  /// space, the shape function space and the operand type.
  template <class NestedDerived, class Range, class Mesh,
    Variational::ShapeFunctionSpaceType Space>
    requires(!FormLanguage::IsMatrixRange<Range>::Value)
  struct Traits<Variational::Grad<
    Variational::ShapeFunction<NestedDerived, Variational::P0<Range, Mesh>, Space>>>
  {
      /// @brief Finite element space type.
      using FESType = Variational::P0<Range, Mesh>;
      /// @brief Shape function space the expression belongs to, trial or test.
      static constexpr Variational::ShapeFunctionSpaceType SpaceType = Space;
      /// @brief Operand type.
      using OperandType = Variational::ShapeFunction<NestedDerived, FESType, SpaceType>;
  };
}

namespace Rodin::Variational
{
  /**
   * @ingroup GradSpecializations
   * @brief Gradient of a P0 GridFunction.
   *
   * Since P0 functions are piecewise constant on each element,
   * their gradient is identically zero: @f$ \nabla u|_K = 0 @f$.
   */
  template <class Range, class Data, class Mesh>
    requires(!FormLanguage::IsMatrixRange<Range>::Value)
  class Grad<GridFunction<P0<Range, Mesh>, Data>> final
    : public GradBase<GridFunction<P0<Range, Mesh>, Data>,
        Grad<GridFunction<P0<Range, Mesh>, Data>>>
  {
    public:
      /// @brief Range (evaluation value) type.
      using RangeType = Range;

      /// @brief Mesh type.
      using MeshType = Mesh;

      /// @brief Finite element space type.
      using FESType = P0<RangeType, MeshType>;

      /// @brief Operand type.
      using OperandType = GridFunction<FESType, Data>;

      /// @brief Parent class type.
      using Parent = GradBase<OperandType, Grad<OperandType>>;

      /**
       * @brief Constructs the gradient of a P0 function @f$ u @f$.
       * @param[in] u P0 GridFunction (piecewise constant)
       *
       * @note The gradient of a P0 function is zero on each element.
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

      /// @brief Interpolates at a geometric point.
      /// @param out Storage for the computed result.
      /// @param p Point at which the operation is evaluated.
      void interpolate(Math::SpatialVector<Real>& out, const Geometry::Point& p) const
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
                << ". Grad at an interface with no trace domain is undefined."
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
          const auto& fe = fes.getFiniteElement(d, i);
          const auto& rc = p.getReferenceCoordinates();
          Math::SpatialVector<Real> grad(d);
          Math::SpatialVector<Real> res(d);
          res.setZero();
          for (size_t local = 0; local < fe.getCount(); local++)
          {
            fe.getGradient(local)(grad, rc);
            res += gf.getValue(fes.getGlobalIndex({d, i}, local)) * grad;
          }
          out = p.getJacobianInverse().transpose() * res;
        }
      }

      /// @brief Creates a polymorphic copy.
      /// @returns Pointer to a newly allocated copy; the caller owns the returned object.
      Grad* copy() const noexcept override
      {
        return new Grad(*this);
      }
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for Grad of a P0 GridFunction
   */
  template <class Range, class Data, class Mesh>
  Grad(const GridFunction<P0<Range, Mesh>, Data>&) -> Grad<GridFunction<P0<Range, Mesh>, Data>>;

  /**
   * @ingroup GradSpecializations
   * @brief Gradient of a P0 ShapeFunction
   */
  template <class NestedDerived, class Scalar, class Mesh,
    ShapeFunctionSpaceType SpaceType>
    requires(!FormLanguage::IsMatrixRange<Scalar>::Value)
  class Grad<ShapeFunction<NestedDerived, P0<Scalar, Mesh>, SpaceType>> final
    : public ShapeFunctionBase<
        Grad<ShapeFunction<NestedDerived, P0<Scalar, Mesh>, SpaceType>>>
  {
    public:
      /// @brief Range (evaluation value) type.
      using RangeType = Scalar;

      /// @brief Mesh type.
      using MeshType = Mesh;

      /// @brief Finite element space type.
      using FESType = P0<RangeType, MeshType>;

      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      /// @brief Shape function space the expression belongs to, trial or test.
      static constexpr ShapeFunctionSpaceType Space = SpaceType;

      /// @brief Operand type.
      using OperandType = ShapeFunction<NestedDerived, FESType, Space>;

      /// @brief Parent class type.
      using Parent = ShapeFunctionBase<Grad<OperandType>, FESType, Space>;

      /// @brief Constructs the expression from its operand.
      /// @param u Operand expression.
      Grad(const OperandType& u)
        : Parent(u.getFiniteElementSpace()),
          m_u(u)
      {}

      /// @brief Copy constructor.
      /// @param other Object to copy from.
      Grad(const Grad& other)
        : Parent(other),
          m_u(other.m_u)
      {}

      /// @brief Move constructor.
      /// @param other Object to move from.
      Grad(Grad&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u))
      {}

      /// @brief Gets the operand function.
      /// @returns The operand function.
      constexpr
      const OperandType& getOperand() const
      {
        return m_u.get();
      }

      /// @brief Gets the operand in the shape function expression.
      /// @returns The operand in the shape function expression.
      constexpr
      const auto& getLeaf() const
      {
        return getOperand().getLeaf();
      }

      /// @brief Returns the number of local basis functions for a polytope.
      /// @param element Finite element used by the operation.
      /// @returns Number of local basis functions on the selected entity.
      constexpr
      size_t getDOFs(const Geometry::Polytope& element) const
      {
        return getOperand().getDOFs(element);
      }

      /// @brief Gets the integration point the expression is evaluated at.
      /// @returns The integration point the expression is evaluated at.
      const IntegrationPoint& getIntegrationPoint() const
      {
        assert(m_ip);
        return *m_ip;
      }

      /// @brief Sets the integration point the expression is evaluated at.
      /// @param ip Integration point at which the expression is evaluated.
      /// @returns Reference to this object after the operation.
      Grad& setIntegrationPoint(const IntegrationPoint& ip)
      {
        m_ip = &ip;
        return *this;
      }

      /// @brief Gets the basis function of a local degree of freedom.
      /// @param local Index in the local numbering.
      /// @returns Value of the selected local basis function at the evaluation point.
      auto getBasis(size_t local) const
      {
        const size_t sdim = getIntegrationPoint().getPoint().getPolytope().getMesh().getSpaceDimension();
        return Math::SpatialVector<ScalarType>::Zero(sdim);
      }

      Grad* copy() const noexcept override
      {
        return new Grad(*this);
      }

    private:
      std::reference_wrapper<const OperandType> m_u;
      const IntegrationPoint* m_ip;
  };

  /// @brief Deduction guide for @c Grad.
  template <class NestedDerived, class Range, class Mesh, ShapeFunctionSpaceType Space>
  Grad(const ShapeFunction<NestedDerived, P0<Range, Mesh>, Space>&)
    -> Grad<ShapeFunction<NestedDerived, P0<Range, Mesh>, Space>>;
}

namespace Rodin::FormLanguage
{
  /// @brief Rank-three physical gradient range for the matrix P0 solution.
  template <class Scalar, class Mesh, class Data>
  struct Traits<Variational::Grad<
    Variational::GridFunction<Variational::P0<Math::SpatialMatrix<Scalar>, Mesh>, Data>>>
  {
      /// @brief Finite element space of this family.
      using FES = Variational::P0<Math::SpatialMatrix<Scalar>, Mesh>;
      /// @brief Finite element space type.
      using FESType = FES;

      /// @brief Operand type.
      using OperandType = Variational::GridFunction<FESType, Data>;

      /// @brief Range (evaluation value) type.
      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;
      /// @brief Rank-three derivative tensor.
      using RangeType = Math::SpatialTensor<ScalarType>;
  };

  /// @brief Rank-three physical gradient range for matrix P0 shape functions.
  template <class Scalar, class Mesh, class Derived,
    Variational::ShapeFunctionSpaceType Space>
  struct Traits<Variational::Grad<Variational::ShapeFunction<Derived,
    Variational::P0<Math::SpatialMatrix<Scalar>, Mesh>, Space>>>
  {
      /// @brief Finite element space of this family.
      using FES = Variational::P0<Math::SpatialMatrix<Scalar>, Mesh>;
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
   * @brief Physical gradient of a matrix P0 solution, with derivative axis last.
   * Uses this family's scalar basis and component ordering.
   */
  template <class Scalar, class Mesh, class Data>
  class Grad<GridFunction<P0<Math::SpatialMatrix<Scalar>, Mesh>, Data>> final
    : public FunctionBase<Grad<GridFunction<P0<Math::SpatialMatrix<Scalar>, Mesh>, Data>>>
  {
    public:
      /// @brief Finite element space of this family.
      using FES = P0<Math::SpatialMatrix<Scalar>, Mesh>;
      /// @brief CRTP or finite element base class.
      using Parent = FunctionBase<Grad>;
      /// @brief Cloned or referenced expression operand type.
      using OperandType = GridFunction<FES, Data>;
      /// @brief Scalar type of matrix or tensor entries.
      using ScalarType = typename FormLanguage::Traits<FES>::ScalarType;
      /// @brief Evaluated matrix, tensor, or scalar range type.
      using RangeType = Math::SpatialTensor<ScalarType>;
      /// @brief Constructs the matrix gradient with the derivative axis last.
      /// @param operand Operand expression.
      explicit Grad(const OperandType& operand)
        : m_operand(operand)
      {}
      /// @brief Constructs the matrix gradient with the derivative axis last.
      /// @param other Object to copy from.
      Grad(const Grad& other)
        : Parent(other),
          m_operand(other.m_operand)
      {}
      /// @brief Constructs the matrix gradient with the derivative axis last.
      /// @param other Object to move from.
      Grad(Grad&& other)
        : Parent(std::move(other)),
          m_operand(other.m_operand)
      {}
      /// @brief Returns the differentiated or indexed operand.
      /// @returns The differentiated or indexed operand.
      const OperandType& getOperand() const
      {
        return m_operand.get();
      }
      /// @brief Evaluates the expression at the supplied physical or integration point.
      /// @param point Point at which the operation is evaluated.
      /// @returns Value of the expression at the supplied evaluation point.
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
      /// @param point Point at which the operation is evaluated.
      /// @returns Value of the expression at the supplied evaluation point.
      RangeType getValue(const IntegrationPoint& point) const
      {
        return getValue(point.getPoint());
      }
      /// @brief Returns the polynomial order when it is known.
      /// @param poly Mesh entity used by this operation.
      /// @returns Polynomial order on the entity, or an empty optional when no order is available.
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
   * @brief Physical gradients of matrix P0 trial or test bases.
   */
  template <class Scalar, class Mesh, class Derived, ShapeFunctionSpaceType Space>
  class Grad<ShapeFunction<Derived, P0<Math::SpatialMatrix<Scalar>, Mesh>, Space>> final
    : public ShapeFunctionBase<
        Grad<ShapeFunction<Derived, P0<Math::SpatialMatrix<Scalar>, Mesh>, Space>>,
        P0<Math::SpatialMatrix<Scalar>, Mesh>, Space>
  {
    public:
      /// @brief Finite element space of this family.
      using FES = P0<Math::SpatialMatrix<Scalar>, Mesh>;
      /// @brief CRTP or finite element base class.
      using Parent = ShapeFunctionBase<Grad, FES, Space>;
      /// @brief Cloned or referenced expression operand type.
      using OperandType = ShapeFunction<Derived, FES, Space>;
      /// @brief Scalar type of matrix or tensor entries.
      using ScalarType = typename FormLanguage::Traits<FES>::ScalarType;
      /// @brief Evaluated matrix, tensor, or scalar range type.
      using RangeType = Math::SpatialTensor<ScalarType>;
      /// @brief Constructs the matrix gradient with the derivative axis last.
      /// @param operand Operand expression.
      explicit Grad(const OperandType& operand)
        : Parent(operand.getFiniteElementSpace()),
          m_operand(operand.copy())
      {}
      /// @brief Constructs the matrix gradient with the derivative axis last.
      /// @param other Object to copy from.
      Grad(const Grad& other)
        : Parent(other),
          m_operand(other.m_operand->copy())
      {}
      /// @brief Constructs the matrix gradient with the derivative axis last.
      /// @param other Object to move from.
      Grad(Grad&& other)
        : Parent(std::move(other)),
          m_operand(std::move(other.m_operand))
      {}
      /// @brief Returns the differentiated or indexed operand.
      /// @returns The differentiated or indexed operand.
      const OperandType& getOperand() const
      {
        return *m_operand;
      }
      /// @brief Returns the leaf shape function used for assembly.
      /// @returns The leaf shape function used for assembly.
      const auto& getLeaf() const
      {
        return getOperand().getLeaf();
      }
      /// @brief Returns the local basis count for the selected polytope.
      /// @param poly Mesh entity used by this operation.
      /// @returns Number of local basis functions on the selected entity.
      size_t getDOFs(const Geometry::Polytope& poly) const
      {
        return getOperand().getDOFs(poly);
      }
      /// @brief Returns the currently bound integration point.
      /// @returns The currently bound integration point.
      const IntegrationPoint& getIntegrationPoint() const
      {
        assert(m_point);
        return *m_point;
      }
      /// @brief Binds the integration point and prepares local basis values.
      /// @param point Point at which the operation is evaluated.
      /// @returns Reference to this object after the operation.
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
      /// @param local Index in the local numbering.
      /// @returns Value of the selected local basis function at the evaluation point.
      const RangeType& getBasis(size_t local) const
      {
        return m_basis.at(local);
      }
      /// @brief Returns the polynomial order when it is known.
      /// @param poly Mesh entity used by this operation.
      /// @returns Polynomial order on the entity, or an empty optional when no order is available.
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
