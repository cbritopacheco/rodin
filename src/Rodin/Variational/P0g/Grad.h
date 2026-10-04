/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P0G_GRAD_H
#define RODIN_VARIATIONAL_P0G_GRAD_H

/**
 * @file
 * @brief Gradient operator specialization for P0g (global constant) functions.
 *
 * For P0g functions (globally constant), the gradient is identically zero:
 *   ∇u = 0.
 *
 * This holds for all element types and at all points (cells/faces/boundary).
 */

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Variational/Grad.h"
#include "Rodin/Variational/ShapeFunction.h"

#include "Rodin/Variational/P0g/ForwardDecls.h"

namespace Rodin::Variational
{
  template <class Operand, class Derived>
  class GradBase;

  // ---------------------------------------------------------------------------
  // Grad of a GridFunction in P0g
  // ---------------------------------------------------------------------------
  /// @brief Gradient of a P0g grid function.
  template <class Scalar, class Mesh, class Data>
    requires(!FormLanguage::IsMatrixRange<Scalar>::Value)
  class Grad<GridFunction<P0g<Scalar, Mesh>, Data>> final
    : public GradBase<GridFunction<P0g<Scalar, Mesh>, Data>,
        Grad<GridFunction<P0g<Scalar, Mesh>, Data>>>
  {
    public:
      /// @brief Finite element space type.
      using FESType = P0g<Scalar, Mesh>;
      /// @brief Operand type.
      using OperandType = GridFunction<FESType, Data>;

      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;
      /// @brief Small spatial vector value type.
      using SpatialVectorType = Math::SpatialVector<ScalarType>;

      /// @brief Parent class type.
      using Parent = GradBase<OperandType, Grad<OperandType>>;

      /// @brief Constructs the expression from its operand.
      explicit Grad(const OperandType& u)
        : Parent(u)
      {}

      /// @brief Copy constructor.
      Grad(const Grad& other)
        : Parent(other)
      {}

      /// @brief Move constructor.
      Grad(Grad&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Interpolates ∇u at point p (always zero for P0g).
       */
      void interpolate(SpatialVectorType& out, const Geometry::Point& p) const
      {
        const auto& poly = p.getPolytope();
        const size_t d = poly.getDimension(); // works for cell or face
        out.resize(d);
        out.setZero();
      }

      /**
       * @brief Polynomial order of the gradient (zero).
       *
       * P0g has order 0; its gradient is identically 0.
       * Returning 0 is consistent with "zero polynomial".
       */
      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope&) const noexcept
      {
        return 0;
      }

      /// @brief Creates a polymorphic copy.
      Grad* copy() const noexcept override
      {
        return new Grad(*this);
      }
  };

  // ---------------------------------------------------------------------------
  // Grad of a ShapeFunction in P0g
  // ---------------------------------------------------------------------------
  /// @brief Gradient of a P0g shape function.
  template <class NestedDerived, class Scalar, class Mesh,
    ShapeFunctionSpaceType SpaceType>
    requires(!FormLanguage::IsMatrixRange<Scalar>::Value)
  class Grad<ShapeFunction<NestedDerived, P0g<Scalar, Mesh>, SpaceType>> final
    : public ShapeFunctionBase<
        Grad<ShapeFunction<NestedDerived, P0g<Scalar, Mesh>, SpaceType>>,
        P0g<Scalar, Mesh>, SpaceType>
  {
    public:
      /// @brief Finite element space type.
      using FESType = P0g<Scalar, Mesh>;
      /// @brief Shape function space the expression belongs to, trial or test.
      static constexpr ShapeFunctionSpaceType Space = SpaceType;

      /// @brief Operand type.
      using OperandType = ShapeFunction<NestedDerived, FESType, Space>;

      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;
      /// @brief Small spatial vector value type.
      using SpatialVectorType = Math::SpatialVector<ScalarType>;

      /// @brief Parent class type.
      using Parent = ShapeFunctionBase<Grad<OperandType>, FESType, Space>;

      /// @brief Constructs the expression from its operand.
      explicit Grad(const OperandType& u)
        : Parent(u.getFiniteElementSpace()),
          m_u(u),
          m_ip(nullptr)
      {}

      /// @brief Copy constructor.
      Grad(const Grad& other)
        : Parent(other),
          m_u(other.m_u),
          m_ip(nullptr)
      {}

      /// @brief Move constructor.
      Grad(Grad&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)),
          m_ip(std::exchange(other.m_ip, nullptr))
      {}

      /// @brief Gets the operand function.
      constexpr
      const OperandType& getOperand() const
      {
        return m_u.get();
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
        // Keep operand aligned (even though basis is constant).
        m_u.get().setIntegrationPoint(ip);
        m_ip = &ip;

        // Dimension comes from the polytope we are integrating on.
        const auto& poly = ip.getPoint().getPolytope();
        const size_t d = poly.getDimension();

        m_zero.resize(d);
        m_zero.setZero();

        return *this;
      }

      /**
       * @brief Number of local basis functions in the gradient object.
       *
       * Same as the operand basis count: scalar P0g has 1 basis, vector P0g has vdim bases.
       */
      constexpr
      size_t getDOFs(const Geometry::Polytope& element) const
      {
        return getOperand().getDOFs(element);
      }

      /**
       * @brief Gradient of any P0g basis function is zero.
       */
      constexpr
      const SpatialVectorType& getBasis(size_t local) const
      {
        (void) local;
        assert(m_ip);
        return m_zero;
      }

      /// @brief Returns the polynomial order used on a mesh entity.
      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope&) const noexcept
      {
        return 0;
      }

      Grad* copy() const noexcept override
      {
        return new Grad(*this);
      }

    private:
      std::reference_wrapper<const OperandType> m_u;
      const IntegrationPoint* m_ip;

      // Cached zero vector sized to the current integration polytope dimension.
      SpatialVectorType m_zero;
  };
}

namespace Rodin::FormLanguage
{
  /// @brief Rank-three physical gradient range for the matrix P0g solution.
  template <class Scalar, class Mesh, class Data>
  struct Traits<Variational::Grad<
    Variational::GridFunction<Variational::P0g<Math::SpatialMatrix<Scalar>, Mesh>, Data>>>
  {
      /// @brief Finite element space of this family.
      using FES = Variational::P0g<Math::SpatialMatrix<Scalar>, Mesh>;
      /// @brief Finite element space type.
      using FESType = FES;

      /// @brief Operand type.
      using OperandType = Variational::GridFunction<FESType, Data>;

      /// @brief Range (evaluation value) type.
      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;
      /// @brief Rank-three derivative tensor.
      using RangeType = Math::SpatialTensor<ScalarType>;
  };

  /// @brief Rank-three physical gradient range for matrix P0g shape functions.
  template <class Scalar, class Mesh, class Derived,
    Variational::ShapeFunctionSpaceType Space>
  struct Traits<Variational::Grad<Variational::ShapeFunction<Derived,
    Variational::P0g<Math::SpatialMatrix<Scalar>, Mesh>, Space>>>
  {
      /// @brief Finite element space of this family.
      using FES = Variational::P0g<Math::SpatialMatrix<Scalar>, Mesh>;
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
   * @brief Physical gradient of a matrix P0g solution, with derivative axis last.
   * Uses this family's scalar basis and component ordering.
   */
  template <class Scalar, class Mesh, class Data>
  class Grad<GridFunction<P0g<Math::SpatialMatrix<Scalar>, Mesh>, Data>> final
    : public FunctionBase<
        Grad<GridFunction<P0g<Math::SpatialMatrix<Scalar>, Mesh>, Data>>>
  {
    public:
      /// @brief Finite element space of this family.
      using FES = P0g<Math::SpatialMatrix<Scalar>, Mesh>;
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
   * @brief Physical gradients of matrix P0g trial or test bases.
   */
  template <class Scalar, class Mesh, class Derived, ShapeFunctionSpaceType Space>
  class Grad<ShapeFunction<Derived, P0g<Math::SpatialMatrix<Scalar>, Mesh>, Space>> final
    : public ShapeFunctionBase<
        Grad<ShapeFunction<Derived, P0g<Math::SpatialMatrix<Scalar>, Mesh>, Space>>,
        P0g<Math::SpatialMatrix<Scalar>, Mesh>, Space>
  {
    public:
      /// @brief Finite element space of this family.
      using FES = P0g<Math::SpatialMatrix<Scalar>, Mesh>;
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
