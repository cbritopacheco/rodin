/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P0G_P0GELEMENT_H
#define RODIN_VARIATIONAL_P0G_P0GELEMENT_H

/**
 * @file
 * @brief P0g (global constant) finite element implementation.
 *
 * IMPORTANT:
 * - P0g is a *space* with 1 global DOF (scalar) or vdim global DOFs (vector),
 *   but its *element* is still the constant basis function on each cell.
 *
 * This element is therefore essentially identical to P0Element:
 * - Scalar basis:  phi(x) = 1
 * - Vector basis:  phi_i(x) = e_i (constant unit vectors)
 * - All derivatives are zero.
 *
 * Supported specializations:
 *
 * | Specialization | Description |
 * |----------------|-------------|
 * | @ref Rodin::Variational::P0gElement "P0gElement<Scalar>" | Scalar globally constant reference element. |
 * | @ref Rodin::Variational::P0gElement "P0gElement<SpatialVector<Scalar>>" | Vector globally constant reference element. |
 * | @ref Rodin::Variational::P0gElement "P0gElement<SpatialMatrix<Scalar>>" | Full rectangular matrix globally constant reference element. |
 *
 * Having a dedicated P0gElement is mostly for naming/clarity and to decouple
 * P0g from P0Element headers if you prefer.
 */

#include <cassert>
#include <utility>
#include <vector>

#include "Rodin/Types.h"
#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Math/SpatialMatrix.h"
#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Variational/FiniteElement.h"

#include "ForwardDecls.h"

namespace Rodin::FormLanguage
{
  /// @brief Type traits for @c P0gElement: exposes the scalar type and the range type.
  template <class Range>
  struct Traits<Variational::P0gElement<Range>>
  {
      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<Range>::ScalarType;
      /// @brief Range (evaluation value) type.
      using RangeType = Range;
  };
}

namespace Rodin::Variational
{
  // ----------------------------------------------------------------------------
  // Scalar P0gElement
  // ----------------------------------------------------------------------------
  template <class Scalar>
  class P0gElement final : public FiniteElementBase<P0gElement<Scalar>>
  {
    using G = Geometry::Polytope::Type;

  public:
    /// @brief Parent class type.
    using Parent = FiniteElementBase<P0gElement<Scalar>>;
    /// @brief Scalar value type.
    using ScalarType = Scalar;
    /// @brief Range (evaluation value) type.
    using RangeType = Scalar;

    /// @brief Degree-of-freedom functional of the P0g element.
    class LinearForm
    {
    public:
      /// @brief Constructs the LinearForm from the given arguments.
      /// @param g Function operand.
      constexpr explicit LinearForm(G g) : m_g(g) {}
      /// @brief Copy constructor.
      constexpr LinearForm(const LinearForm&) = default;

      /// @brief Applies the functional to a callable.
      /// @returns Value of the expression at the supplied evaluation point.
      /// @param v Object whose identifier is hashed.
      template <class T>
      constexpr ScalarType operator()(const T& v) const
      {
        return v(P0gElement<ScalarType>(m_g).getNode(0));
      }

    private:
      const G m_g;
    };

    /// @brief Basis function of the P0g element.
    class BasisFunction
    {
    public:
      /// @brief Type returned by the callable.
      using ReturnType = Scalar;

      /// @brief Derivative of a P0g basis function, identically zero.
      template <size_t Order>
      class DerivativeFunction
      {
      public:
        constexpr DerivativeFunction() = default;

        /// @brief Copy constructor.
        constexpr DerivativeFunction(const DerivativeFunction&) = default;

        /// @brief Evaluates at a point on the reference element.
        /// @returns Value of the expression at the supplied evaluation point.
        constexpr ReturnType operator()(const Math::SpatialVector<Real>&) const
        {
          return ReturnType(0);
        }
      };

      constexpr BasisFunction() = default;

      /// @brief Copy constructor.
      constexpr BasisFunction(const BasisFunction&) = default;

      /// @brief Evaluates at a point on the reference element.
      /// @returns Value of the expression at the supplied evaluation point.
      constexpr ReturnType operator()(const Math::SpatialVector<Real>&) const
      {
        return ReturnType(1);
      }

      /// @brief Gets the derivative of the basis function.
      /// @returns Reference-coordinate derivative function for the basis.
      template <size_t Order>
      constexpr DerivativeFunction<Order> getDerivative(size_t) const
      {
        return DerivativeFunction<Order>();
      }
    };

    constexpr P0gElement() = default;

    /// @brief Constructs the P0gElement from the given arguments.
    /// @param geometry Reference geometry.
    constexpr explicit P0gElement(G geometry)
      : Parent(geometry)
    {}

    /// @brief Copy constructor.
    constexpr P0gElement(const P0gElement&) = default;

    /// @brief Move constructor.
    /// @param other Object to move from.
    constexpr P0gElement(P0gElement&& other)
      : Parent(std::move(other))
    {}

    /// @brief Copy assignment.
    /// @param other Object to copy from.
    /// @returns Reference to this object after the operation.
    constexpr P0gElement& operator=(const P0gElement& other)
    {
      Parent::operator=(other);
      return *this;
    }

    constexpr ~P0gElement() override = default;

    /// @brief Gets the number of degrees of freedom of the element.
    /// @returns The number of degrees of freedom of the element.
    constexpr size_t getCount() const { return 1; }

    /// @brief Gets the node of a local degree of freedom.
    /// @param i Index of the requested entry.
    /// @returns The node of a local degree of freedom.
    const Math::SpatialVector<Real>& getNode(size_t i) const
    {
      assert(i == 0);
      switch (this->getGeometry())
      {
        case G::Point:
        {
          static const Math::SpatialVector<Real> s_node{};
          return s_node;
        }
        case G::Segment:
        {
          static const Math::SpatialVector<Real> s_node{ 0.5 };
          return s_node;
        }
        case G::Triangle:
        {
          static const Math::SpatialVector<Real> s_node{
            { Real(1) / Real(3), Real(1) / Real(3) }
          };
          return s_node;
        }
        case G::Quadrilateral:
        {
          static const Math::SpatialVector<Real> s_node{ { 0.5, 0.5 } };
          return s_node;
        }
        case G::Tetrahedron:
        {
          static const Math::SpatialVector<Real> s_node{ { 0.25, 0.25, 0.25 } };
          return s_node;
        }
        case G::Pyramid:
        {
          static const Math::SpatialVector<Real> s_node{{0.375, 0.375, 0.25}};
          return s_node;
        }
        case G::Wedge:
        {
          static const Math::SpatialVector<Real> s_node{
            { Real(1) / Real(3), Real(1) / Real(3), 0.5 }
          };
          return s_node;
        }
        case G::Hexahedron:
        {
          static const Math::SpatialVector<Real> s_node{ { 0.5, 0.5, 0.5 } };
          return s_node;
        }
      }
      assert(false);
      static const Math::SpatialVector<Real> s_null{};
      return s_null;
    }

    /// @brief Gets the degree-of-freedom functional of a local degree of freedom.
    /// @returns The degree-of-freedom functional of a local degree of freedom.
    constexpr LinearForm getLinearForm(size_t) const
    {
      return LinearForm(this->getGeometry());
    }

    /// @brief Gets the basis function of a local degree of freedom.
    /// @returns Value of the selected local basis function at the evaluation point.
    constexpr BasisFunction getBasis(size_t) const
    {
      return BasisFunction();
    }

    /// @brief Returns the polynomial order.
    /// @returns Polynomial order of the finite element.
    constexpr size_t getOrder() const { return 0; }
  };

  // ----------------------------------------------------------------------------
  // Vector P0gElement
  // ----------------------------------------------------------------------------
  /// @brief Vector-valued cellwise-constant finite element.
  template <class Scalar>
  class P0gElement<Math::SpatialVector<Scalar>> final
    : public FiniteElementBase<P0gElement<Math::SpatialVector<Scalar>>>
  {
    using G = Geometry::Polytope::Type;

  public:
    /// @brief Parent class type.
    using Parent = FiniteElementBase<P0gElement<Math::SpatialVector<Scalar>>>;
    /// @brief Scalar value type.
    using ScalarType = Scalar;
    /// @brief Range (evaluation value) type.
    using RangeType = Math::SpatialVector<Scalar>;

    /// @brief Degree-of-freedom functional of the vector-valued P0g element.
    class LinearForm
    {
    public:
      constexpr LinearForm()
        : m_vdim(0), m_local(0), m_g(G::Point)
      {}

      /// @brief Constructs the functional of a local degree of freedom.
      /// @param vdim Number of components in the value range.
      /// @param local Index in the local numbering.
      /// @param g Function operand.
      constexpr LinearForm(size_t vdim, size_t local, G g)
        : m_vdim(vdim), m_local(local), m_g(g)
      {}

      /// @brief Copy constructor.
      constexpr LinearForm(const LinearForm&) = default;

      /// @brief Applies the functional to a callable.
      /// @returns Reference to the entry at the supplied indices.
      /// @param v Object whose identifier is hashed.
      template <class T>
      ScalarType operator()(const T& v) const
      {
        const auto& xi = P0gElement<ScalarType>(m_g).getNode(m_local / m_vdim);
        const auto value = v(xi);
        return value(static_cast<std::uint8_t>(m_local % m_vdim));
      }

    private:
      const size_t m_vdim;
      const size_t m_local;
      const G m_g;
    };

    /// @brief Basis function of the vector-valued P0g element.
    class BasisFunction
    {
    public:
      /// @brief Type returned by the callable.
      using ReturnType = RangeType;

      /// @brief Derivative of a vector-valued P0g basis function, identically zero.
      template <size_t Order>
      class DerivativeFunction
      {
      public:
        /// @brief Constructs the derivative of a local basis function.
        constexpr DerivativeFunction(size_t, size_t, size_t, size_t, G) {}
        /// @brief Copy constructor.
        constexpr DerivativeFunction(const DerivativeFunction&) = default;

        /// @brief Evaluates at a point on the reference element.
        /// @returns Reference to the entry at the supplied indices.
        constexpr ScalarType operator()(const Math::SpatialVector<Real>&) const
        {
          return ScalarType(0);
        }
      };

      constexpr BasisFunction()
        : m_vdim(0), m_local(0), m_g(G::Point)
      {}

      /// @brief Constructs the basis function of a local degree of freedom.
      /// @param vdim Number of components in the value range.
      /// @param local Index in the local numbering.
      /// @param g Function operand.
      constexpr BasisFunction(size_t vdim, size_t local, G g)
        : m_vdim(vdim), m_local(local), m_g(g)
      {}

      /// @brief Copy constructor.
      constexpr BasisFunction(const BasisFunction&) = default;

      /// @brief Evaluates at a point on the reference element.
      /// @returns Reference to the entry at the supplied indices.
      const ReturnType& operator()(const Math::SpatialVector<Real>&) const
      {
        static thread_local ReturnType s_out;
        s_out = ReturnType(static_cast<std::uint8_t>(m_vdim));
        s_out.setZero();
        s_out[static_cast<std::uint8_t>(m_local % m_vdim)] = ScalarType(1);
        return s_out;
      }

      /// @brief Gets the derivative of the basis function.
      /// @param i Index of the requested entry.
      /// @param j Index of the second coordinate.
      /// @returns Reference-coordinate derivative function for the basis.
      template <size_t Order>
      constexpr DerivativeFunction<Order> getDerivative(size_t i, size_t j) const
      {
        return DerivativeFunction<Order>(i, j, m_vdim, m_local, m_g);
      }

    private:
      const size_t m_vdim;
      const size_t m_local;
      const G m_g;
    };

    P0gElement()
      : Parent(G::Point), m_vdim(0)
    {}

    /// Backward-compatible: vdim defaults to spatial dimension of geometry
    /// @param geometry Reference geometry.
    constexpr explicit P0gElement(G geometry)
      : P0gElement(geometry, Geometry::Polytope::Traits(geometry).getDimension())
    {}

    /// @brief Constructs the P0gElement from the given arguments.
    /// @param geometry Reference geometry.
    /// @param vdim Number of components in the value range.
    constexpr P0gElement(G geometry, size_t vdim)
      : Parent(geometry), m_vdim(vdim)
    {
      const size_t count = getCount();
      m_lfs.reserve(count);
      m_bs.reserve(count);
      for (size_t i = 0; i < count; ++i)
      {
        m_lfs.emplace_back(vdim, i, geometry);
        m_bs.emplace_back(vdim, i, geometry);
      }
    }

    /// @brief Copy constructor.
    /// @param other Object to copy from.
    constexpr P0gElement(const P0gElement& other)
      : Parent(other)
      , m_vdim(other.m_vdim)
      , m_lfs(other.m_lfs)
      , m_bs(other.m_bs)
    {}

    /// @brief Move constructor.
    /// @param other Object to move from.
    constexpr P0gElement(P0gElement&& other)
      : Parent(std::move(other))
      , m_vdim(std::exchange(other.m_vdim, 0))
      , m_lfs(std::move(other.m_lfs))
      , m_bs(std::move(other.m_bs))
    {}

    constexpr ~P0gElement() override = default;

    /// @brief Copy assignment.
    /// @param other Object to copy from.
    /// @returns Reference to this object after the operation.
    constexpr P0gElement& operator=(const P0gElement& other)
    {
      Parent::operator=(other);
      m_vdim = other.m_vdim;
      m_lfs  = other.m_lfs;
      m_bs   = other.m_bs;
      return *this;
    }

    /// @brief Move assignment.
    /// @param other Object to move from.
    /// @returns Reference to this object after the operation.
    constexpr P0gElement& operator=(P0gElement&& other)
    {
      Parent::operator=(std::move(other));
      m_vdim = std::exchange(other.m_vdim, 0);
      m_lfs  = std::move(other.m_lfs);
      m_bs   = std::move(other.m_bs);
      return *this;
    }

    /// @brief Gets the number of degrees of freedom of the element.
    /// @returns The number of degrees of freedom of the element.
    constexpr size_t getCount() const
    {
      return m_vdim;
    }

    /// @brief Gets the degree-of-freedom functional of a local degree of freedom.
    /// @param local Index in the local numbering.
    /// @returns The degree-of-freedom functional of a local degree of freedom.
    constexpr auto getLinearForm(size_t local) const
    {
      return m_lfs.at(local);
    }

    /// @brief Gets the basis function of a local degree of freedom.
    /// @param local Index in the local numbering.
    /// @returns Value of the selected local basis function at the evaluation point.
    constexpr BasisFunction getBasis(size_t local) const
    {
      return m_bs.at(local);
    }

    /// @brief Gets the node of a local degree of freedom.
    /// @param local Index in the local numbering.
    /// @returns The node of a local degree of freedom.
    constexpr const Math::SpatialVector<Real>& getNode(size_t local) const
    {
      // All components share the same barycentric node
      return P0gElement<ScalarType>(this->getGeometry()).getNode(local / m_vdim);
    }

    /// @brief Evaluates the integrand into the output argument.
    /// @param out Storage for the computed result.
    /// @param coefficient Coefficient multiplying the expression.
    template <class Coefficient>
    constexpr void evaluate(
      RangeType& out, Coefficient&& coefficient, const Math::SpatialPoint&) const
    {
      out.resize(m_vdim);
      for (size_t component = 0; component < m_vdim; ++component)
        out(component) = coefficient(component);
    }

    /// @brief Returns the polynomial order.
    /// @returns Polynomial order of the finite element.
    constexpr size_t getOrder() const { return 0; }

  private:
    size_t m_vdim;
    std::vector<LinearForm>    m_lfs;
    std::vector<BasisFunction> m_bs;
  };
}

namespace Rodin::Variational
{
  /// @brief Matrix-valued reference element built from scalar nodal functionals.
  template <class Scalar>
  class P0gElement<Math::SpatialMatrix<Scalar>> final
    : public FiniteElementBase<P0gElement<Math::SpatialMatrix<Scalar>>>
  {
    public:
      /// @brief Scalar reference element for this family.
      using ScalarElement = P0gElement<Scalar>;
      /// @brief Scalar type of matrix or tensor entries.
      using ScalarType = typename FormLanguage::Traits<ScalarElement>::ScalarType;
      /// @brief Evaluated matrix, tensor, or scalar range type.
      using RangeType = Math::SpatialMatrix<ScalarType>;
      /// @brief CRTP or finite element base class.
      using Parent = FiniteElementBase<P0gElement<Math::SpatialMatrix<Scalar>>>;
      /// @brief Basis type of the underlying scalar element.
      using ScalarBasis = decltype(std::declval<ScalarElement>().getBasis(0));
      /// @brief Nodal functional type of the underlying scalar element.
      using ScalarLinearForm = decltype(std::declval<ScalarElement>().getLinearForm(0));

      /// @brief Replicates a scalar reference element for a rectangular matrix range.
      P0gElement(Geometry::Polytope::Type geometry, size_t rows, size_t cols)
        : P0gElement(ScalarElement(geometry), rows, cols)
      {}

      /// @brief Replicates a scalar reference element for a rectangular matrix range.
      P0gElement(const ScalarElement& scalar, size_t rows, size_t cols)
        : Parent(scalar.getGeometry()),
          m_scalar(scalar),
          m_rows(rows),
          m_cols(cols)
      {
        if (rows == 0 || cols == 0 || rows > RODIN_MAXIMAL_SPACE_DIMENSION ||
          cols > RODIN_MAXIMAL_SPACE_DIMENSION)
          Alert::Exception() << "SpatialMatrix ranges require 1 to 3 rows and columns."
                             << Alert::Raise;
      }

      /// @brief Copies the element and its component extents.
      P0gElement(const P0gElement&) = default;
      /// @brief Moves the element and its component extents.
      P0gElement(P0gElement&&) = default;
      /// @brief Copies the element and its component extents.
      P0gElement& operator=(const P0gElement&) = default;
      /// @brief Moves the element and its component extents.
      P0gElement& operator=(P0gElement&&) = default;

      /// @brief Matrix basis with one nonzero entry.
      class BasisFunction
      {
        public:
          /// @brief Selects a matrix unit multiplied by a scalar basis function.
          /// @param basis Basis used by the operation.
          /// @param rows Number of rows.
          /// @param cols Number of columns.
          /// @param component Component of the value range.
          BasisFunction(ScalarBasis basis, size_t rows, size_t cols, size_t component)
            : m_basis(std::move(basis)),
              m_rows(rows),
              m_cols(cols),
              m_component(component)
          {}

          /// @brief Evaluates the selected matrix basis or its component nodal functional.
          /// @param point Point at which the operation is evaluated.
          /// @returns Reference to the entry at the supplied indices.
          RangeType operator()(const Math::SpatialPoint& point) const
          {
            RangeType value(m_rows, m_cols);
            value.setZero();
            value(m_component / m_cols, m_component % m_cols) = m_basis(point);
            return value;
          }

          /// @brief Returns a reference-coordinate derivative of the basis.
          /// @param direction Direction in which the derivative is evaluated.
          /// @returns Reference-coordinate derivative function for the basis.
          template <size_t Order>
          auto getDerivative(size_t direction) const
          {
            return [derivative = m_basis.template getDerivative<Order>(direction),
                     rows = m_rows, cols = m_cols, component = m_component](
                     const Math::SpatialPoint& point) -> RangeType {
              RangeType value(rows, cols);
              value.setZero();
              value(component / cols, component % cols) = derivative(point);
              return value;
            };
          }

        private:
          ScalarBasis m_basis;
          size_t m_rows, m_cols, m_component;
      };

      /// @brief Matrix-range finite element or expression specialization.
      class LinearForm
      {
        public:
          /// @brief Selects a matrix entry for the scalar nodal functional.
          /// @param rows Number of rows.
          /// @param cols Number of columns.
          /// @param component Component of the value range.
          /// @param form Scalar linear form used to construct the component form.
          LinearForm(ScalarLinearForm form, size_t rows, size_t cols, size_t component)
            : m_form(std::move(form)),
              m_rows(rows),
              m_cols(cols),
              m_component(component)
          {}

          /// @brief Evaluates the selected matrix basis or its component nodal functional.
          /// @param function Function to evaluate.
          /// @returns Reference to the entry at the supplied indices.
          template <class Callable>
          ScalarType operator()(const Callable& function) const
          {
            return m_form([&](const Math::SpatialPoint& point) -> ScalarType {
              const auto value = function(point);
              if (value.rows() != m_rows || value.cols() != m_cols)
                Alert::Exception()
                  << "Matrix value does not match the finite element range."
                  << Alert::Raise;
              return value(m_component / m_cols, m_component % m_cols);
            });
          }

        private:
          ScalarLinearForm m_form;
          size_t m_rows, m_cols, m_component;
      };

      /// @brief Contracts each matrix entry with the scalar reference element.
      template <class Coefficient>
      void evaluate(
        RangeType& out, Coefficient&& coefficient, const Math::SpatialPoint& point) const
      {
        out.resize(m_rows, m_cols);
        const size_t components = m_rows * m_cols;
        for (size_t c = 0; c < components; ++c)
          m_scalar.evaluate(
            out(c / m_cols, c % m_cols),
            [&](size_t a) { return coefficient(a * components + c); }, point);
      }

      /// @brief Returns the number of local matrix basis functions.
      size_t getCount() const
      {
        return m_scalar.getCount() * m_rows * m_cols;
      }
      /// @brief Returns the polynomial order when it is known.
      size_t getOrder() const
      {
        return m_scalar.getOrder();
      }
      /// @brief Returns the number of matrix rows.
      size_t getRows() const
      {
        return m_rows;
      }
      /// @brief Returns the number of matrix columns.
      size_t getColumns() const
      {
        return m_cols;
      }
      /// @brief Returns the scalar element whose basis is replicated for matrix entries.
      const ScalarElement& getScalarElement() const
      {
        return m_scalar;
      }

      /// @brief Returns the scalar reference node for the selected component DOF.
      decltype(auto) getNode(size_t local) const
      {
        assert(local < getCount());
        return m_scalar.getNode(local / (m_rows * m_cols));
      }

      /// @brief Returns a basis value at the bound integration point.
      BasisFunction getBasis(size_t local) const
      {
        assert(local < getCount());
        return {m_scalar.getBasis(local / (m_rows * m_cols)), m_rows, m_cols,
          local % (m_rows * m_cols)};
      }

      /// @brief Returns the scalar nodal functional applied to the selected matrix entry.
      LinearForm getLinearForm(size_t local) const
      {
        assert(local < getCount());
        return {m_scalar.getLinearForm(local / (m_rows * m_cols)), m_rows, m_cols,
          local % (m_rows * m_cols)};
      }

    private:
      ScalarElement m_scalar;
      size_t m_rows, m_cols;
  };
}

#endif
