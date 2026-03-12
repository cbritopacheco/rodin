/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_H1_H1ELEMENT_H
#define RODIN_VARIATIONAL_H1_H1ELEMENT_H

/**
 * @file
 * @brief Defines the H1Element class for H1-conforming Lagrange finite elements.
 *
 * This file provides the implementation of high-order H1-conforming Lagrange
 * finite elements. These elements are characterized by:
 * - Continuous basis functions (H¹-conforming)
 * - Lagrange nodal property: @f$ \phi_i(x_j) = \delta_{ij} @f$
 * - Arbitrary polynomial degree K
 * - Support for scalar and vector-valued fields
 *
 * @see H1Element
 * @see RealH1Element
 * @see VectorH1Element
 */

#include <cstddef>
#include <array>
#include <span>
#include <vector>
#include <utility>

#include <boost/serialization/access.hpp>

#include "Rodin/Types.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Geometry/Polytope.h"
#include "Rodin/QF/ForwardDecls.h"

#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Variational/FiniteElement.h"
#include "Rodin/Variational/H1/Fekete.h"
#include "Rodin/Variational/H1/GLL.h"
#include "Rodin/Math/Traits.h"

#include "ForwardDecls.h"

/**
 * @ingroup RodinDirectives
 * @brief Indicates the maximum vector dimension a H1Element can support.
 */
#define RODIN_VARIATIONAL_H1ELEMENT_MAX_VECTOR_DIMENSION 16

/**
 * @ingroup RodinDirectives
 * @brief Tolerance for numerical computations in H1Element.
 *
 * This tolerance is used in the collapsed coordinate transformations
 * (Duffy transformation) to avoid division by zero near element boundaries.
 * Used when computing basis function gradients on triangles and tetrahedra.
 */
#define RODIN_VARIATIONAL_H1ELEMENT_TOLERANCE 1e-14

namespace Rodin::FormLanguage
{
  /**
   * @ingroup TraitsSpecializations
   */
  template <size_t K, class Range>
  struct Traits<Variational::H1Element<K, Range>>
  {
    using ScalarType = typename FormLanguage::Traits<Range>::ScalarType;
    using RangeType = Range;
  };
}

namespace Rodin::Variational
{
  /**
   * @defgroup H1ElementSpecializations H1Element Template Specializations
   * @brief Template specializations of the H1Element class.
   * @see H1Element
   */

  /**
   * @ingroup FiniteElements
   * @ingroup H1ElementSpecializations
   * @brief Continuous H1-conforming piecewise polynomial (degree k) scalar Lagrange element.
   *
   * The H1Element provides a k-th order finite element with Lagrange basis functions.
   * This is a high-order-stable H1-conforming element family.
   *
   * ## Mathematical Properties
   * - **DOF count**: Depends on geometry and polynomial degree K
   *   - Segment: K+1
   *   - Triangle: (K+1)(K+2)/2
   *   - Quadrilateral: (K+1)²
   *   - Tetrahedron: (K+1)(K+2)(K+3)/6
   *   - Wedge: (K+1)·(K+1)(K+2)/2
   * - **Basis functions**: Lagrange polynomials of degree K satisfying the Lagrange property:
   *   @f$ \phi_i(x_j) = \delta_{ij} @f$ where @f$ x_j @f$ are the Lagrange nodes
   * - **Gradient**: Polynomial of degree K-1
   * - **Partition of unity**: @f$ \sum_{i=1}^{n} \phi_i(x) = 1 @f$ for all @f$ x @f$
   * - **Continuity**: Global C⁰ continuity across element interfaces (H¹-conforming)
   *
   * ## Convergence Rates
   * For smooth solutions, H1 elements of degree k provide k-th order convergence:
   * - **L² norm**: @f$ \|u - u_h\|_{L^2} = O(h^{k+1}) @f$
   * - **H¹ norm (energy)**: @f$ \|u - u_h\|_{H^1} = O(h^k) @f$
   *
   * ## Basis Function Construction
   * - **1D (Segment)**: Classical Lagrange interpolation @f$ L_i^K(x) = \prod_{j \neq i} \frac{x - x_j}{x_i - x_j} @f$
   * - **2D/3D Simplices**: Barycentric Lagrange polynomials @f$ L_n^K(\lambda) = \prod_{m=0}^{n-1} \frac{K\lambda - m}{m+1} @f$
   * - **Tensor Products**: For quadrilaterals and wedges, tensor product of lower-dimensional bases
   *
   * ## Usage Example
   * ```cpp
   * // H1 element of degree 3 on triangle - 10 DOFs
   * RealH1Element<3> h1_tri(Polytope::Type::Triangle);
   * std::cout << h1_tri.getCount() << std::endl;  // Output: 10
   *
   * // Evaluate basis function at a point
   * Math::Vector<Real> pt{{0.25, 0.25}};
   * Real value = h1_tri.getBasis(0)(pt);
   *
   * // Get gradient of basis function
   * auto grad = h1_tri.getBasis(0).getGradient()(pt);
   * ```
   *
   * @tparam K Polynomial degree (0, 1, 2, 3, ...). Higher degrees provide better approximation
   *           for smooth functions but increase computational cost.
   * @tparam Scalar Type of scalar range (e.g., Real, Complex<Real>)
   *
   * @see H1Element<K, Math::Vector<Scalar>> for vector-valued version
   * @see P0Element for piecewise constant elements
   * @see P1Element for piecewise linear elements
   */
  template <size_t K, class Scalar>
  class H1Element final : public FiniteElementBase<H1Element<K, Scalar>>
  {
    using G = Geometry::Polytope::Type;

    public:
      static_assert(K > 0, "Polynomial degree K must be greater than 0.");

      struct Tabulation
      {
        size_t nqp  = 0;
        size_t ndof = 0;
        size_t dim  = 0;

        // qp-major storage
        // phi[(qp*ndof) + a]
        std::vector<Scalar> phi;

        // derivative in reference coordinates
        // dphi[((qp*ndof + a)*dim) + i]
        std::vector<Scalar> dphi;

        // ---------- fast path (tight loops) ----------
        const Scalar& getBasis(size_t qp, size_t a) const noexcept
        {
          return phi[qp * ndof + a];
        }

        template <size_t Order>
        const Scalar& getDerivative(size_t qp, size_t a, size_t i) const noexcept
        {
          if constexpr (Order == 0)
          {
            assert(i == 0);
            return phi[qp * ndof + a];
          }
          else if constexpr (Order == 1)
          {
            assert(i < dim);
            return dphi[(qp * ndof + a) * dim + i];
          }
          else
          {
            // Higher-order derivatives not implemented
            static const Scalar zero = Scalar(0);
            return zero;
          }
        }

        std::span<const Scalar> getGradient(size_t qp, size_t a) const noexcept
        {
          return std::span<const Scalar>(&dphi[(qp * ndof + a) * dim], dim);
        }
      };

      friend class boost::serialization::access;

      /// Parent class
      using Parent = FiniteElementBase<H1Element<K, Scalar>>;

      using ScalarType = Scalar;

      /// Type of range
      using RangeType = ScalarType;

      /**
       * @brief Evaluates linear form at a given function/vector field.
       *
       * A linear form represents the evaluation functional at a Lagrange node.
       * For a scalar function v, this returns v(x_i) where x_i is the i-th Lagrange node.
       *
       * The linear form is used in assembling the right-hand side of finite element systems.
       */
      class LinearForm
      {
        public:
          /**
           * @brief Constructs a linear form for the i-th degree of freedom.
           * @param i Index of the degree of freedom (0 to getCount()-1)
           * @param g Geometry type of the element
           */
          constexpr
          LinearForm(size_t i, Geometry::Polytope::Type g)
            : m_i(i), m_g(g)
          {}

          constexpr
          LinearForm(const LinearForm&) = default;

          /**
           * @brief Evaluates the linear form on a function/vector field.
           * @tparam T Type of the function/field to evaluate
           * @param v Function or vector field to evaluate at the node
           * @return Value of v at the i-th Lagrange node
           */
          template <class T>
          ScalarType operator()(const T& v) const
          {
            return v(H1Element<K, Scalar>::getNodes(m_g)[m_i]);
          }

        private:
          const size_t m_i;              ///< Index of the degree of freedom
          const Geometry::Polytope::Type m_g;  ///< Geometry type
      };

      /**
       * @brief Represents a Lagrange basis function of degree K.
       *
       * A basis function φ_i satisfies the Lagrange property:
       * @f$ \phi_i(x_j) = \delta_{ij} @f$ where x_j are the Lagrange nodes.
       *
       * The basis functions form a partition of unity: @f$ \sum_i \phi_i(x) = 1 @f$.
       *
       * ## Evaluation Methods
       * - `operator()(point)`: Evaluates the basis function at a spatial point
       * - `getDerivative<Order>(i)`: Gets the partial derivative of order `Order` w.r.t. coordinate i
       * - `getGradient()`: Gets the gradient function (vector of first-order partial derivatives)
       *
       * ## Implementation
       * The basis functions are constructed using:
       * - **Segments**: Classical Lagrange interpolation
       * - **Triangles/Tetrahedra**: Barycentric Lagrange polynomials
       * - **Quadrilaterals/Wedges**: Tensor product of lower-dimensional bases
       */
      class BasisFunction
      {
        public:
          using ReturnType = ScalarType;  ///< Type returned by basis function evaluation

          /**
           * @brief Represents a partial derivative of a basis function.
           *
           * Computes @f$ \frac{\partial^{Order} \phi}{\partial x_i^{Order}} @f$ where
           * φ is the basis function and x_i is a spatial coordinate.
           *
           * @tparam Order Order of differentiation (1 for first derivative, 2 for second, etc.)
           */
          template <size_t Order>
          class DerivativeFunction
          {
            public:
              /**
               * @brief Constructs a derivative function.
               * @param i Coordinate index for differentiation (0 for x, 1 for y, 2 for z)
               * @param local Local index of the basis function
               * @param g Geometry type
               */
              constexpr
              DerivativeFunction(size_t i, size_t local, Geometry::Polytope::Type g)
                : m_i(i), m_local(local), m_g(g)
              {}

              constexpr
              DerivativeFunction(const DerivativeFunction&) = default;

              /**
               * @brief Evaluates the derivative at a spatial point.
               * @param r Reference point in the element
               * @return Value of the derivative at point r
               */
              ReturnType operator()(const Math::SpatialPoint& r) const;

            private:
              const size_t m_i;      ///< Coordinate index for differentiation
              const size_t m_local;  ///< Local basis function index
              const Geometry::Polytope::Type m_g;  ///< Geometry type
          };

          /**
           * @brief Represents the gradient of a basis function.
           *
           * The gradient is the vector of first-order partial derivatives:
           * @f$ \nabla \phi = \left( \frac{\partial \phi}{\partial x}, \frac{\partial \phi}{\partial y}, \frac{\partial \phi}{\partial z} \right) @f$
           */
          class GradientFunction
          {
            public:
              using ReturnType = Math::SpatialVector<ScalarType>;  ///< Gradient vector type

              /**
               * @brief Constructs a gradient function.
               * @param local Local index of the basis function
               * @param g Geometry type
               */
              constexpr
              GradientFunction(size_t local, Geometry::Polytope::Type g)
                : m_local(local), m_g(g)
              {}

              constexpr
              GradientFunction(const GradientFunction&) = default;

              /**
               * @brief Evaluates the gradient at a spatial point.
               * @param r Reference point in the element
               * @return Gradient vector at point r
               */
              const ReturnType& operator()(const Math::SpatialPoint& r) const
              {
                static thread_local ReturnType s_out;
                const size_t dim = Geometry::Polytope::Traits(m_g).getDimension();
                s_out.resize(dim);
                for (size_t i = 0; i < dim; ++i)
                  s_out(i) = DerivativeFunction<1>(i, m_local, m_g)(r);
                return s_out;
              }

            private:
              const size_t m_local;  ///< Local basis function index
              const Geometry::Polytope::Type m_g;  ///< Geometry type
          };

          /**
           * @brief Constructs a basis function.
           * @param local Local index of the basis function (0 to getCount()-1)
           * @param g Geometry type
           */
          constexpr
          BasisFunction(size_t local, Geometry::Polytope::Type g)
            : m_local(local), m_g(g)
          {}

          constexpr
          BasisFunction(const BasisFunction&) = default;

          /**
           * @brief Evaluates the basis function at a spatial point.
           * @param r Reference point in the element (typically in [0,1]^d)
           * @return Value of basis function at point r
           */
          ReturnType operator()(const Math::SpatialPoint& r) const;

          /**
           * @brief Gets the derivative function of specified order.
           * @tparam Order Order of differentiation
           * @param i Coordinate index (0=x, 1=y, 2=z)
           * @return Derivative function object
           */
          template <size_t Order>
          constexpr
          DerivativeFunction<Order> getDerivative(size_t i) const
          {
            return DerivativeFunction<Order>(i, m_local, m_g);
          }

          /**
           * @brief Gets the gradient function.
           * @return Gradient function object
           */
          constexpr
          GradientFunction getGradient() const
          {
            return GradientFunction(m_local, m_g);
          }

        private:
          const size_t m_local;  ///< Local basis function index
          const Geometry::Polytope::Type m_g;  ///< Geometry type
      };

      const Tabulation& getTabulation(const QF::QuadratureFormulaBase& qf) const;

      /**
       * @brief Builds the high-order stable nodes for the element geometry.
       *
       * Constructs the positions of DOF nodes based on the geometry type and
       * polynomial degree. Uses:
       * - Segment: Gauss-Lobatto-Legendre (GLL) nodes on [0,1]
       * - Quadrilateral: Tensor product of 1D GLL nodes
       * - Triangle: Fekete nodes on reference triangle
       * - Tetrahedron: Fekete nodes on reference tetrahedron
       * - Wedge: Tensor product of triangle Fekete nodes with 1D GLL nodes
       */
      static
      const std::vector<Math::SpatialPoint>&
      getNodes(Geometry::Polytope::Type g)
      {
        switch (g)
        {
          case Geometry::Polytope::Type::Point:
          {
            static thread_local const std::vector<Math::SpatialPoint> s_nodes = [] {
              std::vector<Math::SpatialPoint> n;
              n.emplace_back(Math::SpatialPoint{ 0 });
              return n;
            }();
            return s_nodes;
          }

          case Geometry::Polytope::Type::Segment:
          {
            static thread_local std::vector<Math::SpatialPoint> s_nodes;
            if (s_nodes.empty())
            {
              const auto& xi = GLL01<K>::getNodes();
              s_nodes.reserve(K + 1);
              for (size_t i = 0; i <= K; ++i)
                s_nodes.emplace_back(Math::SpatialPoint{{xi[i]}});
            }
            return s_nodes;
          }

          case Geometry::Polytope::Type::Triangle:
          {
            static thread_local std::vector<Math::SpatialPoint> s_nodes;
            if (s_nodes.empty())
            {
              const auto& tri = FeketeTriangle<K>::getNodes();
              s_nodes.assign(tri.begin(), tri.end());
            }
            return s_nodes;
          }

          case Geometry::Polytope::Type::Quadrilateral:
          {
            static thread_local std::vector<Math::SpatialPoint> s_nodes;
            if (s_nodes.empty())
            {
              const auto& xi = GLL01<K>::getNodes();
              s_nodes.reserve((K + 1) * (K + 1));
              for (size_t j = 0; j <= K; ++j)
                for (size_t i = 0; i <= K; ++i)
                  s_nodes.emplace_back(Math::SpatialPoint{{xi[i], xi[j]}});
            }
            return s_nodes;
          }

          case Geometry::Polytope::Type::Tetrahedron:
          {
            static thread_local std::vector<Math::SpatialPoint> s_nodes;
            if (s_nodes.empty())
            {
              const auto& tet = FeketeTetrahedron<K>::getNodes();
              s_nodes.assign(tet.begin(), tet.end());
            }
            return s_nodes;
          }

          case Geometry::Polytope::Type::Hexahedron:
          {
            static thread_local std::vector<Math::SpatialPoint> s_nodes;
            if (s_nodes.empty())
            {
              const auto& xi = GLL01<K>::getNodes();
              s_nodes.reserve((K + 1) * (K + 1) * (K + 1));
              for (size_t k = 0; k <= K; ++k)
              {
                for (size_t j = 0; j <= K; ++j)
                {
                  for (size_t i = 0; i <= K; ++i)
                  {
                    s_nodes.emplace_back(
                      Math::SpatialPoint{{xi[i], xi[j], xi[k]}});
                  }
                }
              }
            }
            return s_nodes;
          }

          case Geometry::Polytope::Type::Wedge:
          {
            static thread_local std::vector<Math::SpatialPoint> s_nodes;
            if (s_nodes.empty())
            {
              const auto& tri = FeketeTriangle<K>::getNodes();
              const auto& z = GLL01<K>::getNodes();

              s_nodes.reserve(tri.size() * (K + 1));
              for (size_t k = 0; k <= K; ++k)
              {
                for (const auto& p : tri)
                  s_nodes.emplace_back(Math::SpatialPoint{{p.x(), p.y(), z[k]}});
              }
            }
            return s_nodes;
          }
        }

        // Should be unreachable if all enum values are handled
        assert(false && "Unsupported Polytope type.");
        static thread_local const std::vector<Math::SpatialPoint> s_empty;
        return s_empty;
      }

      /**
       * @brief Default constructor. Creates an H1 element on a Point geometry.
       */
      H1Element()
        : Parent(Geometry::Polytope::Type::Point)
      {}

      /**
       * @brief Constructs an H1 element for the specified geometry.
       * @param geometry Type of element geometry (Segment, Triangle, Quadrilateral, etc.)
       *
       * Initializes the element and builds the Lagrange nodes for the specified geometry.
       * The number and positions of nodes depend on both the geometry and polynomial degree K.
       */
      H1Element(Geometry::Polytope::Type geometry)
        : Parent(geometry)
      {}

      /**
       * @brief Copy constructor.
       * @param other Element to copy from
       */
      constexpr
      H1Element(const H1Element& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor.
       * @param other Element to move from
       */
      constexpr
      H1Element(H1Element&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Copy assignment operator.
       * @param other Element to copy from
       * @return Reference to this element
       */
      constexpr
      H1Element& operator=(const H1Element& other)
      {
        Parent::operator=(other);
        return *this;
      }

      /**
       * @brief Move assignment operator.
       * @param other Element to move from
       * @return Reference to this element
       */
      constexpr
      H1Element& operator=(H1Element&& other)
      {
        Parent::operator=(std::move(other));
        return *this;
      }

      /**
       * @brief Gets the number of degrees of freedom in the finite element.
       *
       * The DOF count depends on geometry and polynomial degree K:
       * - Segment: K+1
       * - Triangle: (K+1)(K+2)/2
       * - Quadrilateral: (K+1)²
       * - Tetrahedron: (K+1)(K+2)(K+3)/6
       * - Wedge: (K+1)·(K+1)(K+2)/2
       *
       * @return Number of degrees of freedom
       */
      constexpr
      size_t getCount() const
      {
        switch (this->getGeometry())
        {
          case Geometry::Polytope::Type::Point:
            return 1;
          case Geometry::Polytope::Type::Segment:
            return K + 1;
          case Geometry::Polytope::Type::Triangle:
            return (K + 1) * (K + 2) / 2;
          case Geometry::Polytope::Type::Quadrilateral:
            return (K + 1) * (K + 1);
          case Geometry::Polytope::Type::Tetrahedron:
            return (K + 1) * (K + 2) * (K + 3) / 6;
          case Geometry::Polytope::Type::Wedge:
            return (K + 1) * (K + 1) * (K + 2) / 2;
          case Geometry::Polytope::Type::Hexahedron:
            return (K + 1) * (K + 1) * (K + 1);
        }
        return 0;
      }

      /**
       * @brief Gets the spatial coordinates of the i-th Lagrange node.
       * @param i Node index (0 to getCount()-1)
       * @return Reference coordinates of the node (typically in [0,1]^d)
       *
       * Lagrange nodes are the interpolation points where the basis functions
       * satisfy φ_i(x_j) = δ_ij. Their positions are determined by the polynomial
       * degree and geometry type.
       */
      constexpr
      const Math::SpatialPoint& getNode(size_t i) const;

      /**
       * @brief Gets the linear form (evaluation functional) for the i-th DOF.
       * @param i DOF index (0 to getCount()-1)
       * @return Linear form object that evaluates functions at the i-th node
       */
      const LinearForm& getLinearForm(size_t i) const;

      /**
       * @brief Gets the i-th basis function.
       * @param i Basis function index (0 to getCount()-1)
       * @return Basis function object
       *
       * The basis functions satisfy the Lagrange property and form a partition of unity.
       */
      const BasisFunction& getBasis(size_t i) const;

      /**
       * @brief Gets the total maximum polynomial order of the element.
       * @return Polynomial degree K
       */
      constexpr
      size_t getOrder() const
      {
        using G = Geometry::Polytope::Type;
        switch (this->getGeometry())
        {
          case G::Point:
            // No actual polynomial variation
            return 0;

          case G::Segment:
          case G::Triangle:
          case G::Tetrahedron:
            return K;

          case G::Quadrilateral:
          case G::Wedge:
            // Tensor-product type: max total degree is 2K
            return 2 * K;

          case G::Hexahedron:
            // 3D tensor product: max total degree is 3K
            return 3 * K;
        }

        assert(false && "Unsupported geometry.");
        return 0;
      }

      /**
       * @brief Serializes the element (for boost::serialization).
       * @param ar Archive to serialize to/from
       * @param version Serialization version
       */
      template<class Archive>
      void serialize(Archive& ar, const unsigned int)
      {
        ar & boost::serialization::base_object<Parent>(*this);
      }
  };

  /**
   * @ingroup FiniteElements
   * @ingroup H1ElementSpecializations
   * @brief Continuous H1-conforming piecewise polynomial (degree k) vector Lagrange element.
   *
   * Vector-valued H1 element for approximating vector fields (e.g., displacement, velocity).
   *
   * ## Mathematical Properties
   * - **DOF count**: @f$ d \times n_s @f$ where:
   *   - @f$ d @f$ is the vector dimension (vdim)
   *   - @f$ n_s @f$ is the number of scalar DOFs for the underlying H1 element
   * - **Basis functions**: @f$ \boldsymbol{\phi}_{i,c}(x) = \phi_i(x) \mathbf{e}_c @f$ where:
   *   - @f$ \phi_i(x) @f$ is the scalar H1 basis function
   *   - @f$ \mathbf{e}_c @f$ is the unit vector in direction c
   * - **Jacobian**: @f$ J_{ij} = \frac{\partial u_i}{\partial x_j} @f$ where @f$ \mathbf{u} @f$ is the vector field
   * - **Continuity**: C⁰ continuous vector field (each component is C⁰ continuous, H¹-conforming)
   *
   * ## Component Structure
   * Each vector basis function is non-zero in only one component:
   * - Basis function `i*vdim + c` has value @f$ \phi_i(x) @f$ in component c, zero elsewhere
   * - This allows efficient assembly and clear physical interpretation
   *
   * ## Applications
   * - **Linear elasticity**: Displacement field @f$ \mathbf{u}(x) @f$
   * - **Fluid mechanics**: Velocity field @f$ \mathbf{v}(x) @f$
   * - **Electromagnetics**: Electric field @f$ \mathbf{E}(x) @f$
   * - **General vector PDEs**: Any vector-valued unknown
   *
   * ## Usage Example
   * ```cpp
   * // 2D vector H1 element of degree 2 on triangle - 6 nodes × 2 components = 12 DOFs
   * VectorH1Element<2> vec_h1(2, Polytope::Type::Triangle);
   * std::cout << vec_h1.getCount() << std::endl;  // Output: 12
   *
   * // Evaluate vector basis function
   * Math::Vector<Real> pt{{0.25, 0.25}};
   * auto vec_value = vec_h1.getBasis(0)(pt);  // Returns a 2D vector
   *
   * // Get Jacobian matrix
   * auto jac = vec_h1.getBasis(0).getJacobian()(pt);  // Returns 2×2 matrix
   * ```
   *
   * @tparam K Polynomial degree
   * @tparam Scalar Type of scalar components (e.g., Real, Complex<Real>)
   *
   * @see H1Element<K, Scalar> for scalar-valued version
   */
  template <size_t K, class Scalar>
  class H1Element<K, Math::Vector<Scalar>> final
    : public FiniteElementBase<H1Element<K, Math::Vector<Scalar>>>
  {
    using G = Geometry::Polytope::Type;

    public:
      static_assert(K > 0, "Polynomial degree K must be greater than 0.");

      friend class boost::serialization::access;

      /// Parent class
      using Parent = FiniteElementBase<H1Element<K, Math::Vector<Scalar>>>;

      using ScalarType = Scalar;

      /// Type of range
      using RangeType = Math::Vector<Scalar>;

      class LinearForm
      {
        public:
          constexpr
          LinearForm(size_t vdim, size_t local, Geometry::Polytope::Type g)
            : m_vdim(vdim), m_local(local), m_g(g)
          {}

          constexpr
          LinearForm(const LinearForm&) = default;

          constexpr
          LinearForm(LinearForm&&) = default;

          template <class T>
          ScalarType operator()(const T& v) const
          {
            static thread_local RangeType s_out;
            s_out = v(H1Element<K, ScalarType>::getNodes(m_g)[m_local / m_vdim]);
            return s_out.coeff(m_local % m_vdim);
          }

        private:
          const size_t m_vdim;
          const size_t m_local;
          const Geometry::Polytope::Type m_g;
      };

      class BasisFunction
      {
        public:
          using ReturnType = Math::Vector<ScalarType>;

          /**
           * @brief Represents a derivative function of a Pk vector element.
           * @tparam Order Order of the derivative (0 for function, 1 for first
           * derivative, etc.)
           */
          template <size_t Order>
          class DerivativeFunction
          {
            public:
              constexpr
              DerivativeFunction(size_t i, size_t j, size_t vdim, size_t local, Geometry::Polytope::Type g)
                : m_i(i), m_j(j), m_vdim(vdim), m_local(local), m_g(g)
              {}

              constexpr
              DerivativeFunction(const DerivativeFunction&) = default;

              constexpr
              void operator()(Scalar& out, const Math::SpatialPoint& r) const
              {
                out = this->operator()(r);
              }

              constexpr
              Scalar operator()(const Math::SpatialPoint& rc) const
              {
                if constexpr (Order == 0)
                {
                  if (m_i == m_local % m_vdim)
                  {
                    return H1Element<K, ScalarType>(m_g).getBasis(m_local / m_vdim)(rc);
                  }
                  else
                  {
                    return 0;
                  }
                }
                else if constexpr (Order == 1)
                {
                  if (m_i == m_local % m_vdim)
                  {
                    return H1Element<K, ScalarType>(m_g).getBasis(m_local / m_vdim).template getDerivative<1>(m_j)(rc);
                  }
                  else
                  {
                    return 0;
                  }
                }
                else
                {
                  return 0;
                }
              }

            private:
              const size_t m_i, m_j;
              const size_t m_vdim;
              const size_t m_local;
              const Geometry::Polytope::Type m_g;
          };

          class JacobianFunction
          {
            public:
              using ReturnType = Math::SpatialMatrix<ScalarType>;

              constexpr
              JacobianFunction(size_t vdim, size_t local, Geometry::Polytope::Type g)
                : m_vdim(vdim), m_local(local), m_g(g)
              {}

              constexpr
              JacobianFunction(const JacobianFunction&) = default;

              constexpr
              JacobianFunction(JacobianFunction&&) = default;

              const ReturnType& operator()(const Math::SpatialPoint& r) const
              {
                static thread_local ReturnType s_out;
                const size_t dim = Geometry::Polytope::Traits(m_g).getDimension();
                s_out.resize(m_vdim, dim);
                for (size_t i = 0; i < m_vdim; ++i)
                {
                  for (size_t j = 0; j < dim; ++j)
                    s_out(i, j) = DerivativeFunction<1>(i, j, m_vdim, m_local, m_g)(r);
                }
                return s_out;
              }

            private:
              const size_t m_vdim;
              const size_t m_local;
              const Geometry::Polytope::Type m_g;
          };

          constexpr
          BasisFunction(size_t vdim, size_t local, Geometry::Polytope::Type g)
            : m_vdim(vdim), m_local(local), m_g(g), m_jac(vdim, local, g)
          {}

          constexpr
          BasisFunction(const BasisFunction&) = default;

          constexpr
          BasisFunction(BasisFunction&&) = default;

          const ReturnType& operator()(const Math::SpatialPoint& rc) const
          {
            static thread_local ReturnType s_out;
            s_out.resize(m_vdim);
            s_out.setZero();
            s_out.coeffRef(m_local % m_vdim) = H1Element<K, ScalarType>(m_g).getBasis(m_local / m_vdim)(rc);
            return s_out;
          }

          template <size_t Order>
          constexpr
          DerivativeFunction<Order> getDerivative(size_t i, size_t j) const
          {
            return DerivativeFunction<Order>(i, j, m_vdim, m_local, m_g);
          }

          constexpr
          const JacobianFunction& getJacobian() const
          {
            return m_jac;
          }

        private:
          const size_t m_vdim;
          const size_t m_local;
          const Geometry::Polytope::Type m_g;
          const JacobianFunction m_jac;
      };

      H1Element()
        : Parent(Geometry::Polytope::Type::Point), m_vdim(0)
      {}

      constexpr
      H1Element(Geometry::Polytope::Type geometry, size_t vdim)
        : Parent(geometry), m_vdim(vdim)
      {
        const size_t count = this->getCount();
        m_lfs.reserve(count);
        m_bs.reserve(count);
        for (size_t i = 0; i < count; ++i)
        {
          m_lfs.emplace_back(vdim, i, geometry);
          m_bs.emplace_back(vdim, i, geometry);
        }
      }

      constexpr
      H1Element(const H1Element& other)
        : Parent(other), m_vdim(other.m_vdim), m_lfs(other.m_lfs), m_bs(other.m_bs)
      {}

      constexpr
      H1Element(H1Element&& other)
        : Parent(std::move(other)), m_vdim(std::move(other.m_vdim)),
          m_lfs(std::move(other.m_lfs)), m_bs(std::move(other.m_bs))
      {}

      constexpr
      H1Element& operator=(const H1Element& other)
      {
        Parent::operator=(other);
        m_vdim = other.m_vdim;
        m_lfs = other.m_lfs;
        m_bs = other.m_bs;
        return *this;
      }

      constexpr
      H1Element& operator=(H1Element&& other)
      {
        Parent::operator=(std::move(other));
        m_vdim = std::exchange(other.m_vdim, 0);
        m_lfs = std::move(other.m_lfs);
        m_bs = std::move(other.m_bs);
        return *this;
      }

      constexpr
      size_t getCount() const
      {
        switch (this->getGeometry())
        {
          case Geometry::Polytope::Type::Point:
            return m_vdim * 1;
          case Geometry::Polytope::Type::Segment:
            return m_vdim * (K + 1);
          case Geometry::Polytope::Type::Triangle:
            return m_vdim * ((K + 1) * (K + 2) / 2);
          case Geometry::Polytope::Type::Quadrilateral:
            return m_vdim * ((K + 1) * (K + 1));
          case Geometry::Polytope::Type::Tetrahedron:
            return m_vdim * ((K + 1) * (K + 2) * (K + 3) / 6);
          case Geometry::Polytope::Type::Wedge:
            return m_vdim * ((K + 1) * (K + 1) * (K + 2) / 2);
          case Geometry::Polytope::Type::Hexahedron:
            return m_vdim * ((K + 1) * (K + 1) * (K + 1));
        }
        assert(false);
        return 0;
      }

      constexpr
      const LinearForm& getLinearForm(size_t local) const
      {
        return m_lfs[local];
      }

      constexpr
      const BasisFunction& getBasis(size_t local) const
      {
        return m_bs[local];
      }

      constexpr
      const Math::SpatialPoint& getNode(size_t local) const
      {
        return H1Element<K, ScalarType>::getNodes(this->getGeometry())[local / m_vdim];
      }

      constexpr
      size_t getOrder() const
      {
        using G = Geometry::Polytope::Type;
        switch (this->getGeometry())
        {
          case G::Point:
            // No actual polynomial variation
            return 0;

          case G::Segment:
          case G::Triangle:
          case G::Tetrahedron:
            return K;

          case G::Quadrilateral:
          case G::Wedge:
            // Tensor-product type: max total degree is 2K
            return 2 * K;

          case G::Hexahedron:
            // 3D tensor product: max total degree is 3K
            return 3 * K;
        }
        assert(false);
        return 0;
      }

      template<class Archive>
      void serialize(Archive& ar, const unsigned int)
      {
        ar & boost::serialization::base_object<Parent>(*this);
        ar & m_vdim;
      }

    private:
      size_t m_vdim;

      std::vector<LinearForm> m_lfs;
      std::vector<BasisFunction> m_bs;
  };
}

#include "H1Element.hpp"

#endif

