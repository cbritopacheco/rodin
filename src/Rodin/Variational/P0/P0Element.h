/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P0_P0ELEMENT_H
#define RODIN_VARIATIONAL_P0_P0ELEMENT_H

/**
 * @file
 * @brief P0 (piecewise constant) finite element implementation.
 *
 * This file provides the P0Element class template for piecewise constant
 * finite elements. P0 elements have:
 * - One degree of freedom per element (at the barycenter)
 * - Constant basis function: @f$ \phi(x) = 1 @f$
 * - Zero gradient: @f$ \nabla \phi = 0 @f$
 *
 * ## Supported Specializations
 *
 * **Scalar P0Element<Scalar>:**
 * - Works with Real and Complex scalar types
 * - Single DOF at element barycenter
 * - Constant value throughout element
 *
 * **Vector P0Element<Math::Vector<Scalar>>:**
 * - Component-wise constant vector fields
 * - @f$ d @f$ DOFs per element (one per vector component)
 * - Basis functions: @f$ \boldsymbol{\phi}_i = \mathbf{e}_j @f$ where @f$ j = i \mod d @f$
 *
 * ## Usage Example
 * @code{.cpp}
 * using namespace Rodin::Variational;
 *
 * // Scalar P0 element (real or complex)
 * RealP0Element p0_real(Polytope::Type::Triangle);
 * ComplexP0Element p0_complex(Polytope::Type::Triangle);
 *
 * // Vector P0 element (2D)
 * VectorP0Element p0_vec(2, Polytope::Type::Triangle);  // 2 DOFs
 * @endcode
 *
 * P0 elements are commonly used in discontinuous Galerkin (DG) methods,
 * mixed finite element formulations, and as test spaces for finite volume
 * schemes.
 */

#include "Rodin/Types.h"

#include "Rodin/Math/Traits.h"

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/Connectivity.h"

#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Variational/FiniteElement.h"
#include "Rodin/Variational/FiniteElementSpace.h"

#include "ForwardDecls.h"

namespace Rodin::FormLanguage
{
  /**
   * @ingroup TraitsSpecializations
   */
  template <class Range>
  struct Traits<Variational::P0Element<Range>>
  {
    using ScalarType = typename FormLanguage::Traits<Range>::ScalarType;
    using RangeType = Range;
  };
}

namespace Rodin::Variational
{
  /**
   * @defgroup P0ElementSpecializations P0Element Template Specializations
   * @brief Template specializations of the P0Element class.
   * @see P0Element
   */

  /**
   * @ingroup FiniteElements
   * @ingroup P0ElementSpecializations
   * @brief Piecewise constant (degree 0) scalar Lagrange element.
   *
   * The P0Element represents a piecewise constant finite element with:
   * - **DOF count**: 1 per element (located at barycenter)
   * - **Basis function**: @f$ \phi(x) = 1 @f$ for all @f$ x @f$ in element
   * - **Derivatives**: @f$ \nabla \phi = 0 @f$ (constant function has zero gradient)
   *
   * P0 elements are discontinuous across element interfaces, making them
   * suitable for DG methods, flux computations, and element-wise constant
   * approximations.
   *
   * @tparam Scalar Type of scalar range (e.g., Real, Complex)
   */
  template <class Scalar>
  class P0Element final : public FiniteElementBase<P0Element<Scalar>>
  {
    using G = Geometry::Polytope::Type;

    public:
      /// Parent class
      using Parent = FiniteElementBase<P0Element<Scalar>>;

      using ScalarType = Scalar;

      /// Type of range
      using RangeType = Scalar;

      /**
       * @brief Represents a linear form of a P0 scalar element.
       */
      class LinearForm
      {
        public:
          constexpr
          LinearForm(Geometry::Polytope::Type g)
            : m_g(g)
          {}

          constexpr
          LinearForm(const LinearForm&) = default;

          template <class T>
          constexpr
          ScalarType operator()(const T& v) const
          {
            return v(P0Element<ScalarType>(m_g).getNode(0));
          }

        private:
          const Geometry::Polytope::Type m_g;
      };

      /**
       * @brief Represents the constant basis function of a P0 scalar element.
       *
       * For P0 elements, the basis function is simply @f$ \phi(x) = 1 @f$ for all x
       * in the element. This represents a piecewise constant approximation.
       *
       * ## Properties
       * - **Value**: @f$ \phi(x) = 1 @f$ everywhere in the element
       * - **Derivatives**: @f$ \frac{\partial \phi}{\partial x_i} = 0 @f$ (constant has zero gradient)
       * - **Partition of unity**: Trivially satisfied (single basis function equals 1)
       *
       * P0 basis functions are discontinuous across element boundaries, making them
       * suitable for discontinuous Galerkin (DG) methods and flux-based formulations.
       */
      class BasisFunction
      {
        public:
          using ReturnType = Scalar;  ///< Scalar return type

          /**
           * @brief Represents a derivative of the P0 basis function (always zero).
           *
           * Since P0 elements are piecewise constant, all derivatives are zero.
           *
           * @tparam Order Order of differentiation (1, 2, ...)
           */
          template <size_t Order>
          class DerivativeFunction
          {
            public:
              constexpr
              DerivativeFunction() = default;

              constexpr
              DerivativeFunction(const DerivativeFunction&) = default;

              /**
               * @brief Evaluates the derivative (always returns 0).
               * @return Zero (constant function has zero derivative)
               */
              constexpr
              ReturnType operator()(const Math::SpatialVector<Real>&) const
              {
                return 0;
              }
          };

          constexpr
          BasisFunction() = default;

          constexpr
          BasisFunction(const BasisFunction&) = default;

          constexpr
          ReturnType operator()(const Math::SpatialVector<Real>&) const
          {
            return 1;
          }

          template <size_t Order>
          constexpr
          DerivativeFunction<Order> getDerivative(size_t) const
          {
            return DerivativeFunction<Order>();
          }
      };

      constexpr
      P0Element() = default;

      constexpr
      P0Element(Geometry::Polytope::Type geometry)
        : Parent(geometry)
      {}

      constexpr
      P0Element(const P0Element& other)
        : Parent(other)
      {}

      constexpr
      P0Element(P0Element&& other)
        : Parent(std::move(other))
      {}

      virtual constexpr ~P0Element() override = default;

      /**
       * @brief Gets the number of degrees of freedom in the finite element.
       * @returns Number of degrees of freedom
       */
      constexpr
      size_t getCount() const
      {
        return 1;
      }

      const Math::SpatialVector<Real>& getNode(size_t i) const
      {
        assert(i == 0); // P0 element has only one node
        switch (this->getGeometry())
        {
          case Geometry::Polytope::Type::Point:
          {
            static thread_local const Math::SpatialVector<Real> s_node{};
            return s_node;
          }
          case Geometry::Polytope::Type::Segment:
          {
            static thread_local const Math::SpatialVector<Real> s_node{{ 0.5 }};
            return s_node;
          }
          case Geometry::Polytope::Type::Triangle:
          {
            static thread_local const Math::SpatialVector<Real> s_node{{ Real(1) / Real(3), Real(1) / Real(3) }};
            return s_node;
          }
          case Geometry::Polytope::Type::Quadrilateral:
          {
            static thread_local const Math::SpatialVector<Real> s_node{{ 0.5, 0.5 }};
            return s_node;
          }
          case Geometry::Polytope::Type::Tetrahedron:
          {
            static thread_local const Math::SpatialVector<Real> s_node{{ 0.25, 0.25, 0.25 }};
            return s_node;
          }
          case Geometry::Polytope::Type::Wedge:
          {
            static thread_local const Math::SpatialVector<Real> s_node{{ Real(1) / Real(3), Real(1) / Real(3), 0.5 }};
            return s_node;
          }
          case Geometry::Polytope::Type::Hexahedron:
          {
            static thread_local const Math::SpatialVector<Real> s_node{{ 0.5, 0.5, 0.5 }};
            return s_node;
          }
        }
        assert(false); // Unsupported geometry
        static thread_local const Math::SpatialVector<Real> s_null{};
        return s_null;
      }

      constexpr
      LinearForm getLinearForm(size_t) const
      {
        return LinearForm(this->getGeometry());
      }

      constexpr
      BasisFunction getBasis(size_t) const
      {
        return BasisFunction();
      }

      constexpr
      size_t getOrder() const
      {
        return 0;
      }
  };

  /**
   * @ingroup FiniteElements
   * @ingroup P0ElementSpecializations
   * @brief Piecewise constant (degree 0) vector Lagrange element.
   *
   * Vector-valued P0 element with:
   * - **DOF count**: @f$ d @f$ per element (one for each vector component at barycenter)
   * - **Basis functions**: @f$ \boldsymbol{\phi}_i(x) = \mathbf{e}_j @f$ where @f$ j = i \mod d @f$ (unit vectors)
   * - **Derivatives**: @f$ \nabla \boldsymbol{\phi}_i = \mathbf{0} @f$ (zero gradient)
   * - **Structure**: Component-wise constant, with each basis having one non-zero component
   *
   * Used for vector-valued DG approximations (e.g., velocity in compressible flow, element-wise
   * constant vector fields).
   *
   * ## Example
   * For a 2D element with vdim=2:
   * - @f$ \boldsymbol{\phi}_0(x) = (1, 0)^T @f$ - constant in x-direction
   * - @f$ \boldsymbol{\phi}_1(x) = (0, 1)^T @f$ - constant in y-direction
   *
   * @tparam Scalar Type of scalar components
   */
  template <class Scalar>
  class P0Element<Math::Vector<Scalar>> final
    : public FiniteElementBase<P0Element<Math::Vector<Scalar>>>
  {
    using G = Geometry::Polytope::Type;

    public:
      /// Parent class
      using Parent = FiniteElementBase<P0Element<Math::Vector<Scalar>>>;

      using ScalarType = Scalar;

      /// Type of range
      using RangeType = Math::Vector<Scalar>;

      /**
       * @brief Represents a linear form (evaluation functional) for vector P0 elements.
       *
       * Evaluates a vector field at the barycenter and extracts a specific component.
       */
      class LinearForm
      {
        public:
          constexpr
          LinearForm()
            : m_vdim(0), m_local(0), m_g(Geometry::Polytope::Type::Point)
          {}

          constexpr
          LinearForm(size_t vdim, size_t local, Geometry::Polytope::Type g)
            : m_vdim(vdim), m_local(local), m_g(g)
          {}

          constexpr
          LinearForm(const LinearForm&) = default;

          constexpr
          LinearForm(LinearForm&&) = default;

          /**
           * @brief Evaluates the linear form on a vector field.
           * @tparam T Type of the input vector field
           * @param v Vector field to evaluate
           * @return Component value at the barycenter
           */
          template <class T>
          ScalarType operator()(const T& v) const
          {
            static thread_local RangeType s_out;
            const auto& vtx = P0Element<ScalarType>(m_g).getNode(m_local / m_vdim);
            s_out = v(vtx);
            return s_out.coeff(m_local % m_vdim);
          }

        private:
          const size_t m_vdim;               ///< Vector dimension
          const size_t m_local;              ///< Local DOF index
          const Geometry::Polytope::Type m_g;///< Geometry type
      };

      /**
       * @brief Represents a piecewise constant vector basis function.
       *
       * Each basis function is a constant unit vector in one component direction:
       * @f$ \boldsymbol{\phi}_i(x) = \mathbf{e}_j @f$ where @f$ j = i \mod d @f$.
       */
      class BasisFunction
      {
        public:
          using ReturnType = Math::Vector<ScalarType>;  ///< Vector return type

          /**
           * @brief Represents a derivative of the vector P0 basis function (always zero).
           *
           * Since P0 elements are piecewise constant, all derivatives are zero for all components.
           *
           * @tparam Order Order of differentiation (1, 2, ...)
           */
          template <size_t Order>
          class DerivativeFunction
          {
            public:
              /**
               * @brief Constructs a derivative function.
               * @param i Vector component index (0 to vdim-1)
               * @param j Spatial coordinate index for differentiation (0=x, 1=y, 2=z)
               * @param vdim Vector dimension
               * @param local Local DOF index
               * @param g Geometry type
               */
              constexpr
              DerivativeFunction(size_t i, size_t j, size_t vdim, size_t local, Geometry::Polytope::Type g)
                : m_i(i), m_j(j), m_vdim(vdim), m_local(local), m_g(g)
              {}

              constexpr
              DerivativeFunction(const DerivativeFunction&) = default;

              /**
               * @brief Evaluates the derivative (always returns 0).
               * @return Zero (constant function has zero derivative)
               */
              constexpr
              ScalarType operator()(const Math::SpatialVector<Real>&) const
              {
                return ScalarType(0);
              }

            private:
              const size_t m_i, m_j;    ///< Component and coordinate indices
              const size_t m_vdim;      ///< Vector dimension
              const size_t m_local;     ///< Local DOF index
              const Geometry::Polytope::Type m_g;  ///< Geometry type
          };

          constexpr
          BasisFunction()
            : m_vdim(0), m_local(0), m_g(Geometry::Polytope::Type::Point)
          {}

          constexpr
          BasisFunction(size_t vdim, size_t local, Geometry::Polytope::Type g)
            : m_vdim(vdim), m_local(local), m_g(g)
          {}

          constexpr
          BasisFunction(const BasisFunction&) = default;

          constexpr
          BasisFunction(BasisFunction&&) = default;

          /**
           * @brief Evaluates the vector basis function at a spatial point.
           * @return Constant unit vector: e_j where j = local % vdim
           */
          const ReturnType& operator()(const Math::SpatialVector<ScalarType>&) const
          {
            static thread_local ReturnType s_out;
            s_out.resize(m_vdim);
            s_out.setZero();
            s_out.coeffRef(m_local % m_vdim) = ScalarType(1);
            return s_out;
          }

          /**
           * @brief Gets the derivative function (always zero for P0).
           * @tparam Order Order of differentiation
           * @param i Vector component index
           * @param j Spatial coordinate index
           * @return Derivative function (returns zero)
           */
          template <size_t Order>
          constexpr
          DerivativeFunction<Order> getDerivative(size_t i, size_t j) const
          {
            return DerivativeFunction<Order>(i, j, m_vdim, m_local, m_g);
          }

        private:
          const size_t m_vdim;               ///< Vector dimension
          const size_t m_local;              ///< Local DOF index
          const Geometry::Polytope::Type m_g;///< Geometry type
      };

      P0Element()
        : Parent(Geometry::Polytope::Type::Point)
        , m_vdim(0)
      {}

      /**
       * @brief Constructor with geometry, defaulting vdim to spatial dimension.
       *
       * Keeps backward compatibility with previous behaviour.
       */
      constexpr
      P0Element(Geometry::Polytope::Type geometry)
        : P0Element(geometry, Geometry::Polytope::Traits(geometry).getDimension())
      {}

      /**
       * @brief Constructor with explicit vector dimension.
       *
       * @param geometry Polytope geometry type
       * @param vdim     Vector dimension (number of components)
       */
      constexpr
      P0Element(Geometry::Polytope::Type geometry, size_t vdim)
        : Parent(geometry)
        , m_vdim(vdim)
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

      constexpr
      P0Element(const P0Element& other)
        : Parent(other)
        , m_vdim(other.m_vdim)
        , m_lfs(other.m_lfs)
        , m_bs(other.m_bs)
      {}

      constexpr
      P0Element(P0Element&& other)
        : Parent(std::move(other))
        , m_vdim(std::exchange(other.m_vdim, 0))
        , m_lfs(std::move(other.m_lfs))
        , m_bs(std::move(other.m_bs))
      {}

      constexpr
      P0Element& operator=(const P0Element& other)
      {
        Parent::operator=(other);
        m_vdim = other.m_vdim;
        m_lfs  = other.m_lfs;
        m_bs   = other.m_bs;
        return *this;
      }

      constexpr
      P0Element& operator=(P0Element&& other)
      {
        Parent::operator=(std::move(other));
        m_vdim = std::exchange(other.m_vdim, 0);
        m_lfs  = std::move(other.m_lfs);
        m_bs   = std::move(other.m_bs);
        return *this;
      }

      constexpr
      size_t getCount() const
      {
        // One DOF per vector component (all share the same barycenter)
        return m_vdim;
      }

      constexpr
      auto getLinearForm(size_t local) const
      {
        return m_lfs[local];
      }

      constexpr
      BasisFunction getBasis(size_t local) const
      {
        return m_bs[local];
      }

      constexpr
      const Math::SpatialVector<Real>& getNode(size_t local) const
      {
        return P1Element<ScalarType>(this->getGeometry()).getNode(local / m_vdim);
      }

      constexpr
      size_t getOrder() const
      {
        return 0;
      }

    private:
      size_t m_vdim;                        ///< Vector dimension

      std::vector<LinearForm>  m_lfs;       ///< Linear forms per DOF
      std::vector<BasisFunction> m_bs;      ///< Basis functions per DOF
  };
}

#endif
