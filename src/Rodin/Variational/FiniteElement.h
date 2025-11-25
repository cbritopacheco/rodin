/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file FiniteElement.h
 * @brief Base classes for finite element definitions.
 *
 * This file defines the FiniteElementBase class, which provides the foundation
 * for all finite element types. A finite element is characterized by the triple
 * @f$ (K, P, \Sigma) @f$ defining its geometry, polynomial space, and DOFs.
 *
 * ## Mathematical Foundation
 * A finite element is defined by:
 * - @f$ K @f$: Reference element geometry (interval, triangle, tetrahedron, etc.)
 * - @f$ P @f$: Polynomial space (e.g., @f$ \mathbb{P}_k @f$, @f$ \mathbb{Q}_k @f$)
 * - @f$ \Sigma @f$: Degrees of freedom (nodal values, moments, etc.)
 *
 * ## Element Types
 * Common finite element families:
 * - **Lagrange**: @f$ C^0 @f$ continuous, nodal DOFs at vertices/edges/faces
 * - **Hermite**: @f$ C^1 @f$ continuous, includes derivative DOFs
 * - **Raviart-Thomas**: @f$ H(\text{div}) @f$ conforming
 * - **Nédélec**: @f$ H(\text{curl}) @f$ conforming
 * - **Discontinuous**: @f$ L^2 @f$ elements for DG methods
 *
 * ## Polynomial Spaces
 * - @f$ \mathbb{P}_k @f$: Complete polynomials of degree @f$ \leq k @f$
 * - @f$ \mathbb{Q}_k @f$: Tensor product polynomials of degree @f$ \leq k @f$ per direction
 * - Serendipity: Subset of @f$ \mathbb{Q}_k @f$ with reduced DOFs
 *
 * ## Reference Elements
 * Standard geometries:
 * - 1D: Interval [0,1]
 * - 2D: Triangle (vertices at (0,0), (1,0), (0,1)) or quadrilateral [0,1]²
 * - 3D: Tetrahedron, hexahedron, prism, pyramid
 *
 * @see P0, P1, FiniteElementSpace
 */
#ifndef RODIN_VARIATIONAL_FINITEELEMENT_H
#define RODIN_VARIATIONAL_FINITEELEMENT_H

#include <boost/serialization/access.hpp>

#include "Rodin/Geometry/Polytope.h"

#include "ForwardDecls.h"

namespace Rodin::Variational
{
  /**
   * @defgroup FiniteElements Supported finite elements
   * @brief List of finite elements already implemented.
   */

  /**
   * @brief Base class for finite elements.
   *
   * This class provides a minimal CRTP interface shared by all finite elements.
   * Derived classes are expected to implement the following methods:
   * - `size_t getCount() const`
   * - `auto getBasis(size_t) const`
   * - `auto getLinearForm(size_t) const`
   * - `size_t getOrder() const`
   *
   * The base class stores only the reference geometry and forwards all
   * finite-element-specific queries to the Derived class.
   *
   * @tparam Derived Concrete finite element type implementing the CRTP interface.
   */
  template <class Derived>
  class FiniteElementBase
  {
    public:
      /**
       * @brief Scalar type associated with the finite element.
       *
       * This type is obtained from the @ref FormLanguage::Traits specialization
       * of the Derived finite element. It corresponds to the scalar field
       * underlying the finite element (e.g. @c Real, @c Complex<Real>, etc.).
       */
      using ScalarType = typename FormLanguage::Traits<Derived>::ScalarType;

      /**
       * @brief Default constructor.
       *
       * Initializes the finite element with a Point geometry. This is mainly
       * provided for default-constructibility and is typically overridden
       * by Derived classes that expose more specific constructors.
       */
      constexpr
      FiniteElementBase()
        : m_g(Geometry::Polytope::Type::Point)
      {}

      /**
       * @brief Constructs a finite element base with a given reference geometry.
       *
       * @param g Reference geometry of the element (segment, triangle, etc.).
       *
       * The geometry is stored and later queried via getGeometry(). Derived
       * classes may interpret this geometry to build their DOFs and basis.
       */
      constexpr
      FiniteElementBase(Geometry::Polytope::Type g)
        : m_g(g)
      {}

      /**
       * @brief Copy constructor.
       *
       * Copies the reference geometry from another base object.
       *
       * @param other Object to copy from.
       */
      constexpr
      FiniteElementBase(const FiniteElementBase& other)
        : m_g(other.m_g)
      {}

      /**
       * @brief Move constructor.
       *
       * Moves the reference geometry from another base object.
       *
       * @param other Object to move from.
       */
      constexpr
      FiniteElementBase(FiniteElementBase&& other)
        : m_g(std::move(other.m_g))
      {}

      /**
       * @brief Copy assignment operator.
       *
       * Assigns the reference geometry from another base object. The Derived
       * part is handled by the concrete finite element type.
       *
       * @param other Object to copy from.
       * @return Reference to this object.
       */
      constexpr
      FiniteElementBase& operator=(const FiniteElementBase& other)
      {
        if (this != &other)
        {
          m_g = other.m_g;
        }
        return static_cast<Derived&>(*this);
      }

      /**
       * @brief Move assignment operator.
       *
       * Moves the reference geometry from another base object. The Derived
       * part is handled by the concrete finite element type.
       *
       * @param other Object to move from.
       * @return Reference to this object.
       */
      constexpr
      FiniteElementBase& operator=(FiniteElementBase&& other)
      {
        if (this != &other)
        {
          m_g = std::move(other.m_g);
        }
        return static_cast<Derived&>(*this);
      }

      /**
       * @brief Virtual destructor.
       *
       * Provided to allow polymorphic destruction through a base-class pointer
       * (e.g. when storing different finite elements in a common container).
       */
      virtual ~FiniteElementBase() = default;

      /**
       * @brief Returns the reference geometry of the finite element.
       *
       * The geometry identifies the reference domain @f$ K @f$ used to define
       * the element (segment, triangle, quadrilateral, etc.). It is typically
       * used to:
       * - Select appropriate quadrature rules.
       * - Interpret reference coordinates.
       * - Build geometry-dependent DOFs and basis functions in Derived classes.
       *
       * @return Geometry::Polytope::Type representing the reference geometry.
       */
      constexpr
      Geometry::Polytope::Type getGeometry() const
      {
        return m_g;
      }

      /**
       * @brief Gets the number of degrees of freedom in the finite element.
       *
       * This forwards to the Derived class via CRTP:
       * @code
       * static_cast<const Derived&>(*this).getCount();
       * @endcode
       *
       * For polynomial elements, this is typically the number of interpolation
       * nodes (scalar DOFs) on the reference element.
       *
       * @return Total number of local DOFs associated with the element.
       */
      constexpr
      size_t getCount() const
      {
        return static_cast<const Derived&>(*this).getCount();
      }

      /**
       * @brief Gets the i-th basis function of the finite element.
       *
       * The basis functions @f$ \{\phi_i\} @f$ span the local approximation
       * space on the reference element. For Lagrange elements, they satisfy
       * the nodal property:
       * @f[
       *   \phi_i(x_j) = \delta_{ij},
       * @f]
       * where @f$ x_j @f$ are the interpolation nodes.
       *
       * The exact return type depends on the Derived class (scalar-valued,
       * vector-valued, etc.) and is deduced automatically.
       *
       * @param i Local basis function index.
       * @return Basis function object usable as a callable functor.
       *
       * @note CRTP method: implemented by the Derived class.
       */
      constexpr
      decltype(auto) getBasis(size_t i) const
      {
        return static_cast<const Derived&>(*this).getBasis(i);
      }

      /**
       * @brief Gets the i-th linear functional (degree of freedom) on the element.
       *
       * A linear form represents one element of the DOF set @f$ \Sigma @f$, such
       * as:
       * - Point evaluation at a node (for nodal elements),
       * - Moments against test polynomials (for non-nodal elements).
       *
       * Applying the linear form to a function @f$ v @f$ typically produces
       * the coefficient associated with the i-th degree of freedom.
       *
       * @param i Local DOF index.
       * @return Linear functional object representing the i-th degree of freedom.
       *
       * @note CRTP method: implemented by the Derived class.
       */
      constexpr
      decltype(auto) getLinearForm(size_t i) const
      {
        return static_cast<const Derived&>(*this).getLinearForm(i);
      }

      /**
       * @brief Returns a characteristic polynomial order of the finite element.
       *
       * For polynomial-based elements, this is defined as the maximum total
       * polynomial degree of the shape functions on the reference element,
       * or a safe upper bound for it. It is typically used to select
       * appropriate quadrature rules (e.g. rules exact for polynomials
       * up to degree getOrder(), or 2*getOrder() for bilinear forms).
       *
       * The precise meaning may be refined in the Derived class
       * documentation (e.g. @f$ \mathbb{P}_k @f$ vs @f$ \mathbb{Q}_k @f$
       * elements, tensor-product spaces, etc.).
       *
       * @return Polynomial order (total degree or upper bound).
       *
       * @note CRTP method: implemented by the Derived class.
       */
      constexpr
      size_t getOrder() const
      {
        return static_cast<const Derived&>(*this).getOrder();
      }

      /**
       * @brief Serializes the finite element base (Boost.Serialization).
       *
       * This function stores or loads the reference geometry. Derived classes
       * are expected to provide their own @c serialize method and call this
       * one to ensure the geometry is correctly handled.
       *
       * @tparam Archive Boost.Serialization archive type.
       * @param ar Archive to serialize to or from.
       * @param version Serialization version (unused).
       */
      template<class Archive>
      void serialize(Archive& ar, const unsigned int version)
      {
        (void) version;
        ar & m_g;
      }

    private:
      /**
       * @brief Reference geometry of the finite element.
       *
       * This encodes the type of the reference domain @f$ K @f$ (point, segment,
       * triangle, quadrilateral, tetrahedron, wedge, etc.) on which the finite
       * element is defined.
       */
      Geometry::Polytope::Type m_g;
  };
}

#endif
