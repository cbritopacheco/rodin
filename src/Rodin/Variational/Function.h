/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Function.h
 * @brief Base class hierarchy for functions in variational formulations.
 *
 * This file defines the fundamental FunctionBase template, which serves as the
 * foundation for all function types (scalar, vector, matrix) in Rodin's
 * variational formulation framework.
 */
#ifndef RODIN_VARIATIONAL_FUNCTION_H
#define RODIN_VARIATIONAL_FUNCTION_H

#include "Rodin/Cast.h"

#include "Rodin/Geometry/Point.h"

#include "Rodin/Variational/Traits.h"

#include "Rodin/FormLanguage/Base.h"
#include "Rodin/FormLanguage/Traits.h"

#include "ForwardDecls.h"

namespace Rodin::FormLanguage
{
  template <class Derived>
  struct Traits<Variational::FunctionBase<Derived>>
  {
    using ResultType = typename ResultOf<Variational::FunctionBase<Derived>>::Type;

    using RangeType = typename RangeOf<Variational::FunctionBase<Derived>>::Type;

    using ScalarType = typename FormLanguage::Traits<RangeType>::ScalarType;
  };
}

namespace Rodin::Variational
{
  /**
   * @ingroup RodinVariational
   * @brief Base class for functions defined on finite element meshes.
   *
   * FunctionBase provides the foundation for representing mathematical functions
   * defined over geometric domains discretized by finite element meshes. This
   * includes both analytical functions and discrete finite element functions
   * (grid functions) that arise in variational formulations.
   *
   * @tparam Derived Derived class following CRTP (Curiously Recurring Template Pattern)
   *
   * ## Mathematical Foundation
   * A function @f$ f : \Omega \to \mathbb{R}^n @f$ maps points in the domain 
   * @f$ \Omega @f$ to values in @f$ \mathbb{R}^n @f$. In finite element analysis,
   * functions serve various roles:
   * - **Analytical functions**: Exact solutions, boundary conditions, source terms
   * - **Discrete functions**: Finite element approximations @f$ u_h = \sum_i u_i \phi_i @f$
   * - **Test/Trial functions**: Basis functions @f$ \phi_i @f$ spanning finite element spaces
   *
   * ## Key Features
   * - **Domain support**: Functions can be restricted to subdomains or boundaries
   * - **Trace operations**: Support for function traces on mesh boundaries
   * - **Point evaluation**: Evaluation at arbitrary points within the domain
   * - **Polymorphic design**: CRTP enables compile-time polymorphism for efficiency
   */
  template <class Derived>
  class FunctionBase : public FormLanguage::Base
  {
    public:
      /// @brief Parent class type
      using Parent = FormLanguage::Base;

      /// @brief Domain type for trace operations (set of mesh attributes)
      using TraceDomain = FlatSet<Geometry::Attribute>;

      /// @brief Default constructor
      FunctionBase() = default;

      /// @brief Copy constructor
      FunctionBase(const FunctionBase& other)
        : Parent(other),
          m_traceDomain(other.m_traceDomain)
      {}

      /// @brief Move constructor
      FunctionBase(FunctionBase&& other)
        : Parent(std::move(other)),
          m_traceDomain(std::move(other.m_traceDomain))
      {}

      /// @brief Virtual destructor
      virtual ~FunctionBase() = default;

      FunctionBase& operator=(FunctionBase&& other)
      {
        m_traceDomain = std::move(other.m_traceDomain);
        return *this;
      }

      /**
       * @brief Evaluates the function on a Point belonging to the mesh.
       *
       * This operator provides convenient function call syntax for evaluation.
       * Delegates to the derived class's getValue() method via CRTP.
       *
       * @param[in] p Point at which to evaluate the function
       * @returns Function value at the given point
       */
      constexpr
      decltype(auto) operator()(const Geometry::Point& p) const
      {
        return static_cast<const Derived&>(*this).getValue(p);
      }

      /**
       * @brief Extracts the i-th component of the function.
       *
       * For vector or matrix-valued functions, this returns a scalar function
       * representing the specified component.
       *
       * @param[in] i Component index
       * @returns Component function object
       * @see Component
       */
      auto operator()(size_t i) const
      {
        return Component(*this, i);
      }

      /**
       * @brief Extracts the (i,j)-th component of a matrix function.
       *
       * @param[in] i Row index
       * @param[in] j Column index
       * @returns Component function object
       * @see Component
       */
      auto operator()(size_t i, size_t j) const
      {
        return Component(*this, i, j);
      }

      /**
       * @brief Convenience accessor for the first component (x-component).
       * @returns Component function for index 0
       */
      constexpr
      auto x() const
      {
        return Component(*this, 0);
      }

      /**
       * @brief Convenience accessor for the second component (y-component).
       * @returns Component function for index 1
       */
      constexpr
      auto y() const
      {
        return Component(*this, 1);
      }

      /**
       * @brief Convenience accessor for the third component (z-component).
       * @returns Component function for index 2
       */
      constexpr
      auto z() const
      {
        return Component(*this, 2);
      }

      /**
       * @brief Returns the transpose of the function.
       *
       * For matrix-valued functions @f$ A(x) @f$, returns @f$ A^T(x) @f$.
       *
       * @returns Transposed function object
       * @see Transpose
       */
      constexpr
      auto T() const
      {
        return Transpose(*this);
      }

      /**
       * @brief Sets a single attribute as the trace domain.
       *
       * Convenience function to call traceOf(FlatSet<int>) with only one
       * attribute. The trace operation restricts function evaluation to the
       * specified boundary or interface regions.
       *
       * @param[in] attr Mesh attribute defining the trace domain
       * @returns Reference to self (for method chaining)
       * @see getTraceDomain
       */
      constexpr
      Derived& traceOf(const Geometry::Attribute& attr)
      {
        return this->traceOf(FlatSet<Geometry::Attribute>{ attr });
      }

      /**
       * @brief Sets multiple attributes as the trace domain.
       *
       * @param[in] a1 First attribute
       * @param[in] a2 Second attribute
       * @param[in] as Additional attributes (variadic)
       * @returns Reference to self (for method chaining)
       */
      template <class A1, class A2, class ... As>
      constexpr
      Derived& traceOf(const A1& a1, const A2& a2, const As& ... as)
      {
        return this->traceOf(FlatSet<Geometry::Attribute>{ a1, a2, as... });
      }

      /**
       * @brief Sets a set of attributes as the trace domain.
       *
       * The trace domain specifies regions (typically boundaries or interfaces)
       * where the function should be evaluated via continuous extension from
       * adjacent elements.
       *
       * @param[in] attr Set of mesh attributes defining the trace domain
       * @returns Reference to self (for method chaining)
       * @see getTraceDomain
       */
      constexpr
      Derived& traceOf(const FlatSet<Geometry::Attribute>& attr)
      {
        m_traceDomain = attr;
        return static_cast<Derived&>(*this);
      }

      /**
       * @brief Gets the set of attributes which will be interpreted as the
       * domains to "trace".
       *
       * The domains to trace are interpreted as the domains where there
       * shall be a continuous extension from values to the interior
       * boundaries. If the trace domain is empty, then this has the
       * semantic value that it has not been specified yet.
       */
      constexpr
      const TraceDomain& getTraceDomain() const
      {
        return m_traceDomain;
      }

      /**
       * @brief Evaluates the function on a Point belonging to the mesh.
       * @note CRTP function to be overriden in Derived class.
       */
      constexpr
      decltype(auto) getValue(const Geometry::Point& p) const
      {
        return static_cast<const Derived&>(*this).getValue(p);
      }

      /**
       * @brief Returns a geometry-dependent polynomial order bound of the expression
       *        on the reference element.
       *
       * The returned value is a **safe upper bound** on the total polynomial degree
       * of the expression in reference coordinates, ignoring the geometry map.
       *
       * - Used for quadrature selection and composition rules.
       * - May depend on the reference geometry (simplex, tensor-product, wedge).
       * - Returns std::nullopt for non-polynomial expressions.
       * - The value is not guaranteed to be sharp.
       *
       * @param geom Reference geometry type.
       * @return Polynomial order bound, or std::nullopt if not polynomial.
       */
      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope& geom) const noexcept
      {
        return static_cast<const Derived&>(*this).getOrder(geom);
      }

      Derived& getDerived() noexcept
      {
        return static_cast<Derived&>(*this);
      }

      const Derived& getDerived() const noexcept
      {
        return static_cast<const Derived&>(*this);
      }

      virtual FunctionBase* copy() const noexcept override
      {
        return static_cast<const Derived&>(*this).copy();
      }

    private:
      FlatSet<Geometry::Attribute> m_traceDomain;
  };
}

#endif
