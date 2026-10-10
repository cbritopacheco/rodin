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

#include <functional>
#include <optional>
#include <type_traits>

#include <Eigen/Core>

#include "Rodin/Cast.h"

#include "Rodin/Geometry/Point.h"

#include "Rodin/Variational/Traits.h"

#include "Rodin/FormLanguage/Base.h"
#include "Rodin/FormLanguage/Traits.h"

#include "ForwardDecls.h"
#include "IntegrationPoint.h"

namespace Rodin::FormLanguage
{
  /// @brief Form-language traits for variational functions.
  template <class Derived>
  struct Traits<Variational::FunctionBase<Derived>>
  {
      /// @brief Result type of the evaluation.
      using ResultType = typename ResultOf<Variational::FunctionBase<Derived>>::Type;

      /// @brief Range (evaluation value) type.
      using RangeType = typename RangeOf<Variational::FunctionBase<Derived>>::Type;

      /// @brief Scalar value type.
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

      /**
       * @brief Copy constructor
       * @param other Object to copy from.
       */
      FunctionBase(const FunctionBase& other)
        : Parent(other),
          m_traceDomain(other.m_traceDomain)
      {}

      /**
       * @brief Move constructor
       * @param other Object to move from.
       */
      FunctionBase(FunctionBase&& other)
        : Parent(std::move(other)),
          m_traceDomain(std::move(other.m_traceDomain))
      {}

      /// @brief Virtual destructor
      virtual ~FunctionBase() = default;

      /**
       * @brief Move assignment operator.
       * @param other Object to move from.
       * @returns Reference to this object after the operation.
       */
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
      auto operator()(const Geometry::Point& p) const
      {
        return static_cast<const Derived&>(*this).getValue(p);
      }

      /**
       * @brief Evaluates the function at an integration point.
       * @param ip Integration point at which the expression is evaluated.
       * @returns Value of the expression at the supplied evaluation point.
       */
      constexpr
      auto operator()(const IntegrationPoint& ip) const
      {
        return getValue(ip);
      }

      /**
       * @brief Extracts the i-th component of the function.
       *
       * For vector or matrix-valued functions, this returns a scalar function
       * representing the specified component.
       *
       * @param[in] i Component index
       * @returns Component function object
       * @see <a href="_component_8h.html">Component</a>
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
       * @see <a href="_component_8h.html">Component</a>
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
       * @see <a href="_transpose_8h.html">Transpose</a>
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
       * @see getTraceDomain()
       */
      constexpr
      Derived& traceOf(const Geometry::Attribute& attr)
      {
        return this->traceOf(FlatSet<Geometry::Attribute>{ attr });
      }

      /**
       * @brief Sets multiple attributes as the trace domain.
       *
       * The attributes are collected into the trace domain of the function.
       *
       * @returns Reference to self (for method chaining)
       * @param a1 Mesh attributes selecting the region.
       * @param a2 Mesh attributes selecting the region.
       * @param as Mesh attributes selecting the region.
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
       * @see getTraceDomain()
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
       * @returns The set of attributes which will be interpreted as the domains to "trace".
       */
      constexpr
      const TraceDomain& getTraceDomain() const
      {
        return m_traceDomain;
      }

      /**
       * @brief Evaluates the function on a Point belonging to the mesh.
       * @note CRTP function to be overriden in Derived class.
       * @param p Point at which the operation is evaluated.
       * @returns Value of the expression at the supplied evaluation point.
       */
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return static_cast<const Derived&>(*this).getValue(p);
      }

      /**
       * @brief Evaluates the function at an integration point.
       * @param ip Integration point at which the expression is evaluated.
       * @returns Value of the expression at the supplied evaluation point.
       */
      constexpr
      auto getValue(const IntegrationPoint& ip) const
      {
        if constexpr (requires (const Derived& f, const IntegrationPoint& q) { f.getValue(q); })
          return static_cast<const Derived&>(*this).getValue(ip);
        else
          return static_cast<const Derived&>(*this).getValue(ip.getPoint());
      }

      /**
       * @brief Function value retained for the current quadrature binding.
       *
       * Holds the function value at the quadrature point last passed to
       * setIntegrationPoint(), allowing shape expressions to reuse it across basis indices.
       * The owning shape expression refreshes this snapshot on every point
       * binding, including reassembly at an unchanged point. Direct pointwise
       * bindings clear it and retain normal function evaluation. Only owning,
       * copyable values are held; lazy Eigen expressions retain direct evaluation.
       * Each shape expression owns a separate cache so evaluation passes do not
       * share mutable state.
       */
      class Cache
      {
        public:
          /// @brief Type returned when evaluating the enclosing function.
          using Value =
            std::decay_t<decltype(std::declval<const FunctionBase&>().getValue(
              std::declval<const IntegrationPoint&>()))>;

          /// @brief Declared range of the enclosing function.
          using Range = typename FormLanguage::Traits<FunctionBase>::RangeType;

          /// @brief Whether an owning, assignable snapshot can be retained.
          static constexpr bool Enabled =
            (std::is_same_v<Value, Range> ||
              std::is_base_of_v<Eigen::PlainObjectBase<Value>, Value>) &&
            std::is_copy_constructible_v<Value> && std::is_copy_assignable_v<Value>;

          /**
           * @brief Constructs an empty cache bound to a function.
           * @param[in] f Function to evaluate on subsequent point bindings.
           * @note The function is not owned and must outlive this cache.
           */
          explicit Cache(const FunctionBase& f)
            : m_function(f)
          {}

          /**
           * @brief Prevents binding a cache to a temporary function.
           * @param[in] f Temporary function whose lifetime cannot cover the cache.
           */
          Cache(const FunctionBase&& f) = delete;

          /**
           * @brief Constructs an empty cache when copying an evaluation pass.
           * @param[in] other Source cache whose function binding is copied, but whose snapshot is not.
           */
          Cache(const Cache& other)
            : Cache(other.m_function.get())
          {}

          /**
           * @brief Transfers the function binding and current snapshot.
           * @param[in] other Source cache to move from.
           */
          Cache(Cache&& other) = default;

          /**
           * @brief Clears the snapshot when copying another evaluation pass.
           * @param[in] other Source cache whose function binding is copied, but whose snapshot is not.
           * @returns This cache with an empty snapshot.
           */
          Cache& operator=(const Cache& other)
          {
            m_function = other.m_function;
            m_value.reset();
            return *this;
          }

          /**
           * @brief Transfers the function binding and snapshot on move assignment.
           * @param[in] other Source cache to move from.
           * @returns This cache.
           */
          Cache& operator=(Cache&& other) = default;

          /**
           * @brief Evaluates the bound function at an integration point.
           * @param[in] ip Evaluation point; no reference to it is retained.
           * @returns This cache.
           *
           * When @ref Enabled is true and @p ip has quadrature metadata, evaluates
           * the bound function and owns its value. Every call evaluates again,
           * even at the same point, so reassembly observes changed function data. A point without
           * quadrature metadata clears the snapshot. With @ref Enabled false,
           * leaves the cache empty and does not evaluate the function.
           */
          Cache& setIntegrationPoint(const IntegrationPoint& ip)
          {
            if constexpr (Enabled)
            {
              if (!ip.getQuadratureFormula())
              {
                m_value.reset();
                return *this;
              }
              if (m_value)
                *m_value = m_function.get().getValue(ip);
              else
                m_value.emplace(m_function.get().getValue(ip));
            }
            return *this;
          }

          /**
           * @brief Gets the snapshot for the current binding.
           * @returns Pointer to the owned value, or @c nullptr when empty or disabled.
           * @note The pointer is valid until this cache is rebound, assigned,
           * moved from or destroyed. Consumers must not retain it across bindings.
           */
          const Value* get() const
          {
            if constexpr (Enabled)
            {
              if (m_value)
                return &*m_value;
            }
            return nullptr;
          }

        private:
          std::reference_wrapper<const FunctionBase> m_function;
          std::optional<Value> m_value;
      };

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

      /**
       * @brief Returns this object as the CRTP-derived type.
       * @returns This object as the CRTP-derived type.
       */
      Derived& getDerived() noexcept
      {
        return static_cast<Derived&>(*this);
      }

      /**
       * @brief Returns this object as the CRTP-derived type.
       * @returns This object as the CRTP-derived type.
       */
      const Derived& getDerived() const noexcept
      {
        return static_cast<const Derived&>(*this);
      }

      /**
       * @brief Polymorphically copies the derived function.
       * @returns Pointer to a newly allocated copy; the caller owns the returned object.
       */
      virtual FunctionBase* copy() const noexcept override
      {
        return static_cast<const Derived&>(*this).copy();
      }

    private:
      FlatSet<Geometry::Attribute> m_traceDomain;
  };

  /**
   * @brief Returns the order only when a function is elementwise constant.
   * @param f Function operand.
   * @param polytope Mesh entity used by this operation.
   * @returns Zero when the function is known to be elementwise constant, or an empty optional otherwise.
   */
  template <class Derived>
  inline Optional<size_t>
  GetOrderIfConstant(const FunctionBase<Derived>& f, const Geometry::Polytope& polytope) noexcept
  {
    const auto o = f.getOrder(polytope);
    if (o && *o == 0)
      return o;
    return std::nullopt;
  }
}

#endif
