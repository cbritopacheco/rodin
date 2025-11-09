/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file VectorFunction.h
 * @brief Vector-valued functions for variational formulations.
 *
 * This file defines VectorFunctionBase and VectorFunction for representing
 * functions mapping points to vectors: @f$ \mathbf{f}: \Omega \to \mathbb{R}^d @f$.
 * These are used for vector fields such as velocity, displacement, or force fields.
 */
#ifndef RODIN_VARIATIONAL_VECTORFUNCTION_H
#define RODIN_VARIATIONAL_VECTORFUNCTION_H

#include <memory>
#include <optional>
#include <type_traits>

#include "ForwardDecls.h"

#include "Rodin/Alert.h"
#include "Rodin/Utility/ForConstexpr.h"

#include "Function.h"
#include "RealFunction.h"

namespace Rodin::FormLanguage
{
  template <class Scalar, class Derived>
  struct Traits<Variational::VectorFunctionBase<Scalar, Derived>>
  {
    using ScalarType = Scalar;
    using DerivedType = Derived;
  };
}

namespace Rodin::Variational
{
  /**
   * @defgroup VectorFunctionSpecializations VectorFunction Template Specializations
   * @brief Template specializations of the VectorFunction class.
   * @see VectorFunction
   */

  /**
   * @brief Base class for vector-valued functions defined on a mesh.
   *
   * VectorFunctionBase extends FunctionBase to represent vector-valued functions:
   * @f[
   *    \mathbf{f}: \Omega \to \mathbb{R}^d
   * @f]
   * where @f$ d @f$ is the dimension of the vector space.
   *
   * These functions are fundamental in finite element analysis for representing:
   * - **Vector fields**: Velocity @f$ \mathbf{v}(x) @f$, displacement @f$ \mathbf{u}(x) @f$
   * - **Force fields**: Body forces, surface tractions
   * - **Gradients**: @f$ \nabla u @f$ for scalar fields
   * - **Fluxes**: Heat flux, mass flux
   *
   * @tparam Scalar The scalar component type (typically Real or Complex)
   * @tparam Derived The derived class following CRTP pattern
   *
   * ## Component Access
   * Vector components can be accessed via:
   * - Index notation: `f(i)` for the i-th component
   * - Coordinate accessors: `f.x()`, `f.y()`, `f.z()` for first three components
   *
   * @see FunctionBase, RealFunction, MatrixFunction
   */
  template <class Scalar, class Derived>
  class VectorFunctionBase : public FunctionBase<VectorFunctionBase<Scalar, Derived>>
  {
    public:
      /// @brief Type of scalar components
      using ScalarType = Scalar;

      /// @brief Parent class type
      using Parent = FunctionBase<VectorFunctionBase<Scalar, Derived>>;

      /// @brief Import operator() from parent
      using Parent::operator();

      /// @brief Default constructor
      VectorFunctionBase() = default;

      /// @brief Copy constructor
      /// @param[in] other Vector function to copy from
      VectorFunctionBase(const VectorFunctionBase& other)
        : Parent(other)
      {}

      /// @brief Move constructor
      /// @param[in] other Vector function to move from
      VectorFunctionBase(VectorFunctionBase&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Convenience accessor for the x-component (first component).
       * @returns Component function for index 0
       * @pre getDimension() >= 1
       */
      constexpr
      auto x() const
      {
        assert(getDimension() >= 1);
        return operator()(0);
      }

      /**
       * @brief Convenience accessor for the y-component (second component).
       * @returns Component function for index 1
       * @pre getDimension() >= 2
       */
      constexpr
      auto y() const
      {
        assert(getDimension() >= 2);
        return operator()(1);
      }

      /**
       * @brief Convenience accessor for the z-component (third component).
       * @returns Component function for index 2
       * @pre getDimension() >= 3
       */
      constexpr
      auto z() const
      {
        assert(getDimension() >= 3);
        return operator()(2);
      }

      /**
       * @brief Accesses the i-th component of the vector function.
       *
       * Returns a scalar function representing the specified component:
       * @f$ f_i(x) = (\mathbf{f}(x))_i @f$
       *
       * @param[in] i Component index (0-based)
       * @returns Component function object
       * @pre i < getDimension()
       * @see Component
       */
      constexpr
      auto operator()(size_t i) const
      {
        assert(i < getDimension());
        return Component(*this, i);
      }

      virtual ~VectorFunctionBase() = default;

      constexpr
      decltype(auto) getValue(const Geometry::Point& p) const
      {
        return static_cast<const Derived&>(*this).getValue(p);
      }

      /**
       * @brief Gets the dimension of the vector object.
       * @returns Dimension of vector.
       */
      constexpr
      size_t getDimension() const
      {
        return static_cast<const Derived&>(*this).getDimension();
      }

      virtual VectorFunctionBase* copy() const noexcept override
      {
        return static_cast<const Derived&>(*this).copy();
      }
  };

  template <class Scalar>
  class VectorFunction<Math::Vector<Scalar>> final
    : public VectorFunctionBase<Scalar, VectorFunction<Math::Vector<Scalar>>>
  {
    public:
      using ScalarType = Scalar;

      using VectorType = Math::Vector<ScalarType>;

      using Parent = VectorFunctionBase<ScalarType, VectorFunction<VectorType>>;

      /**
       * @brief Constructs a vector with the given values.
       * @param[in] values Parameter pack of values
       *
       * Each value passed must be convertible to any specialization of
       * RealFunction.
       */
      VectorFunction(const VectorType& v)
        : m_vector(v)
      {}

      VectorFunction(const VectorFunction& other)
        : Parent(other),
          m_vector(other.m_vector)
      {}

      VectorFunction(VectorFunction&& other)
        : Parent(std::move(other)),
          m_vector(std::move(other.m_vector))
      {}

      const VectorType& getValue(const Geometry::Point& p) const
      {
        return m_vector.get();
      }

      constexpr
      size_t getDimension() const
      {
        return m_vector.get().size();
      }

      constexpr
      VectorFunction& traceOf(const FlatSet<Geometry::Attribute>& attr)
      {
        return *this;
      }

      VectorFunction* copy() const noexcept override
      {
        return new VectorFunction(*this);
      }

    private:
      std::reference_wrapper<const VectorType> m_vector;
  };

  template <class Scalar>
  VectorFunction(const Math::Vector<Scalar>&) -> VectorFunction<Math::Vector<Scalar>>;

  /**
   * @ingroup VectorFunctionSpecializations
   * @tparam V Type of first value
   * @tparam Values Parameter pack of remaining values
   * @brief Represents a vector function which may be constructed from values
   * which can be converted to objects of type RealFunction.
   *
   * In general one may construct any VectorFunction by specifying its values
   * in a uniform initialization manner. For example, to construct a
   * VectorFunction with constant entries (1, 2, 3) :
   * @code{.cpp}
   * auto v = VectorFunction{1, 2, 3};
   * @endcode
   * Alternatively, we may construct instances of VectorFunction from any type
   * which is convertible to specializations of RealFunction:
   * @code{.cpp}
   * auto s = RealFunction(3.1416);
   * auto v = VectorFunction{Dx(s), 42, s};
   * @endcode
   */
  template <class V, class ... Values>
  class VectorFunction<V, Values...> final
    : public VectorFunctionBase<Real, VectorFunction<V, Values...>>
  {
    public:
      using ScalarType = Real;

      using VectorType = Math::Vector<ScalarType>;

      using FixedSizeVectorType = Math::FixedSizeVector<ScalarType, 1 + sizeof...(Values)>;

      using Parent = VectorFunctionBase<ScalarType, VectorFunction<V, Values...>>;
      /**
       * @brief Constructs a vector with the given values.
       * @param[in] values Parameter pack of values
       *
       * Each value passed must be convertible to any specialization of
       * RealFunction.
       */
      VectorFunction(const V& v, const Values&... values)
        : m_fs(RealFunction(v), RealFunction(values)...)
      {}

      VectorFunction(const VectorFunction& other)
        : Parent(other),
          m_fs(other.m_fs)
      {}

      VectorFunction(VectorFunction&& other)
        : Parent(std::move(other)),
          m_fs(std::move(other.m_fs))
      {}

      decltype(auto) getValue(const Geometry::Point& p) const
      {
        static thread_local Math::FixedSizeVector<ScalarType, 1 + sizeof...(Values)> s_res;
        Utility::ForIndex<1 + sizeof...(Values)>(
          [&](auto i)
          {
            s_res.coeffRef(static_cast<Eigen::Index>(i)) = std::get<i>(m_fs).getValue(p);
          });
        return s_res;
      }

      constexpr
      size_t getDimension() const
      {
        return 1 + sizeof...(Values);
      }

      template <class ... Args>
      constexpr
      VectorFunction& traceOf(const Args&... attrs)
      {
        std::apply([&](auto&... s) { (s.traceOf(attrs...), ...); }, m_fs);
        return *this;
      }

      VectorFunction* copy() const noexcept override
      {
        return new VectorFunction(*this);
      }

    private:
      std::tuple<RealFunction<V>, RealFunction<Values>...> m_fs;
  };

  template <class V, class ... Values>
  VectorFunction(const V&, const Values&...) -> VectorFunction<V, Values...>;

  template <class F>
  class VectorFunction<F> final : public VectorFunctionBase<Real, VectorFunction<F>>
  {
    public:
      using ScalarType = Real;

      using VectorType = Math::Vector<ScalarType>;

      using Parent = VectorFunctionBase<ScalarType, VectorFunction<F>>;

      using Parent::traceOf;

      VectorFunction(size_t vdim, F f)
        : m_vdim(vdim), m_f(f)
      {}

      VectorFunction(const VectorFunction& other)
        : Parent(other),
          m_vdim(other.m_vdim),
          m_f(other.m_f)
      {}

      VectorFunction(VectorFunction&& other)
        : Parent(std::move(other)),
          m_vdim(std::move(other.m_vdim)),
          m_f(std::move(other.m_f))
      {}

      decltype(auto) getValue(const Geometry::Point& p) const
      {
        return m_f(p);
      }

      constexpr
      VectorFunction& traceOf(const FlatSet<Geometry::Attribute>& attr)
      {
        return *this;
      }

      VectorFunction* copy() const noexcept override
      {
        return new VectorFunction(*this);
      }

    private:
      const size_t m_vdim;
      const F m_f;
  };

  template <class F, typename = std::enable_if_t<std::is_invocable_v<F, const Geometry::Point&>>>
  VectorFunction(size_t, F) -> VectorFunction<F>;
}

#endif
