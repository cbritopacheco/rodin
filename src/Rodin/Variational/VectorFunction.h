/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
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
   */
  template <class Scalar, class Derived>
  class VectorFunctionBase : public FunctionBase<VectorFunctionBase<Scalar, Derived>>
  {
    public:
      using ScalarType = Scalar;

      using Parent = FunctionBase<VectorFunctionBase<Scalar, Derived>>;

      VectorFunctionBase() = default;

      VectorFunctionBase(const VectorFunctionBase& other)
        : Parent(other)
      {}

      VectorFunctionBase(VectorFunctionBase&& other)
        : Parent(std::move(other))
      {}

      inline
      constexpr
      VectorFunctionBase& traceOf(Geometry::Attribute attr)
      {
        Parent::traceOf(attr);
        return *this;
      }

      inline
      constexpr
      VectorFunctionBase& traceOf(const FlatSet<Geometry::Attribute>& attrs)
      {
        Parent::traceOf(attrs);
        return *this;
      }

      /**
       * @brief Convenience function to access the 1st component of the
       * vector.
       */
      inline
      constexpr
      auto x() const
      {
        assert(getDimension() >= 1);
        return operator()(0);
      }

      /**
       * @brief Convenience function to access the 2nd component of the
       * vector.
       */
      inline
      constexpr
      auto y() const
      {
        assert(getDimension() >= 2);
        return operator()(1);
      }

      inline
      constexpr
      auto z() const
      {
        assert(getDimension() >= 3);
        return operator()(2);
      }

      inline
      constexpr
      auto operator()(const Geometry::Point& p) const
      {
        return getValue(p);
      }

      /**
       * @brief Access the ith component of the vector function.
       * @returns Object of type Component<VectorFunctionBase> representing
       * the ith component of the VectorFunction.
       */
      inline
      constexpr
      auto operator()(size_t i) const
      {
        assert(0 <= i);
        assert(i < getDimension());
        return Component(*this, i);
      }

      virtual ~VectorFunctionBase() = default;

      inline
      constexpr
      RangeShape getRangeShape() const
      {
        return { getDimension(), 1 };
      }

      inline
      const Derived& getDerived() const
      {
        return static_cast<const Derived&>(*this);
      }

      inline
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return static_cast<const Derived&>(*this).getValue(p);
      }

      template <class VectorType>
      inline
      constexpr
      void getValue(VectorType& res, const Geometry::Point& p) const
      {
        if constexpr (Internal::HasGetValueMethod<Derived, VectorType&, const Geometry::Point&>::Value)
        {
          return static_cast<const Derived&>(*this).getValue(res, p);
        }
        else
        {
          res = getValue(p);
        }
      }

      /**
       * @brief Gets the dimension of the vector object.
       * @returns Dimension of vector.
       */
      inline
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

      inline
      const VectorType& getValue(const Geometry::Point& p) const
      {
        return m_vector.get();
      }

      inline
      constexpr
      void getValue(VectorType& res, const Geometry::Point&) const
      {
        res = m_vector.get();
      }

      inline
      constexpr
      size_t getDimension() const
      {
        return m_vector.get().size();
      }

      inline
      constexpr
      VectorFunction& traceOf(Geometry::Attribute attr)
      {
        return *this;
      }

      inline VectorFunction* copy() const noexcept override
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

      inline
      auto getValue(const Geometry::Point& p) const
      {
        Math::FixedSizeVector<ScalarType, 1 + sizeof...(Values)> res;
        Utility::ForIndex<1 + sizeof...(Values)>(
            [&](auto i){ res.coeffRef(static_cast<Eigen::Index>(i)) = std::get<i>(m_fs).getValue(p); });
        return res;
      }

      inline
      void getValue(VectorType& res, const Geometry::Point& p) const
      {
        res.resize(1 + sizeof...(Values));
        Utility::ForIndex<1 + sizeof...(Values)>(
            [&](auto i){ res.coeffRef(static_cast<Eigen::Index>(i)) = std::get<i>(m_fs).getValue(p); });
      }

      inline
      void getValue(FixedSizeVectorType& res, const Geometry::Point& p) const
      {
        Utility::ForIndex<1 + sizeof...(Values)>(
            [&](auto i){ res.coeffRef(static_cast<Eigen::Index>(i)) = std::get<i>(m_fs).getValue(p); });
      }

      inline
      constexpr
      size_t getDimension() const
      {
        return 1 + sizeof...(Values);
      }

      inline
      constexpr
      VectorFunction& traceOf(Geometry::Attribute attrs)
      {
        std::apply([&](auto& s) { s.traceOf(attrs); }, m_fs);
        return *this;
      }

      inline VectorFunction* copy() const noexcept override
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

      inline
      constexpr
      VectorFunction& traceOf(Geometry::Attribute)
      {
        return *this;
      }

      inline
      VectorType getValue(const Geometry::Point& p) const
      {
        if constexpr (std::is_invocable_r_v<VectorType, F, const Geometry::Point&>)
        {
          return m_f(p);
        }
        else if constexpr (std::is_invocable_r_v<void, F, VectorType&, const Geometry::Point&>)
        {
          VectorType res;
          m_f(res, p);
          return res;
        }
        else
        {
          assert(false);
          VectorType res;
          res.setConstant(NAN);
        }
      }

      inline
      void getValue(VectorType& res, const Geometry::Point& p) const
      {
        if constexpr (std::is_invocable_r_v<VectorType, F, const Geometry::Point&>)
        {
          res = m_f(p);
        }
        else if constexpr (std::is_invocable_v<F, VectorType&, const Geometry::Point&>)
        {
          m_f(res, p);
        }
        else
        {
          assert(false);
          VectorType res;
          res.setConstant(NAN);
        }
      }

      inline VectorFunction* copy() const noexcept override
      {
        return new VectorFunction(*this);
      }

    private:
      const size_t m_vdim;
      const F m_f;
  };

  template <class F,
           typename = std::enable_if_t<
             std::is_invocable_r_v<Math::Vector<Real>, F, const Geometry::Point&> ||
             std::is_invocable_v<F, Math::Vector<Real>&, const Geometry::Point&>>>
  VectorFunction(size_t, F) -> VectorFunction<F>;
}

#endif
