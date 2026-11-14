/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_MATH_VECTOR_H
#define RODIN_KOKKOS_MATH_VECTOR_H

/**
 * @file
 * @brief Owning vectors in the memory space of a Kokkos execution backend.
 */

#include <cassert>
#include <cmath>
#include <utility>
#include <Kokkos_Core.hpp>
#include <KokkosBlas1_dot.hpp>
#include <KokkosBlas1_nrm2.hpp>
#include "Rodin/Math/Vector.h"
#include "Rodin/Solver/ForwardDecls.h"
#include "Rodin/Kokkos/Context/Kokkos.h"

namespace Rodin::Kokkos::Math
{
  /**
   * @brief Dense coefficient vector with backend storage and value semantics.
   * @tparam Scalar Coefficient type.
   * @tparam Backend Kokkos backend tag or execution space.
   *
   * ## Architecture
   * The owner retains one Kokkos allocation. Copies are deep; getView() borrows
   * its allocation for kernels. Arithmetic, reductions and segment copies run
   * in the execution space. Scalar host access transfers only the requested
   * entry through a proxy; bulk host access is explicit through toHost().
   * Resizing invalidates borrowed views. The Kokkos runtime must outlive owners.
   */
  template <class Scalar, class Backend>
  class Vector
  {
    public:
      /// @brief Coefficient type.
      using ScalarType = Scalar;

      /// @brief Context selecting execution and memory.
      using ContextType = ::Rodin::Context::Kokkos<::Rodin::Context::Local, Backend>;

      /// @brief Execution space.
      using ExecutionSpace = typename ContextType::ExecutionSpaceType;

      /// @brief Allocation memory space.
      using MemorySpace = typename ContextType::MemorySpaceType;

      /// @brief Borrowed kernel view.
      using View = ::Kokkos::View<Scalar*, MemorySpace>;

      /**
       * @brief Scalar host access without a retained host mirror.
       * @details Assignment transfers an entry to the backend; reading transfers
       * an entry to the host. Use getView() for bulk kernel access.
       */
      class Reference
      {
        public:
          /// @param values Allocation. @param i Entry index.
          Reference(View values, Index i)
            : m_values(std::move(values)),
              m_i(i)
          {}

          /// @return Coefficient copied from backend memory.
          operator Scalar() const
          {
            Scalar value;
            ::Kokkos::deep_copy(value, ::Kokkos::subview(m_values, m_i));
            return value;
          }

          /// @param value New coefficient. @return This reference.
          Reference& operator=(Scalar value)
          {
            ::Kokkos::deep_copy(::Kokkos::subview(m_values, m_i), value);
            return *this;
          }

          /// @param other Source reference. @return This reference.
          Reference& operator=(const Reference& other)
          {
            return *this = Scalar(other);
          }

          /// @param value Increment. @return This reference.
          Reference& operator+=(Scalar value)
          {
            return *this = Scalar(*this) + value;
          }

          /// @param value Decrement. @return This reference.
          Reference& operator-=(Scalar value)
          {
            return *this = Scalar(*this) - value;
          }

        private:
          View m_values;
          Index m_i;
      };

      /// @brief Mutable segment supporting backend-to-backend assignment.
      class Segment
      {
        public:
          /// @param values Borrowed contiguous storage.
          explicit Segment(View values)
            : m_values(std::move(values))
          {}

          /// @param rhs Source coefficients. @return This segment.
          Segment& operator=(const Vector& rhs)
          {
            assert(m_values.extent(0) == rhs.size());
            ::Kokkos::deep_copy(m_values, rhs.getView());
            return *this;
          }

          /// @return Independent copy of the segment.
          operator Vector() const
          {
            return Vector::copy(m_values);
          }

        private:
          View m_values;
      };

      /// @brief Constructs an empty vector.
      Vector() = default;

      /// @param count Number of zero coefficients.
      explicit Vector(size_t count)
        : m_values("Rodin vector", count)
      {}

      /// @param rhs Vector to deep-copy.
      Vector(const Vector& rhs)
        : Vector(rhs.size())
      {
        ::Kokkos::deep_copy(m_values, rhs.m_values);
      }

      /// @param rhs Vector whose allocation is moved.
      Vector(Vector&& rhs) noexcept = default;

      /// @param rhs Vector to deep-copy. @return This vector.
      Vector& operator=(const Vector& rhs)
      {
        if (this != &rhs)
        {
          resize(rhs.size());
          ::Kokkos::deep_copy(m_values, rhs.m_values);
        }
        return *this;
      }

      /// @param rhs Vector whose allocation is moved. @return This vector.
      Vector& operator=(Vector&& rhs) noexcept = default;

      /// @param rhs Host expression to copy to the backend.
      template <class Derived>
      Vector(const Eigen::MatrixBase<Derived>& rhs)
      {
        *this = rhs;
      }

      /// @param rhs Host array expression to copy to the backend. @return This vector.
      template <class Derived>
      Vector& operator=(const Eigen::ArrayBase<Derived>& rhs)
      {
        return *this = rhs.matrix();
      }

      /// @param rhs Host expression to copy to the backend. @return This vector.
      template <class Derived>
      Vector& operator=(const Eigen::MatrixBase<Derived>& rhs)
      {
        resize(rhs.size());
        auto host = ::Kokkos::create_mirror(m_values);
        for (Index i = 0; i < size(); ++i)
          host(i) = rhs.derived().coeff(i);
        ::Kokkos::deep_copy(m_values, host);
        return *this;
      }

      /// @return Number of coefficients.
      size_t size() const
      {
        return m_values.extent(0);
      }

      /// @return Number of rows.
      size_t rows() const
      {
        return size();
      }

      /// @return One column.
      size_t cols() const
      {
        return 1;
      }

      /// @return Borrowed allocation handle, invalidated by replacement.
      View getView() const
      {
        return m_values;
      }

      /// @return Pointer in the backend memory space.
      Scalar* data()
      {
        return m_values.data();
      }

      /// @return Pointer in the backend memory space.
      const Scalar* data() const
      {
        return m_values.data();
      }

      /// @param count New coefficient count; unchanged sizes keep the allocation.
      void resize(size_t count)
      {
        if (size() != count)
          ::Kokkos::resize(m_values, count);
      }

      /// @brief Sets all coefficients to zero.
      void setZero()
      {
        ::Kokkos::deep_copy(m_values, Scalar(0));
      }

      /// @param count New size of the zero vector.
      void setZero(size_t count)
      {
        resize(count);
        setZero();
      }

      /// @param value Constant coefficient value.
      void setConstant(Scalar value)
      {
        ::Kokkos::deep_copy(m_values, value);
      }

      /// @param i Entry index. @return Mutable scalar proxy.
      Reference coeffRef(Index i)
      {
        assert(i < size());
        return Reference(m_values, i);
      }

      /// @param i Entry index. @return Coefficient copied to the host.
      Scalar coeff(Index i) const
      {
        assert(i < size());
        return Reference(m_values, i);
      }

      /// @param i Entry index. @return Mutable scalar proxy.
      Reference operator[](Index i)
      {
        return coeffRef(i);
      }

      /// @param i Entry index. @return Coefficient copied to the host.
      Scalar operator[](Index i) const
      {
        return coeff(i);
      }

      /// @param i Entry index. @return Mutable scalar proxy.
      Reference operator()(Index i)
      {
        return coeffRef(i);
      }

      /// @param i Entry index. @return Coefficient copied to the host.
      Scalar operator()(Index i) const
      {
        return coeff(i);
      }

      /**
       * @brief Copies a contiguous source range into this owner's storage.
       * @param source Source coefficients. @param offset First source coefficient.
       * @param count Number of coefficients to copy.
       */
      void assign(const Vector& source, size_t offset, size_t count)
      {
        assert(offset + count <= source.size());
        if (this == &source)
        {
          if (offset == 0 && count == size())
            return;
          Vector snapshot(source);
          assign(snapshot, offset, count);
          return;
        }
        resize(count);
        ::Kokkos::deep_copy(m_values,
          ::Kokkos::subview(source.getView(), std::make_pair(offset, offset + count)));
      }

      /// @param first First entry. @param count Entry count. @return Writable segment.
      Segment segment(size_t first, size_t count)
      {
        assert(first + count <= size());
        return Segment(::Kokkos::subview(m_values, std::make_pair(first, first + count)));
      }

      /// @param first First entry. @param count Entry count. @return Copied segment.
      Vector segment(size_t first, size_t count) const
      {
        assert(first + count <= size());
        return copy(::Kokkos::subview(m_values, std::make_pair(first, first + count)));
      }

      /// @return Host vector containing an independent copy.
      ::Rodin::Math::Vector<Scalar> toHost() const
      {
        const auto host =
          ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, m_values);
        ::Rodin::Math::Vector<Scalar> result(size());
        for (Index i = 0; i < size(); ++i)
          result(i) = host(i);
        return result;
      }

      /// @return Independent host coefficients for a host algebra destination.
      operator ::Rodin::Math::Vector<Scalar>() const
      {
        return toHost();
      }

      /// @return Whether every coefficient is finite.
      bool allFinite() const
      {
        Index invalid = 0;
        ::Kokkos::parallel_reduce("Rodin vector finite",
          ::Kokkos::RangePolicy<ExecutionSpace, ::Kokkos::IndexType<Index>>(0, size()),
          Finite{m_values}, invalid);
        return invalid == 0;
      }

      /// @brief Counts invalid coefficients.
      struct Finite
      {
        /// @brief Coefficients read by the reduction.
        View values;

        /// @param i Entry index. @param invalid Invalid coefficient count.
        KOKKOS_INLINE_FUNCTION void operator()(Index i, Index& invalid) const
        {
          if (::KokkosKernels::ArithTraits<Scalar>::isNan(values(i)) ||
            ::KokkosKernels::ArithTraits<Scalar>::isInf(values(i)))
            ++invalid;
        }
      };

      /// @return Euclidean coefficient norm.
      auto norm() const
      {
        return ::KokkosBlas::nrm2(m_values);
      }

      /// @param rhs Second vector, conjugated in Rodin's convention. @return Inner product.
      auto dot(const Vector& rhs) const
      {
        assert(rhs.size() == size());
        return ::KokkosBlas::dot(rhs.getView(), m_values);
      }

      /// @param tolerance Maximum Euclidean norm. @return Whether the vector is zero within tolerance.
      bool isZero(Scalar tolerance = 0) const
      {
        return norm() <= tolerance;
      }

      /// @return Squared Euclidean coefficient norm.
      auto squaredNorm() const
      {
        const auto n = norm();
        return n * n;
      }

      /// @return Coefficient-wise mutable operations.
      Vector& array()
      {
        return *this;
      }

      /// @return Coefficient-wise read access.
      const Vector& array() const
      {
        return *this;
      }

      /// @param rhs Scalar increment. @return This vector.
      Vector& operator+=(Scalar rhs)
      {
        apply({}, rhs, Operation::Add);
        return *this;
      }

      /// @param rhs Scalar decrement. @return This vector.
      Vector& operator-=(Scalar rhs)
      {
        return *this += -rhs;
      }

      /// @param rhs Scalar multiplier. @return This vector.
      Vector& operator*=(Scalar rhs)
      {
        apply({}, rhs, Operation::Multiply);
        return *this;
      }

      /// @param rhs Scalar divisor. @return This vector.
      Vector& operator/=(Scalar rhs)
      {
        return *this *= Scalar(1) / rhs;
      }

      /// @param rhs Vector increment. @return This vector.
      Vector& operator+=(const Vector& rhs)
      {
        apply(rhs.getView(), 1, Operation::Add);
        return *this;
      }

      /// @param rhs Vector decrement. @return This vector.
      Vector& operator-=(const Vector& rhs)
      {
        apply(rhs.getView(), -1, Operation::Add);
        return *this;
      }

      /// @param rhs Coefficient-wise multiplier. @return This vector.
      Vector& operator*=(const Vector& rhs)
      {
        apply(rhs.getView(), 1, Operation::Multiply);
        return *this;
      }

      /// @param rhs Coefficient-wise divisor. @return This vector.
      Vector& operator/=(const Vector& rhs)
      {
        apply(rhs.getView(), 1, Operation::Divide);
        return *this;
      }

      /// @param rhs Vector increment. @return Independent sum.
      Vector operator+(const Vector& rhs) const
      {
        Vector result(*this);
        result += rhs;
        return result;
      }

      /// @param rhs Vector decrement. @return Independent difference.
      Vector operator-(const Vector& rhs) const
      {
        Vector result(*this);
        result -= rhs;
        return result;
      }

      /// @return Independent negated coefficients.
      Vector operator-() const
      {
        return *this * Scalar(-1);
      }

      /// @param rhs Scalar factor. @return Independent scaled coefficients.
      Vector operator*(Scalar rhs) const
      {
        Vector result(*this);
        result *= rhs;
        return result;
      }

      /// @param rhs Scalar divisor. @return Independent divided coefficients.
      Vector operator/(Scalar rhs) const
      {
        return *this * (Scalar(1) / rhs);
      }

      /// @param index Receives the smallest index attaining the minimum. @return Minimum coefficient.
      Scalar minCoeff(Index* index = nullptr) const
      {
        return extremum<false>(index);
      }

      /// @param index Receives the smallest index attaining the maximum. @return Maximum coefficient.
      Scalar maxCoeff(Index* index = nullptr) const
      {
        return extremum<true>(index);
      }

      /// @brief Coefficient update operation.
      enum class Operation
      {
        Add,
        Multiply,
        Divide
      };

      /// @brief Kernel for coefficient-wise arithmetic.
      struct Transform
      {
        /// @brief Coefficients updated in place.
        View out;

        /// @brief Optional right-hand-side coefficients.
        View rhs;

        /// @brief Multiplicative scalar.
        Scalar factor;

        /// @brief Coefficient-wise arithmetic operation.
        Operation operation;

        /// @param i Coefficient index.
        KOKKOS_INLINE_FUNCTION void operator()(Index i) const
        {
          const Scalar value = rhs.extent(0) ? factor * rhs(i) : factor;
          if (operation == Operation::Add)
            out(i) += value;
          else if (operation == Operation::Multiply)
            out(i) *= value;
          else
            out(i) /= value;
        }
      };

      /// @brief Reduction kernel retaining a deterministic extremum index.
      template <bool Maximum>
      struct Extremum
      {
        /// @brief Coefficients read by the reduction.
        View values;

        /// @param i Coefficient index. @param result Current extremum.
        KOKKOS_INLINE_FUNCTION void operator()(Index i,
          typename ::Kokkos::MinFirstLoc<Scalar, Index>::value_type& result) const
        {
          const Scalar value = values(i);
          if ((Maximum ? value > result.val : value < result.val) ||
            (value == result.val && i < result.loc))
          {
            result.val = value;
            result.loc = i;
          }
        }
      };

    private:
      /// @param source Borrowed backend coefficients. @return Independent vector.
      static Vector copy(View source)
      {
        Vector result(source.extent(0));
        ::Kokkos::deep_copy(result.m_values, source);
        return result;
      }

      /// @param rhs Optional vector operand. @param factor Scalar operand. @param operation Arithmetic operation.
      void apply(View rhs, Scalar factor, Operation operation)
      {
        assert(!rhs.extent(0) || rhs.extent(0) == size());
        ::Kokkos::parallel_for("Rodin vector arithmetic",
          ::Kokkos::RangePolicy<ExecutionSpace, ::Kokkos::IndexType<Index>>(0, size()),
          Transform{m_values, rhs, factor, operation});
        ExecutionSpace().fence("Rodin vector arithmetic");
      }

      /// @tparam Maximum Select maximum rather than minimum. @param index Output index. @return Extremum.
      template <bool Maximum>
      Scalar extremum(Index* index) const
      {
        assert(size() > 0);
        typename ::Kokkos::MinFirstLoc<Scalar, Index>::value_type result;
        using Reducer = std::conditional_t<Maximum, ::Kokkos::MaxFirstLoc<Scalar, Index>,
          ::Kokkos::MinFirstLoc<Scalar, Index>>;
        ::Kokkos::parallel_reduce("Rodin vector extremum",
          ::Kokkos::RangePolicy<ExecutionSpace, ::Kokkos::IndexType<Index>>(0, size()),
          Extremum<Maximum>{m_values}, Reducer(result));
        if (index)
          *index = result.loc;
        return result.val;
      }

      View m_values;
  };

  /// @param factor Scalar factor. @param vector Backend coefficients. @return Scaled coefficients.
  template <class Factor, class Scalar, class Backend>
    requires std::is_convertible_v<Factor, Scalar>
  Vector<Scalar, Backend> operator*(Factor factor, const Vector<Scalar, Backend>& vector)
  {
    return vector * static_cast<Scalar>(factor);
  }

}

namespace Rodin::FormLanguage
{
  /// @brief Scalar traits of a backend vector.
  template <class Scalar, class Backend>
  struct Traits<::Rodin::Kokkos::Math::Vector<Scalar, Backend>>
  {
    /// @brief Coefficient scalar type.
    using ScalarType = Scalar;
  };
}

namespace Rodin::Solver::Internal
{
  /// @brief Backend vectors are read directly when host access is unavailable.
  template <class Scalar, class Backend>
  struct CuDSSVector<::Rodin::Kokkos::Math::Vector<Scalar, Backend>>
  {
    /// @brief Whether the allocation requires backend access.
    static constexpr bool Resident = !::Kokkos::SpaceAccessibility<::Kokkos::HostSpace,
      typename ::Rodin::Kokkos::Math::Vector<Scalar, Backend>::MemorySpace>::accessible;
  };
}

#endif
