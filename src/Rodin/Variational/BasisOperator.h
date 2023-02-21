#ifndef RODIN_VARIATIONAL_TENSORBASIS_H
#define RODIN_VARIATIONAL_TENSORBASIS_H

#include <cassert>
#include <vector>
#include <optional>

#include "RangeShape.h"
#include "RangeType.h"

#include "Rodin/Math/Vector.h"
#include "Rodin/Math/Matrix.h"

namespace Rodin::Variational
{
  /**
   * @brief Represents a tensor basis for functions defined on finite element
   * spaces.
   *
   * Let @f$ u \in V_h @f$ be a function which has a basis representation
   * consisting of @f$ n @f$ degrees of freedom. If the value of @f$ u @f$ at a
   * point is a rank-@f$ n @f$ tensor, then this class represents a rank-@f$ (n
   * + 1) @f$ tensor @f$ T @f$. In this manner, the tensor @f$ T @f$ may be
   * visualized as a multidimensional array:
   * @f[
   *   T =
   *   \begin{bmatrix}
   *     T_1\\
   *     \vdots\\
   *     T_n
   *   \end{bmatrix}
   * @f]
   * where each @f$ T_k @f$ is a tensor of rank-@f$ n @f$ and we call it the
   * _k-th degree of freedom_.
   *
   * @note Currently, @f$ u @f$ is allowed to take rank-2 values only.
   */
  template <class T>
  class TensorBasis
  {
    public:
      using ValueType = T;

      explicit
      constexpr
      TensorBasis(std::initializer_list<T> vs)
        : m_dofs(vs.size()), m_basis(vs)
      {}

      template <class F, typename = std::enable_if_t<std::is_invocable_v<F, size_t>>>
      explicit
      constexpr
      TensorBasis(size_t dofs, F&& f)
        : m_dofs(dofs)
      {
        m_basis.reserve(dofs);
        for (size_t i = 0; i < dofs; i++)
          m_basis.push_back(f(i));
      }

      TensorBasis(const TensorBasis&) = delete;

      constexpr
      TensorBasis(TensorBasis&& other)
        : m_dofs(std::move(other.m_dofs)),
          m_basis(std::move(other.m_basis))
      {}

      void operator=(const TensorBasis&) = delete;

      template <class F, typename = std::enable_if_t<std::is_invocable_v<F, const T&>>>
      inline
      constexpr
      auto apply(F&& f)
      {
        using R = std::invoke_result_t<F, const T&>;
        return TensorBasis(m_dofs, [&](size_t i) -> R { return f(m_basis[i]); });
      }

      TensorBasis& operator=(TensorBasis&& other)
      {
        m_dofs = std::move(other.m_basis);
        m_basis = std::move(other.m_basis);
      }

      inline
      constexpr
      size_t getDOFs() const
      {
        return m_dofs;
      }

      inline
      constexpr
      const T& operator()(size_t i) const
      {
        assert(m_basis.size() > i);
        return m_basis[i];
      }

      inline
      constexpr
      auto begin() const
      {
        return m_basis.begin();
      }

      inline
      constexpr
      auto end() const
      {
        return m_basis.end();
      }

      const std::vector<T>& getData() const
      {
        return m_basis;
      }

    private:
      const size_t m_dofs;
      std::vector<T> m_basis;
  };

  template <>
  class TensorBasis<Scalar> : public Math::Vector
  {
    public:
      using ValueType = Scalar;
      using Parent = Math::Vector;

      explicit
      TensorBasis(std::initializer_list<Scalar> vs)
        : Parent(vs.size())
      {}

      explicit
      TensorBasis(const Math::Vector& other)
        : Parent(other)
      {}

      explicit
      TensorBasis(Math::Vector&& other)
        : Parent(std::move(other))
      {}

      template <class F, typename = std::enable_if_t<std::is_invocable_v<F, size_t>>>
      explicit
      constexpr
      TensorBasis(size_t dofs, F&& f)
        : Parent(dofs)
      {
        static_assert(std::is_same_v<Scalar, std::invoke_result_t<F, size_t>>);
        for (size_t i = 0; i < dofs; i++)
          Parent::coeffRef(static_cast<Math::Vector::Index>(i)) = f(i);
      }

      TensorBasis(const TensorBasis& other)
        : Parent(other)
      {}

      TensorBasis(TensorBasis&& other)
        : Parent(std::move(other))
      {}

      void operator=(const TensorBasis&) = delete;

      template <class F, typename = std::enable_if_t<std::is_invocable_v<F, const Scalar&>>>
      inline
      constexpr
      auto apply(F&& f)
      {
        using R = std::invoke_result_t<F, const Scalar&>;
        return TensorBasis(size(), [&](size_t i) -> R { return f(operator()(i)); });
      }

      TensorBasis& operator=(TensorBasis&& other)
      {
        Parent::operator=(std::move(other));
        return *this;
      }

      inline
      constexpr
      size_t getDOFs() const
      {
        return size();
      }

      inline
      constexpr
      Scalar operator()(size_t i) const
      {
        return Parent::operator()(static_cast<Math::Vector::Index>(i));
      }

      const Scalar* getData() const
      {
        return Parent::data();
      }
  };

  template <class T>
  TensorBasis(std::initializer_list<T>) -> TensorBasis<T>;

  template <class T>
  TensorBasis(size_t, const T&) -> TensorBasis<T>;

  template <class F, typename = std::enable_if_t<std::is_invocable_v<F, size_t>>>
  TensorBasis(size_t, F&&) -> TensorBasis<std::invoke_result_t<F, size_t>>;

  TensorBasis(const Math::Vector&) -> TensorBasis<Scalar>;

  TensorBasis(Math::Vector&&) -> TensorBasis<Scalar>;

  template <class LHS, class RHS>
  inline
  constexpr
  auto operator+(const TensorBasis<LHS>& lhs, const TensorBasis<RHS>& rhs)
  {
    assert(lhs.getDOFs() == rhs.getDOFs());
    return TensorBasis(lhs.getDOFs(), [&](size_t i) { return lhs(i) + rhs(i); });
  }

  template <class LHS, class RHS>
  inline
  constexpr
  auto operator-(const TensorBasis<LHS>& lhs, const TensorBasis<RHS>& rhs)
  {
    assert(lhs.getDOFs() == rhs.getDOFs());
    return TensorBasis(lhs.getDOFs(), [&](size_t i) { return lhs(i) - rhs(i); });
  }

  template <class M>
  inline
  constexpr
  auto operator-(const TensorBasis<M>& op)
  {
    return TensorBasis(op.getDOFs(), [&](size_t i) { return -op(i); });
  }

  template <class LHS, class RHS>
  inline
  constexpr
  auto operator*(const TensorBasis<LHS>& lhs, const TensorBasis<RHS>& rhs)
  {
    assert(lhs.getDOFs() == rhs.getDOFs());
    return TensorBasis(lhs.getDOFs(), [&](size_t i) { return lhs(i) * rhs(i); });
  }

  template <class RHS>
  inline
  constexpr
  auto operator*(Scalar lhs, const TensorBasis<RHS>& rhs)
  {
    return TensorBasis(rhs.getDOFs(), [&](size_t i) { return lhs * rhs(i); });
  }

  template <class LHS>
  inline
  constexpr
  auto operator*(const TensorBasis<LHS>& lhs, Scalar rhs)
  {
    return TensorBasis(lhs.getDOFs(), [&](size_t i) { return lhs(i) * rhs; });
  }

  template <class LHS>
  inline
  constexpr
  auto operator/(const TensorBasis<LHS>& lhs, Scalar rhs)
  {
    return TensorBasis(lhs.getDOFs(), [&](size_t i) { return lhs(i) / rhs; });
  }
}

#endif
