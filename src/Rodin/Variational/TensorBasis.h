#ifndef RODIN_VARIATIONAL_TENSORBASIS_H
#define RODIN_VARIATIONAL_TENSORBASIS_H

#include <cassert>
#include <vector>
#include <optional>

#include "RangeShape.h"
#include "RangeType.h"

#include "Rodin/Math/Vector.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/Tensor.h"

#include "Rodin/FormLanguage/Traits.h"

namespace Rodin::FormLanguage
{
  template <class T>
  struct Traits<Variational::TensorBasis<T>>
  {
    using BasisType = T;
  };
}

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
  class TensorBasis final
  {
    public:
      using BasisType = T;

      template <class F, typename = std::enable_if_t<std::is_invocable_v<F, size_t>>>
      constexpr
      TensorBasis(size_t dofs, F&& f)
        : m_dofs(dofs)
      {
        assert(dofs > 0);
        m_data.reserve(dofs);
        for (size_t i = 0; i < dofs; i++)
          m_data.push_back(f(i));
      }

      constexpr
      TensorBasis(const TensorBasis&) = default;

      constexpr
      TensorBasis(TensorBasis&& other) = default;

      constexpr
      TensorBasis& operator=(const TensorBasis&) = default;

      constexpr
      TensorBasis& operator=(TensorBasis&&) = default;

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
        assert(m_data.size() > i);
        return m_data[i];
      }

      inline
      constexpr
      auto begin() const
      {
        return m_data.begin();
      }

      inline
      constexpr
      auto end() const
      {
        return m_data.end();
      }

    private:
      const size_t m_dofs;
      std::vector<T> m_data;
  };

  template <class F, typename = std::enable_if_t<std::is_invocable_v<F, size_t>>>
  TensorBasis(size_t, F&&) -> TensorBasis<std::invoke_result_t<F, size_t>>;

  template <class LHS, class RHS>
  inline
  constexpr
  auto operator+(const TensorBasis<LHS>& lhs, const TensorBasis<RHS>& rhs)
  {
    assert(lhs.getDOFs() == rhs.getDOFs());
    return TensorBasis(lhs.getDOFs(), [&](size_t i){ return lhs(i) + rhs(i); });
  }

  template <class LHS, class RHS>
  inline
  constexpr
  auto operator-(const TensorBasis<LHS>& lhs, const TensorBasis<RHS>& rhs)
  {
    assert(lhs.getDOFs() == rhs.getDOFs());
    return TensorBasis(lhs.getDOFs(), [&](size_t i){ return lhs(i) - rhs(i); });
  }

  template <class M>
  inline
  constexpr
  auto operator-(const TensorBasis<M>& op)
  {
    return TensorBasis(op.getDOFs(), [&](size_t i){ return -op(i); });
  }

  template <class RHS>
  inline
  constexpr
  auto operator*(Scalar lhs, const TensorBasis<RHS>& rhs)
  {
    return TensorBasis(rhs.getDOFs(), [&](size_t i){ return lhs * rhs(i); });
  }

  template <class LHS>
  inline
  constexpr
  auto operator*(const TensorBasis<LHS>& lhs, Scalar rhs)
  {
    return TensorBasis(lhs.getDOFs(), [&](size_t i){ return lhs(i) * rhs; });
  }

  template <class LHS>
  inline
  constexpr
  auto operator/(const TensorBasis<LHS>& lhs, Scalar rhs)
  {
    return TensorBasis(lhs.getDOFs(), [&](size_t i){ return lhs(i) / rhs; });
  }
}

#endif
