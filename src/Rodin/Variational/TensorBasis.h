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
      using ValueType = T;

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

  // template <>
  // class TensorBasis<Scalar> final
  // {
  //   public:
  //     using ValueType = Scalar;

  //     TensorBasis() = default;

  //     template <class EigenDerived>
  //     TensorBasis(const Eigen::MatrixBase<EigenDerived>& v)
  //       : m_data(v)
  //     {}

  //     template <class EigenDerived>
  //     TensorBasis(Eigen::MatrixBase<EigenDerived>&& v)
  //       : m_data(std::move(v))
  //     {}

  //     template <class F, typename = std::enable_if_t<std::is_invocable_v<F, size_t>>>
  //     constexpr
  //     TensorBasis(size_t dofs, F&& f)
  //       : m_data(dofs)
  //     {
  //       assert(dofs > 0);
  //       for (size_t i = 0; i < dofs; i++)
  //         m_data.coeffRef(static_cast<Math::Vector::Index>(i)) = f(i);
  //     }

  //     TensorBasis(const TensorBasis&) = default;

  //     TensorBasis(TensorBasis&&) = default;

  //     TensorBasis<Scalar>& operator=(const TensorBasis&) = default;

  //     TensorBasis& operator=(TensorBasis&&) = default;

  //     inline
  //     constexpr
  //     size_t getDOFs() const
  //     {
  //       return m_data.size();
  //     }

  //     const Math::Vector& getVector() const
  //     {
  //       return m_data;
  //     }

  //     inline
  //     Scalar operator()(size_t i) const
  //     {
  //       return m_data.coeff(static_cast<Math::Vector::Index>(i));
  //     }

  //     friend inline TensorBasis operator+(const TensorBasis&, const TensorBasis&);
  //     friend inline TensorBasis operator-(const TensorBasis&, const TensorBasis&);
  //     friend inline TensorBasis operator-(const TensorBasis&);
  //     friend inline TensorBasis operator*(Scalar, const TensorBasis& rhs);
  //     friend inline TensorBasis operator*(const TensorBasis&, Scalar);
  //     friend inline TensorBasis operator/(const TensorBasis&, Scalar);

  //   private:
  //     Math::Vector m_data;
  // };

  // inline
  // TensorBasis<Scalar> operator+(const TensorBasis<Scalar>& lhs, const TensorBasis<Scalar>& rhs)
  // {
  //   assert(lhs.getDOFs() == rhs.getDOFs());
  //   return lhs.m_data + rhs.m_data;
  // }

  // inline
  // TensorBasis<Scalar> operator-(const TensorBasis<Scalar>& lhs, const TensorBasis<Scalar>& rhs)
  // {
  //   assert(lhs.getDOFs() == rhs.getDOFs());
  //   return lhs.m_data - rhs.m_data;
  // }

  // inline
  // TensorBasis<Scalar> operator-(const TensorBasis<Scalar>& op)
  // {
  //   return -op.m_data;
  // }

  // inline
  // TensorBasis<Scalar> operator*(Scalar lhs, const TensorBasis<Scalar>& rhs)
  // {
  //   return lhs * rhs.m_data;
  // }

  // inline
  // TensorBasis<Scalar> operator*(const TensorBasis<Scalar>& lhs, Scalar rhs)
  // {
  //   return lhs.m_data * rhs;
  // }

  // inline
  // TensorBasis<Scalar> operator/(const TensorBasis<Scalar>& lhs, Scalar rhs)
  // {
  //   return lhs.m_data / rhs;
  // }

  // template <>
  // class TensorBasis<Math::Vector> final
  // {
  //   public:
  //     using ValueType = Math::Vector;

  //     TensorBasis() = default;

  //     template <class EigenDerived>
  //     TensorBasis(const Eigen::MatrixBase<EigenDerived>& v)
  //       : m_data(v)
  //     {}

  //     template <class EigenDerived>
  //     TensorBasis(Eigen::MatrixBase<EigenDerived>&& v)
  //       : m_data(std::move(v))
  //     {}

  //     template <class F, typename = std::enable_if_t<std::is_invocable_v<F, size_t>>>
  //     constexpr
  //     TensorBasis(size_t dim, size_t dofs, F&& f)
  //       : m_data(dim, dofs)
  //     {
  //       assert(dofs > 0);
  //       for (size_t i = 0; i < dofs; i++)
  //         m_data.col(i) = f(i);
  //     }

  //     TensorBasis(const TensorBasis&) = default;

  //     TensorBasis(TensorBasis&&) = default;

  //     void operator=(const TensorBasis&) = delete;

  //     TensorBasis& operator=(TensorBasis&&) = default;

  //     inline
  //     constexpr
  //     size_t getDimension() const
  //     {
  //       return m_data.rows();
  //     }

  //     inline
  //     constexpr
  //     size_t getDOFs() const
  //     {
  //       return m_data.cols();
  //     }

  //     inline
  //     auto operator()(size_t i) const
  //     {
  //       return m_data.col(i);
  //     }

  //     friend inline TensorBasis operator+(const TensorBasis&, const TensorBasis&);
  //     friend inline TensorBasis operator-(const TensorBasis&, const TensorBasis&);
  //     friend inline TensorBasis operator-(const TensorBasis&);
  //     friend inline TensorBasis operator*(Scalar, const TensorBasis& rhs);
  //     friend inline TensorBasis operator*(const TensorBasis&, Scalar);
  //     friend inline TensorBasis operator/(const TensorBasis&, Scalar);

  //     const Math::Matrix& getMatrix() const
  //     {
  //       return m_data;
  //     }

  //   private:
  //     Math::Matrix m_data;
  // };

  // inline
  // TensorBasis<Math::Vector> operator+(const TensorBasis<Math::Vector>& lhs, const TensorBasis<Math::Vector>& rhs)
  // {
  //   assert(lhs.getDOFs() == rhs.getDOFs());
  //   return lhs.m_data + rhs.m_data;
  // }

  // inline
  // TensorBasis<Math::Vector> operator-(const TensorBasis<Math::Vector>& lhs, const TensorBasis<Math::Vector>& rhs)
  // {
  //   assert(lhs.getDOFs() == rhs.getDOFs());
  //   return lhs.m_data - rhs.m_data;
  // }

  // inline
  // TensorBasis<Math::Vector> operator-(const TensorBasis<Math::Vector>& op)
  // {
  //   return -op.m_data;
  // }

  // inline
  // TensorBasis<Math::Vector> operator*(Scalar lhs, const TensorBasis<Math::Vector>& rhs)
  // {
  //   return lhs * rhs.m_data;
  // }

  // inline
  // TensorBasis<Math::Vector> operator*(const TensorBasis<Math::Vector>& lhs, Scalar rhs)
  // {
  //   return lhs.m_data * rhs;
  // }

  // inline
  // TensorBasis<Math::Vector> operator/(const TensorBasis<Math::Vector>& lhs, Scalar rhs)
  // {
  //   return lhs.m_data / rhs;
  // }

  // template <>
  // class TensorBasis<Math::Matrix> final
  // {
  //   public:
  //     using ValueType = Math::Matrix;

  //     TensorBasis() = default;

  //     template <class ... Args>
  //     constexpr
  //     TensorBasis(Args&&... args)
  //       : m_data(std::forward<Args>(args)...)
  //     {}

  //     template <class F, typename = std::enable_if_t<std::is_invocable_v<F, size_t>>>
  //     constexpr
  //     TensorBasis(size_t dofs, size_t rows, size_t cols, F&& f)
  //       : m_data(rows, cols, dofs)
  //     {
  //       for (size_t i = 0; i < dofs; i++)
  //         operator()(i) = f(i);
  //     }

  //     TensorBasis(const TensorBasis& other) = default;

  //     TensorBasis(TensorBasis&& other) = default;

  //     TensorBasis& operator=(const TensorBasis&) = delete;

  //     TensorBasis& operator=(TensorBasis&&) = default;

  //     inline
  //     size_t getRows() const
  //     {
  //       return m_data.dimension(0);
  //     }

  //     inline
  //     size_t getColumns() const
  //     {
  //       return m_data.dimension(1);
  //     }

  //     inline
  //     size_t getDOFs() const
  //     {
  //       return m_data.dimension(2);
  //     }

  //     inline
  //     Eigen::Map<const Math::Matrix, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>
  //     operator()(size_t i) const
  //     {
  //       return {
  //         m_data.data() + i * m_data.dimension(0) * m_data.dimension(1),
  //         m_data.dimension(0), m_data.dimension(1),
  //         { m_data.dimension(0), 1 }
  //       };
  //     }

  //     friend inline TensorBasis operator+(const TensorBasis&, const TensorBasis&);
  //     friend inline TensorBasis operator-(const TensorBasis&, const TensorBasis&);
  //     friend inline TensorBasis operator-(const TensorBasis&);
  //     friend inline TensorBasis operator*(Scalar, const TensorBasis& rhs);
  //     friend inline TensorBasis operator*(const TensorBasis&, Scalar);
  //     friend inline TensorBasis operator/(const TensorBasis&, Scalar);

  //     const Math::Tensor<3>& getTensor() const
  //     {
  //       return m_data;
  //     }

  //   private:
  //     Math::Tensor<3> m_data;
  // };

  // inline
  // TensorBasis<Math::Matrix>
  // operator+(const TensorBasis<Math::Matrix>& lhs, const TensorBasis<Math::Matrix>& rhs)
  // {
  //   assert(lhs.getDOFs() == rhs.getDOFs());
  //   return lhs.m_data + rhs.m_data;
  // }

  // inline
  // TensorBasis<Math::Matrix>
  // operator-(const TensorBasis<Math::Matrix>& lhs, const TensorBasis<Math::Matrix>& rhs)
  // {
  //   assert(lhs.getDOFs() == rhs.getDOFs());
  //   return lhs.m_data - rhs.m_data;
  // }

  // inline
  // TensorBasis<Math::Matrix>
  // operator-(const TensorBasis<Math::Matrix>& op)
  // {
  //   return -op.m_data;
  // }

  // inline
  // TensorBasis<Math::Matrix>
  // operator*(Scalar lhs, const TensorBasis<Math::Matrix>& rhs)
  // {
  //   return lhs * rhs.m_data;
  // }

  // inline
  // TensorBasis<Math::Matrix>
  // operator*(const TensorBasis<Math::Matrix>& lhs, Scalar rhs)
  // {
  //   return lhs.m_data * rhs;
  // }

  // inline
  // TensorBasis<Math::Matrix>
  // operator/(const TensorBasis<Math::Matrix>& lhs, Scalar rhs)
  // {
  //   return lhs.m_data / rhs;
  // }
}

#endif
