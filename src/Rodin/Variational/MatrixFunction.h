/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_MATRIXFUNCTION_H
#define RODIN_VARIATIONAL_MATRIXFUNCTION_H

#include <set>
#include <optional>

#include "Rodin/Alert.h"

#include "ForwardDecls.h"
#include "Function.h"
#include "RangeShape.h"

namespace Rodin::FormLanguage
{
  template <class Scalar, class Derived>
  struct Traits<Variational::MatrixFunctionBase<Scalar, Derived>>
  {
    using ScalarType = Scalar;
    using DerivedType = Derived;
  };
}

namespace Rodin::Variational
{
  /**
   * @defgroup MatrixFunctionSpecializations MatrixFunction Template Specializations
   * @brief Template specializations of the MatrixFunction class.
   * @see MatrixFunction
   */

  template <class Scalar, class Derived>
  class MatrixFunctionBase : public FunctionBase<MatrixFunctionBase<Scalar, Derived>>
  {
    public:
      using ScalarType = Scalar;

      using MatrixType = Math::Matrix<ScalarType>;

      using Parent = FunctionBase<MatrixFunctionBase<ScalarType, Derived>>;

      MatrixFunctionBase() = default;

      MatrixFunctionBase(const MatrixFunctionBase& other)
        : Parent(other)
      {}

      MatrixFunctionBase(MatrixFunctionBase&& other)
        : Parent(std::move(other))
      {}

      virtual ~MatrixFunctionBase() = default;

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

      inline
      constexpr
      void getValue(MatrixType& res, const Geometry::Point& p) const
      {
        if constexpr (Internal::HasGetValueMethod<Derived, MatrixType&, const Geometry::Point&>::Value)
        {
          return static_cast<const Derived&>(*this).getValue(res, p);
        }
        else
        {
          res = getValue(p);
        }
      }

      inline
      constexpr
      RangeShape getRangeShape() const
      {
        return { getRows(), getColumns() };
      }

      /**
       * @brief Gets the number of rows in the matrix
       * @returns Number of rows
       */
      inline
      constexpr
      size_t getRows() const
      {
        return static_cast<const Derived&>(*this).getRows();
      }

      /**
       * @brief Gets the number of columns in the matrix
       * @returns Number of columns
       */
      inline
      constexpr
      size_t getColumns() const
      {
        return static_cast<const Derived&>(*this).getColumns();
      }

      inline
      constexpr
      MatrixFunctionBase& traceOf(Geometry::Attribute attr)
      {
        Parent::traceOf(attr);
        return *this;
      }

      inline
      constexpr
      MatrixFunctionBase& traceOf(const FlatSet<Geometry::Attribute>& attrs)
      {
        Parent::traceOf(attrs);
        return *this;
      }

      virtual inline MatrixFunctionBase* copy() const noexcept override
      {
        return static_cast<const Derived&>(*this).copy();
      }
  };

  /**
   * @ingroup MatrixFunctionSpecializations
   */
  template <class Scalar>
  class MatrixFunction<Math::Matrix<Scalar>> final
    : public MatrixFunctionBase<Scalar, MatrixFunction<Math::Matrix<Scalar>>>
  {
    public:
      using ScalarType = Scalar;

      using MatrixType = Math::Matrix<ScalarType>;

      using Parent = MatrixFunctionBase<Scalar, MatrixFunction<MatrixType>>;

      MatrixFunction(const MatrixType& matrix)
        : m_matrix(matrix)
      {}

      MatrixFunction(const MatrixFunction& other)
        : Parent(other),
          m_matrix(other.m_matrix)
      {}

      MatrixFunction(MatrixFunction&& other)
        : Parent(std::move(other)),
          m_matrix(std::move(other.m_matrix))
      {}

      inline
      constexpr
      const MatrixType& getValue(const Geometry::Point&) const
      {
        return m_matrix.get();
      }

      inline
      constexpr
      void getValue(MatrixType& res, const Geometry::Point&) const
      {
        res = m_matrix.get();
      }

      inline
      constexpr
      MatrixFunction& traceOf(Geometry::Attribute)
      {
        return *this;
      }

      inline
      constexpr
      MatrixFunction& traceOf(const FlatSet<Geometry::Attribute>& attr)
      {
        return *this;
      }

      inline
      constexpr
      size_t getRows() const
      {
        return m_matrix.get().rows();
      }

      /**
       * @brief Gets the number of columns in the matrix
       * @returns Number of columns
       */
      inline
      constexpr
      size_t getColumns() const
      {
        return m_matrix.get().cols();
      }

      inline MatrixFunction* copy() const noexcept override
      {
        return new MatrixFunction(*this);
      }

    private:
      std::reference_wrapper<const MatrixType> m_matrix;
  };

  MatrixFunction(std::reference_wrapper<const Math::Matrix<Real>>)
    -> MatrixFunction<Math::Matrix<Real>>;
}

#endif
