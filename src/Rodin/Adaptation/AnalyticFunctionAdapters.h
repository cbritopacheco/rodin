/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_ANALYTICFUNCTIONADAPTERS_H
#define RODIN_ADAPTATION_ANALYTICFUNCTIONADAPTERS_H

/**
 * @file
 * @brief Lightweight `FunctionBase` adapters for analytic vector- and
 *        matrix-valued functions defined by a single callable.
 *
 * `Rodin::Variational::RealFunction(F)` already wraps a scalar callable
 * with signature `(const Geometry::Point&) -> Real`. The stock
 * `VectorFunction` and `MatrixFunction` only accept either a constant
 * value or a component-wise variadic list of scalar callables. For
 * analytic level sets the natural shape is a SINGLE callable returning
 * the whole vector or matrix; this header adds the missing adapters.
 *
 * Usage:
 *
 *   Adaptation::AnalyticVectorFunction grad(
 *     [&](const Geometry::Point& p) -> Math::SpatialVector<Real>
 *     {
 *       const auto& c = p.getPhysicalCoordinates();
 *       return analyticLevelSet.grad(c);
 *     });
 *
 *   Adaptation::AnalyticMatrixFunction hess(
 *     [&](const Geometry::Point& p) -> Math::SpatialMatrix<Real>
 *     {
 *       const auto& c = p.getPhysicalCoordinates();
 *       return analyticLevelSet.hess(c);
 *     });
 */

#include "Rodin/Geometry/Point.h"
#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Math/SpatialMatrix.h"
#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Types.h"
#include "Rodin/Variational/Function.h"
#include "Rodin/Variational/MatrixFunction.h"
#include "Rodin/Variational/VectorFunction.h"

namespace Rodin::Adaptation
{
  /**
   * @brief Vector-valued `FunctionBase` adapter built from a single callable
   *        of signature `(const Geometry::Point&) -> Math::SpatialVector<Real>`.
   */
  template <class F>
  class AnalyticVectorFunction final
    : public Variational::VectorFunctionBase<Real, AnalyticVectorFunction<F>>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = Real;
      /// @brief Small spatial vector value type.
      using SpatialVectorType = Math::SpatialVector<ScalarType>;
      /// @brief Range (evaluation value) type.
      using RangeType = SpatialVectorType;
      /// @brief Parent class type.
      using Parent =
        Variational::VectorFunctionBase<ScalarType, AnalyticVectorFunction<F>>;

      /**
       * @brief Constructs the adapter from a callable and vector dimension.
       * @param f Function operand.
       * @param dimension Spatial dimension.
       */
      AnalyticVectorFunction(F f, std::size_t dimension)
        : m_f(std::move(f)),
          m_dimension(dimension)
      {}

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      AnalyticVectorFunction(const AnalyticVectorFunction& other)
        : Parent(other),
          m_f(other.m_f),
          m_dimension(other.m_dimension)
      {}

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      AnalyticVectorFunction(AnalyticVectorFunction&& other)
        : Parent(std::move(other)),
          m_f(std::move(other.m_f)),
          m_dimension(other.m_dimension)
      {}

      /**
       * @brief Evaluates the wrapped vector-valued callable.
       * @param p Point at which the operation is evaluated.
       * @returns Value of the expression at the supplied evaluation point.
       */
      RangeType getValue(const Geometry::Point& p) const
      {
        return m_f(p);
      }

      /**
       * @brief Returns the vector dimension.
       * @returns The vector dimension.
       */
      std::size_t getDimension() const noexcept
      {
        return m_dimension;
      }

      /**
       * @brief Returns no intrinsic polynomial order for analytic callables.
       * @returns Polynomial order on the entity, or an empty optional when no order is available.
       */
      Optional<std::size_t> getOrder(const Geometry::Polytope&) const noexcept
      {
        return std::nullopt;
      }

      AnalyticVectorFunction* copy() const noexcept override
      {
        return new AnalyticVectorFunction(*this);
      }

    private:
      F m_f;
      std::size_t m_dimension;
  };

  template <class F>
  /// @brief Deduction guide for analytic vector functions.
  AnalyticVectorFunction(F, std::size_t) -> AnalyticVectorFunction<F>;

  /**
   * @brief Matrix-valued `FunctionBase` adapter built from a single callable
   *        of signature `(const Geometry::Point&) -> Math::SpatialMatrix<Real>`.
   */
  template <class F>
  class AnalyticMatrixFunction final
    : public Variational::MatrixFunctionBase<Real, AnalyticMatrixFunction<F>>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = Real;
      /// @brief Small spatial matrix value type.
      using SpatialMatrixType = Math::SpatialMatrix<ScalarType>;
      /// @brief Range (evaluation value) type.
      using RangeType = SpatialMatrixType;
      /// @brief Parent class type.
      using Parent =
        Variational::MatrixFunctionBase<ScalarType, AnalyticMatrixFunction<F>>;

      /**
       * @brief Constructs the adapter from a callable and matrix dimensions.
       * @param f Function operand.
       * @param rows Number of rows.
       * @param cols Number of columns.
       */
      AnalyticMatrixFunction(F f, std::size_t rows, std::size_t cols)
        : m_f(std::move(f)),
          m_rows(rows),
          m_cols(cols)
      {}

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      AnalyticMatrixFunction(const AnalyticMatrixFunction& other)
        : Parent(other),
          m_f(other.m_f),
          m_rows(other.m_rows),
          m_cols(other.m_cols)
      {}

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      AnalyticMatrixFunction(AnalyticMatrixFunction&& other)
        : Parent(std::move(other)),
          m_f(std::move(other.m_f)),
          m_rows(other.m_rows),
          m_cols(other.m_cols)
      {}

      /**
       * @brief Evaluates the wrapped matrix-valued callable.
       * @param p Point at which the operation is evaluated.
       * @returns Value of the expression at the supplied evaluation point.
       */
      RangeType getValue(const Geometry::Point& p) const
      {
        return m_f(p);
      }

      /**
       * @brief Returns the row count.
       * @returns The row count.
       */
      std::size_t getRows() const noexcept
      {
        return m_rows;
      }
      /**
       * @brief Returns the column count.
       * @returns The column count.
       */
      std::size_t getColumns() const noexcept
      {
        return m_cols;
      }

      /**
       * @brief Returns no intrinsic polynomial order for analytic callables.
       * @returns Polynomial order on the entity, or an empty optional when no order is available.
       */
      Optional<std::size_t> getOrder(const Geometry::Polytope&) const noexcept
      {
        return std::nullopt;
      }

      /**
       * @brief Clones this adapter.
       * @returns Pointer to a newly allocated copy; the caller owns the returned object.
       */
      AnalyticMatrixFunction* copy() const noexcept override
      {
        return new AnalyticMatrixFunction(*this);
      }

    private:
      F m_f;
      std::size_t m_rows;
      std::size_t m_cols;
  };

  template <class F>
  /// @brief Deduction guide for analytic matrix functions.
  AnalyticMatrixFunction(F, std::size_t, std::size_t) -> AnalyticMatrixFunction<F>;
}

#endif
