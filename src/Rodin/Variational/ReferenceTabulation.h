/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_REFERENCETABULATION_H
#define RODIN_VARIATIONAL_REFERENCETABULATION_H

/**
 * @file
 * @brief Numerical reference-basis values and derivatives, independent of caching.
 */

#include <cassert>
#include <span>
#include <limits>
#include <utility>
#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Alert/Raise.h"
#include <vector>
#include "Rodin/Types.h"

namespace Rodin::Variational
{
  /**
   * @brief Scalar basis values and first derivatives at reference samples.
   *
   * Values are indexed by sample and basis; derivatives additionally carry
   * a reference-coordinate index. This data object owns its storage and has
   * no cache identity, lookup, eviction or finite-element dependency.
   * @tparam Scalar Stored scalar type.
   */
  template <class Scalar>
  class ReferenceTabulation final
  {
    public:
      /// @brief Writable values of all scalar basis functions at one reference point.
      struct Values
      {
        std::span<Scalar> data; ///< Basis-indexed output storage.

        /**
         * @brief Borrows output storage.
         * @param output Basis-indexed values.
         */
        Values(std::span<Scalar> output)
          : data(output)
        {}
      };

      /// @brief Writable reference gradients, indexed by basis then coordinate.
      struct Gradients
      {
        std::span<Scalar> data; ///< Basis-major output storage.
        size_t dimension; ///< Number of reference coordinates per basis function.

        /**
         * @brief Borrows basis-major derivative storage.
         * @param output Derivatives indexed by basis then coordinate.
         * @param dim Number of reference coordinates.
         */
        Gradients(std::span<Scalar> output, size_t dim)
          : data(output),
            dimension(dim)
        {}
      };

      /// @brief Constructs an empty reference table.
      ReferenceTabulation()
        : m_size(0),
          m_count(0),
          m_dimension(0),
          m_values(),
          m_derivatives()
      {}

      /**
       * @brief Copies reference dimensions and numerical storage.
       * @param other Source table.
       */
      ReferenceTabulation(const ReferenceTabulation& other) = default;

      /**
       * @brief Transfers storage and leaves an empty source table.
       * @param other Source table.
       */
      ReferenceTabulation(ReferenceTabulation&& other) noexcept
        : m_size(std::exchange(other.m_size, 0)),
          m_count(std::exchange(other.m_count, 0)),
          m_dimension(std::exchange(other.m_dimension, 0)),
          m_values(std::move(other.m_values)),
          m_derivatives(std::move(other.m_derivatives))
      {
        other.m_values.clear();
        other.m_derivatives.clear();
      }

      /**
       * @brief Copies all extents and storage without publishing a partial table.
       * @param other Source table.
       * @returns This table.
       */
      ReferenceTabulation& operator=(const ReferenceTabulation& other)
      {
        if (this != &other)
        {
          ReferenceTabulation copy(other);
          *this = std::move(copy);
        }
        return *this;
      }

      /**
       * @brief Transfers all storage and leaves an empty source table.
       * @param other Source table.
       * @returns This table.
       */
      ReferenceTabulation& operator=(ReferenceTabulation&& other) noexcept
      {
        if (this != &other)
        {
          m_size = std::exchange(other.m_size, 0);
          m_count = std::exchange(other.m_count, 0);
          m_dimension = std::exchange(other.m_dimension, 0);
          m_values = std::move(other.m_values);
          m_derivatives = std::move(other.m_derivatives);
          other.m_values.clear();
          other.m_derivatives.clear();
        }
        return *this;
      }

      /**
       * @brief Sizes values and first derivatives together.
       *
       * Resizing can invalidate previously borrowed views, including when a
       * later allocation fails. Storage sizes and published extents remain coherent.
       * @param samples Number of reference samples.
       * @param count Number of scalar basis functions.
       * @param dimension Number of reference coordinates, including zero.
       */
      void resize(size_t samples, size_t count, size_t dimension)
      {
        constexpr size_t Maximum = std::numeric_limits<size_t>::max();
        if ((count && samples > Maximum / count) ||
          (dimension && samples * count > Maximum / dimension))
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Reference tabulation extents overflow the storage size." << Alert::Raise;
        }
        const size_t values = samples * count;
        const size_t derivatives = values * dimension;
        // Reserve both arrays before changing extents, retaining reusable capacity.
        m_values.reserve(values);
        m_derivatives.reserve(derivatives);
        m_values.resize(values);
        m_derivatives.resize(derivatives);
        m_size = samples;
        m_count = count;
        m_dimension = dimension;
      }

      /**
       * @brief Gets the number of reference samples.
       * @returns The sample count.
       */
      size_t getSize() const noexcept
      {
        return m_size;
      }

      /**
       * @brief Gets the number of scalar basis functions.
       * @returns The basis count.
       */
      size_t getCount() const noexcept
      {
        return m_count;
      }

      /**
       * @brief Gets the reference dimension.
       * @returns The reference dimension, including zero.
       */
      size_t getDimension() const noexcept
      {
        return m_dimension;
      }

      /**
       * @brief Gets the owned values in sample-major order.
       * @returns Read-only basis storage.
       */
      const std::vector<Scalar>& getValues() const noexcept
      {
        return m_values;
      }

      /**
       * @brief Gets the owned derivatives in sample/basis/coordinate order.
       * @returns Read-only first-derivative storage.
       */
      const std::vector<Scalar>& getDerivatives() const noexcept
      {
        return m_derivatives;
      }

      /**
       * @brief Borrows writable basis values at one sample.
       * @param sample Reference sample index.
       * @returns Basis-indexed output with the established basis count.
       */
      Values getValues(size_t sample) noexcept
      {
        assert(sample < m_size);
        return Values(std::span<Scalar>(m_values).subspan(sample * m_count, m_count));
      }

      /**
       * @brief Borrows writable first derivatives at one sample.
       * @param sample Reference sample index.
       * @returns Basis-major output with the established reference dimension.
       */
      Gradients getGradients(size_t sample) noexcept
      {
        assert(sample < m_size);
        return Gradients(
          std::span<Scalar>(m_derivatives)
            .subspan(sample * m_count * m_dimension, m_count * m_dimension),
          m_dimension);
      }

      /**
       * @brief Gets a reference basis value.
       * @param qp Quadrature-point index.
       * @param a Local basis-function index.
       * @returns The tabulated basis value.
       */
      const Scalar& getBasis(size_t qp, size_t a) const noexcept
      {
        return m_values[qp * m_count + a];
      }

      /**
       * @brief Gets a tabulated reference derivative.
       * @tparam Order Zero selects values; one selects first derivatives.
       * @param qp Quadrature-point index.
       * @param a Local basis-function index.
       * @param i Reference coordinate, zero when selecting a value.
       * @returns The selected value or first derivative.
       */
      template <size_t Order>
        requires(Order <= 1)
      const Scalar& getDerivative(size_t qp, size_t a, size_t i) const noexcept
      {
        if constexpr (Order == 0)
        {
          assert(i == 0);
          return getBasis(qp, a);
        }
        else if constexpr (Order == 1)
        {
          assert(i < m_dimension);
          return m_derivatives[(qp * m_count + a) * m_dimension + i];
        }
      }

      /**
       * @brief Borrows the reference gradient of one basis function.
       * @param qp Quadrature-point index.
       * @param a Local basis-function index.
       * @returns A span of reference dimension, empty on a point.
       */
      std::span<const Scalar> getGradient(size_t qp, size_t a) const noexcept
      {
        if (m_dimension == 0)
          return {};
        return {m_derivatives.data() + (qp * m_count + a) * m_dimension, m_dimension};
      }

    private:
      size_t m_size;
      size_t m_count;
      size_t m_dimension;
      std::vector<Scalar> m_values;
      std::vector<Scalar> m_derivatives;
  };
}
#endif
