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

      size_t nqp; ///< Number of quadrature points.
      size_t ndof; ///< Number of scalar basis functions.
      size_t dim; ///< Reference dimension, including zero.
      std::vector<Scalar> phi; ///< Values indexed by qp*ndof+a.
      std::vector<Scalar> dphi; ///< Derivatives indexed by (qp*ndof+a)*dim+i.

      /// @brief Constructs an empty reference table.
      ReferenceTabulation()
        : nqp(0),
          ndof(0),
          dim(0),
          phi(),
          dphi()
      {}

      /**
       * @brief Gets a reference basis value.
       * @param qp Quadrature-point index.
       * @param a Local basis-function index.
       * @returns The tabulated basis value.
       */
      const Scalar& getBasis(size_t qp, size_t a) const noexcept
      {
        return phi[qp * ndof + a];
      }

      /**
       * @brief Gets a tabulated reference derivative.
       * @tparam Order Zero selects values; one selects first derivatives.
       * Higher orders retain the H1 table's historical zero result.
       * @param qp Quadrature-point index.
       * @param a Local basis-function index.
       * @param i Reference coordinate, zero when selecting a value.
       * @returns The selected value, derivative, or historical higher-order zero.
       */
      template <size_t Order>
      const Scalar& getDerivative(size_t qp, size_t a, size_t i) const noexcept
      {
        if constexpr (Order == 0)
        {
          assert(i == 0);
          return getBasis(qp, a);
        }
        else if constexpr (Order == 1)
        {
          assert(i < dim);
          return dphi[(qp * ndof + a) * dim + i];
        }
        else
        {
          static const Scalar zero = Scalar(0);
          return zero;
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
        if (dim == 0)
          return {};
        return {dphi.data() + (qp * ndof + a) * dim, dim};
      }
  };
}
#endif
