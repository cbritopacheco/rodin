/*
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_VARIATIONAL_BILINEARFORMINTEGRATOR_H
#define RODIN_KOKKOS_VARIATIONAL_BILINEARFORMINTEGRATOR_H

/**
 * @file
 * @brief Backend interface for bilinear contribution preparation and evaluation.
 */
#include "Integration.h"

namespace Rodin::Kokkos::Variational
{
  /**
   * @brief Interface of bilinear form integrators evaluated on the device.
   *
   * Assembly::Kokkos discovers this interface on the local integrators of a
   * bilinear form and asks each one for its contributions as device triplets.
   * An integration is a call to prepare(), which locates the entries of the
   * cells, followed by fill(), which writes them. Preparation updates owned
   * selection and layout state through a non-const integrator; fill() reads
   * that prepared state without changing it.
   *
   * @tparam Backend Kokkos backend of the mesh context.
   */
  template <class Backend>
  class BilinearFormIntegrator
  {
    public:
      /// @brief Mesh context.
      using ContextType = ::Rodin::Context::Kokkos<::Rodin::Context::Local, Backend>;

      /// @brief Memory space of the triplets.
      using MemorySpace = typename ContextType::MemorySpaceType;

      /// @brief Device triplet type.
      using TripletsType = Triplets<MemorySpace>;

      virtual ~BilinearFormIntegrator() = default;

      /// @brief Gets the arrangement of the entries. @return Entry layout.
      virtual Layout getLayout() const = 0;

      /**
       * @brief Tells whether the entries vanish off the diagonal of the
       * component blocks, so that prepare() can list the diagonals only.
       * @return Whether the local matrices are diagonal by component blocks.
       */
      virtual bool isBlockDiagonal() const
      {
        return false;
      }

      /**
       * @brief Locates the entries of every selected cell.
       * @param diagonal Whether to list only the entries on the diagonal of
       * the component blocks, honored when isBlockDiagonal(). Integrators
       * with equal layouts prepared with the same honored choice list the
       * same coordinates in the same order.
       * @return Number of entries.
       */
      virtual size_t prepare(bool diagonal) = 0;

      /**
       * @brief Gets the spacing of the entries located by the last call to
       * prepare().
       * @return Spacing given to getEntryPosition().
       */
      virtual size_t getSpacing() const = 0;

      /**
       * @brief Integrates the form over its region.
       * @param values Entry values, the entry of rank r being written at
       * offset + getEntryPosition(r, count, getSpacing()).
       * @param offset Position of the first entry.
       * @param sign Factor of the contributions.
       * @param accumulate Whether the contributions are added to the values
       * instead of replacing them.
       * @param rows Row indices, test degrees of freedom, written as the
       * values unless empty.
       * @param cols Column indices, trial degrees of freedom, written as the
       * values unless empty.
       *
       * The entries are the ones located by the last call to prepare().
       */
      virtual void fill(const ::Kokkos::View<Real*, MemorySpace>& values, size_t offset,
        Real sign, bool accumulate, const ::Kokkos::View<Index*, MemorySpace>& rows,
        const ::Kokkos::View<Index*, MemorySpace>& cols) const = 0;

      /**
       * @brief Integrates the form over its region into new triplets.
       * @return Local contributions, with rows indexing test and columns
       * indexing trial degrees of freedom, arranged as by fill() at offset
       * zero.
       */
      TripletsType assemble()
      {
        const size_t count = prepare(false);
        TripletsType result;
        result.rows = ::Kokkos::View<Index*, MemorySpace>(
          ::Kokkos::view_alloc(
            ::Kokkos::WithoutInitializing, std::string("Rodin triplet rows")),
          count);
        result.cols = ::Kokkos::View<Index*, MemorySpace>(
          ::Kokkos::view_alloc(
            ::Kokkos::WithoutInitializing, std::string("Rodin triplet columns")),
          count);
        result.values = ::Kokkos::View<Real*, MemorySpace>(
          ::Kokkos::view_alloc(
            ::Kokkos::WithoutInitializing, std::string("Rodin triplet values")),
          count);
        fill(result.values, 0, Real(1), false, result.rows, result.cols);
        ::Kokkos::fence("Rodin bilinear form integration");
        return result;
      }
  };

}

#endif
