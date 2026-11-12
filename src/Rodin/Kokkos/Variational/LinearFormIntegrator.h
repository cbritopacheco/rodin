/*
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_VARIATIONAL_LINEARFORMINTEGRATOR_H
#define RODIN_KOKKOS_VARIATIONAL_LINEARFORMINTEGRATOR_H

/**
 * @file
 * @brief Backend interface for linear contribution preparation and evaluation.
 */
#include "Integration.h"

namespace Rodin::Kokkos::Variational
{
  /**
   * @brief Interface of linear form integrators evaluated on the device.
   *
   * Assembly::Kokkos discovers this interface on the integrators of a
   * linear form and asks each one for its contributions as device entries.
   * An integration is a call to prepare(), which locates the entries of the
   * cells, followed by fill(), which writes them. Preparation updates owned
   * selection and layout state through a non-const integrator; fill() reads
   * that prepared state without changing it.
   *
   * @tparam Backend Kokkos backend of the mesh context.
   */
  template <class Backend>
  class LinearFormIntegrator
  {
    public:
      /// @brief Mesh context.
      using ContextType = ::Rodin::Context::Kokkos<::Rodin::Context::Local, Backend>;

      /// @brief Memory space of the entries.
      using MemorySpace = typename ContextType::MemorySpaceType;

      /// @brief Device entry type.
      using EntriesType = Entries<MemorySpace>;

      virtual ~LinearFormIntegrator() = default;

      /// @brief Gets the arrangement of the entries. @return Entry layout.
      virtual Layout getLayout() const = 0;

      /**
       * @brief Locates the entries of every selected cell.
       * @return Number of entries.
       */
      virtual size_t prepare() = 0;

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
       * @param indices Entry indices, test degrees of freedom, written as the
       * values unless empty.
       *
       * The entries are the ones located by the last call to prepare().
       */
      virtual void fill(const ::Kokkos::View<Real*, MemorySpace>& values, size_t offset,
        Real sign, bool accumulate,
        const ::Kokkos::View<Index*, MemorySpace>& indices) const = 0;

      /**
       * @brief Integrates the form over its region into new entries.
       * @return Local contributions, indexed by test degrees of freedom and
       * arranged as by fill() at offset zero.
       */
      EntriesType assemble()
      {
        const size_t count = prepare();
        EntriesType result;
        result.indices = ::Kokkos::View<Index*, MemorySpace>(
          ::Kokkos::view_alloc(
            ::Kokkos::WithoutInitializing, std::string("Rodin entry indices")),
          count);
        result.values = ::Kokkos::View<Real*, MemorySpace>(
          ::Kokkos::view_alloc(
            ::Kokkos::WithoutInitializing, std::string("Rodin entry values")),
          count);
        fill(result.values, 0, Real(1), false, result.indices);
        ::Kokkos::fence("Rodin linear form integration");
        return result;
      }
  };

}

#endif
