/*
 *          Copyright Carlos BRITO PACHECO 2026.
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_VARIATIONAL_GRIDFUNCTIONDATA_H
#define RODIN_KOKKOS_VARIATIONAL_GRIDFUNCTIONDATA_H

/**
 * @file
 * @brief Prepared access to the Data container of a grid function.
 */

#include <string>
#include <type_traits>
#include <vector>
#include <Kokkos_Core.hpp>
#include "Rodin/FormLanguage/Traits.h"
#include "Rodin/Kokkos/Math/Vector.h"
#include "Rodin/Kokkos/Utility.h"

namespace Rodin::Kokkos::Variational
{

  /// @brief Coefficient vector matching a Kokkos finite element space.
  template <class FES>
  using Data =
    ::Rodin::Kokkos::Math::Vector<typename ::Rodin::FormLanguage::Traits<FES>::ScalarType,
      typename FES::MeshType::ContextType::BackendType>;

  /**
   * @brief Copies the coefficients of a grid function to the memory space
   * of its mesh, in the local numbering of the device degrees of freedom.
   *
   * This primary template reads containers indexed by coeff(), such as
   * ::Rodin::Math::Vector. Modules whose grid functions keep their coefficients in
   * another container specialize it, as Rodin/PETSc/Kokkos/Variational/GridFunctionData.h does
   * for PETSc vectors.
   *
   * @tparam Data Coefficient container of the grid function.
   */
  template <class Data>
  struct GridFunctionData
  {
    /**
     * @param fes Device space of the grid function.
     * @param data Coefficients, numbered as by the space.
     * @returns Backend data ordered on the mesh execution instance.
     */
    template <class FES>
    static ::Kokkos::View<Real*, typename FES::MeshType::MemorySpace> prepare(
      const FES& fes, const Data& data)
    {
      const auto& mesh = fes.getMesh();
      const auto& execution = mesh.getContext().getExecutionSpace();
      const Index size = data.size();
      if constexpr (requires { mesh.getTransferCache(); })
      {
        // Lowered at every assembly, so that the same coefficients stay
        // resident instead of being reallocated and copied each time.
        std::vector<Real> values(size);
        for (Index i = 0; i < size; ++i)
          values[i] = data.coeff(i);
        return mesh.getTransferCache().upload(
          execution, std::move(values), "Rodin grid function");
      }
      else
      {
        ::Kokkos::View<Real*, typename FES::MeshType::MemorySpace> res(
          ::Kokkos::view_alloc(execution, std::string("Rodin grid function")), size);
        auto host = ::Kokkos::create_mirror_view(res);
        for (Index i = 0; i < size; ++i)
          host(i) = data.coeff(i);
        ::Kokkos::deep_copy(execution, res, host);
        execution.fence("Rodin grid function transfer");
        return res;
      }
    }
  };

  /// @brief Borrows coefficients already stored in the mesh memory space.
  template <class Scalar, class Backend>
  struct GridFunctionData<::Rodin::Kokkos::Math::Vector<Scalar, Backend>>
  {
    /// @param fes Owning space. @param data Backend coefficients. @return Borrowed coefficients, without a transfer.
    template <class FES>
    static auto prepare([[maybe_unused]] const FES& fes,
      const ::Rodin::Kokkos::Math::Vector<Scalar, Backend>& data)
    {
      static_assert(std::is_same_v<typename FES::MeshType::MemorySpace,
        typename ::Rodin::Kokkos::Math::Vector<Scalar, Backend>::MemorySpace>);
      return data.getView();
    }
  };
}

#endif
