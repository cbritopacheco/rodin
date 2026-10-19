/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_IO_MFEM_H
#define RODIN_PETSC_IO_MFEM_H

/**
 * @file
 * @brief MFEM serialization of PETSc fields through the native serializer.
 */
#include <petscvec.h>
#include <type_traits>

#include "Rodin/IO/MFEM.h"
#include "Rodin/IO/GridFunctionPrinter.h"
#include "Rodin/Math/Vector.h"

namespace Rodin::IO
{
  /**
   * @brief PETSc storage adapter for supported MFEM finite element collections.
   * @tparam FES Finite element space with an MFEM data-printer interface.
   *
   * ## Architecture
   *
   * Select the shard-local space when present, copy its coefficients through
   * the source field's read view, then restore that view before serialization.
   * The native serializer owns node ordering, component layout and changes of
   * basis. The working field borrows its space from the source field, whose
   * lifetime already bounds this printer.
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref GridFunctionPrinter "MFEM with PETSc Vec storage" | Local and MPI P0, P1 and H1 fields with scalar, vector or matrix ranges. |
   */
  template <class FES>
    requires requires(GridFunctionPrinterBase<FileFormat::MFEM, FES, ::Vec>& printer,
      std::ostream& stream) { printer.printData(stream); }
  class GridFunctionPrinter<FileFormat::MFEM, FES, ::Vec> final
    : public GridFunctionPrinterBase<FileFormat::MFEM, FES, ::Vec>
  {
    public:
      /// @brief Finite element space type.
      using FESType = FES;

      /// @brief PETSc coefficient storage type.
      using DataType = ::Vec;

      /// @brief Printed grid function type.
      using ObjectType = Variational::GridFunction<FESType, DataType>;

      /// @brief Format-specific printer interface.
      using Parent = GridFunctionPrinterBase<FileFormat::MFEM, FESType, DataType>;

      /// @brief Scalar field of the finite element space.
      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      /// @brief Context of the source finite element mesh.
      using FESMeshContextType =
        typename FormLanguage::Traits<typename FESType::MeshType>::ContextType;
      using Parent::Parent;

      static_assert(std::is_same_v<ScalarType, PetscScalar>,
        "PETSc MFEM printer: Range scalar must be PetscScalar.");

    private:
      /**
       * @brief Copies local coefficients without assuming a component layout.
       * @returns Working field borrowing the source's local finite element space.
       */
      auto makeLocalField() const
      {
        const auto& field = this->getObject();
        const auto& fes = field.getFiniteElementSpace();
        const auto& local = [&]() -> const auto& {
          if constexpr (requires { fes.getShard(); })
            return fes.getShard();
          else
            return fes;
        }();
        using LocalFES = std::remove_cvref_t<decltype(local)>;
        Variational::GridFunction<LocalFES, Math::Vector<PetscScalar>> result(local);
        field.acquire();
        for (size_t a = 0; a < local.getSize(); ++a)
        {
          const Index global = [&]() {
            if constexpr (requires { fes.getShard(); })
              return fes.getGlobalIndex(a);
            else
              return static_cast<Index>(a);
          }();
          result.getData()[a] = field[global];
        }
        field.flush();
        return result;
      }

    public:
      /**
       * @brief Writes the MFEM collection header and coefficient payload.
       * @param os Output stream receiving the serialized field.
       */
      void print(std::ostream& os) override
      {
        auto localField = makeLocalField();
        using LocalFES =
          std::remove_cvref_t<decltype(localField.getFiniteElementSpace())>;
        GridFunctionPrinter<FileFormat::MFEM, LocalFES, Math::Vector<PetscScalar>>(
          localField)
          .print(os);
      }

      /**
       * @brief Writes the native collection's coefficient payload.
       * @param os Output stream receiving coefficient data.
       */
      void printData(std::ostream& os) override
      {
        auto localField = makeLocalField();
        using LocalFES =
          std::remove_cvref_t<decltype(localField.getFiniteElementSpace())>;
        GridFunctionPrinter<FileFormat::MFEM, LocalFES, Math::Vector<PetscScalar>>(
          localField)
          .printData(os);
      }
  };

  /// @brief Local PETSc-backed H1 MFEM printer.
  template <size_t K, class Range>
  using MFEM_H1_PETSc_Local_Printer = GridFunctionPrinter<FileFormat::MFEM,
    Variational::H1<K, Range, Geometry::Mesh<Context::Local>>, ::Vec>;

  /// @brief MPI PETSc-backed H1 MFEM printer.
  template <size_t K, class Range>
  using MFEM_H1_PETSc_MPI_Printer = GridFunctionPrinter<FileFormat::MFEM,
    Variational::H1<K, Range, Geometry::Mesh<Context::MPI>>, ::Vec>;
}

#endif
