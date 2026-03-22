/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_IO_HDF5_H
#define RODIN_PETSC_IO_HDF5_H

/**
 * @file HDF5.h
 * @brief HDF5 IO support for PETSc-backed grid functions.
 *
 * Provides `GridFunctionPrinter` and `GridFunctionLoader` specializations
 * that serialise the locally-owned portion of a PETSc `Vec` to/from HDF5
 * files.  This enables checkpoint/restart workflows and post-processing
 * with external HDF5-compatible tools.
 *
 * ## HDF5 File Layout
 *
 * ```
 * /GridFunction/
 * ├── Meta/
 * │   ├── Size       (unsigned long long — number of local DOFs)
 * │   └── Dimension  (unsigned long long — vector dimension)
 * └── Values/
 *     └── Data       (1-D array of double — DOF coefficients)
 * ```
 *
 * @see Rodin::IO::GridFunctionPrinter,
 *      Rodin::IO::GridFunctionLoader,
 *      Rodin::PETSc::Variational::GridFunction
 */

#include <petscvec.h>
#include <petscsys.h>
#include <petscmath.h>
#include <vector>
#include <numeric>

#include "Rodin/Alert/Notation.h"
#include "Rodin/IO/HDF5.h"
#include "Rodin/Alert.h"

#include <hdf5.h>

namespace Rodin::IO
{
  /**
   * @brief Internal helpers for PETSc HDF5 serialization primitives.
   */
  namespace Internal
  {
    /**
     * @brief Writes a single `unsigned long long` scalar to an HDF5 dataset.
     * @param file  Open HDF5 file handle.
     * @param path  Dataset path within the file (e.g. `"/GridFunction/Meta/Size"`).
     * @param value The value to write.
     */
    inline void writeScalarULL(hid_t file, const char* path, unsigned long long value)
    {
      const auto space = H5Screate(H5S_SCALAR);
      assert(space >= 0);
      const auto ds = H5Dcreate2(file, path, H5T_NATIVE_ULLONG, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      assert(ds >= 0);
      const auto writeStatus = H5Dwrite(ds, H5T_NATIVE_ULLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);
      assert(writeStatus >= 0);
      (void) writeStatus;
      H5Dclose(ds);
      H5Sclose(space);
    }

    /**
     * @brief Reads a 1-D array of `double` values from an HDF5 dataset.
     * @param file Open HDF5 file handle.
     * @param path Dataset path within the file.
     * @returns Vector of `double` values read from the dataset.
     */
    inline std::vector<double> readVectorDouble(hid_t file, const char* path)
    {
      const auto ds = H5Dopen2(file, path, H5P_DEFAULT);
      assert(ds >= 0);
      const auto space = H5Dget_space(ds);
      assert(space >= 0);
      hsize_t dims[1] = {0};
      const auto ndims = H5Sget_simple_extent_ndims(space);
      assert(ndims == 1);
      (void) ndims;
      const auto dimStatus = H5Sget_simple_extent_dims(space, dims, nullptr);
      assert(dimStatus >= 0);
      (void) dimStatus;
      std::vector<double> v(static_cast<size_t>(dims[0]));
      if (!v.empty())
      {
        const auto readStatus = H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data());
        assert(readStatus >= 0);
        (void) readStatus;
      }
      H5Sclose(space);
      H5Dclose(ds);
      return v;
    }

    /**
     * @brief Reads a single `unsigned long long` scalar from an HDF5 dataset.
     * @param file Open HDF5 file handle.
     * @param path Dataset path within the file.
     * @returns The scalar value read.
     */
    inline unsigned long long readScalarULL(hid_t file, const char* path)
    {
      const auto ds = H5Dopen2(file, path, H5P_DEFAULT);
      assert(ds >= 0);
      unsigned long long value = 0;
      const auto readStatus = H5Dread(ds, H5T_NATIVE_ULLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);
      assert(readStatus >= 0);
      (void) readStatus;
      H5Dclose(ds);
      return value;
    }

  }

  /**
   * @brief HDF5 loader for PETSc-backed grid functions.
   *
   * Reads per-rank DOF values from an HDF5 file into the local portion
   * of a distributed PETSc vector.
   *
   * ## HDF5 Layout
   *
   * The expected file structure is:
   * - `/GridFunction/Values/Data` — 1-D array of DOF values
   * - `/GridFunction/Meta/Size` — number of DOFs
   * - `/GridFunction/Meta/Dimension` — vector dimension
   *
   * @tparam FES Finite element space type.
   *
   * @note Stream-based loading is not supported; use the path-based overload.
   *
   * @see GridFunctionPrinter
   */
  template <class FES>
  class GridFunctionLoader<FileFormat::HDF5, FES, ::Vec>
    : public GridFunctionLoaderBase<FES, ::Vec>
  {
    public:
      /// @brief PETSc vector data type (`::Vec`).
      using DataType = ::Vec;
      /// @brief Grid function type being loaded.
      using ObjectType = Variational::GridFunction<FES, DataType>;
      /// @brief Parent loader base class.
      using Parent = GridFunctionLoaderBase<FES, DataType>;

      /**
       * @brief Constructs an HDF5 loader bound to a PETSc-backed grid function.
       */
      GridFunctionLoader(ObjectType& gf)
        : Parent(gf)
      {}

      /**
       * @brief Stream-based overload — not supported for HDF5.
       * @throws Alert::MemberFunctionException Always; use load(path) instead.
       */
      void load(std::istream&) override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "HDF5 GridFunction loading is file-path based."
          << "Please use the "
          << Alert::Identifier::Function("load(const boost::filesystem::path&)")
          << " overload."
          << Alert::Raise;
      }

      /**
       * @brief Loads grid function data from an HDF5 file.
       * @param filename Path to the HDF5 file.
       */
      void load(const boost::filesystem::path& filename) override
      {
        auto& gf = this->getObject();
        auto& vec = gf.getData();

        PetscInt rb = 0, re = 0;
        auto ierr = VecGetOwnershipRange(vec, &rb, &re);
        assert(ierr == PETSC_SUCCESS);
        const PetscInt localN = re - rb;

        const auto file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        if (file < 0)
        {
          Alert::Exception()
            << "Failed to open HDF5 file for reading: " << filename
            << Alert::Raise;
        }

        const auto values = Internal::readVectorDouble(file, "/GridFunction/Values/Data");

        if (static_cast<PetscInt>(values.size()) != localN)
        {
          H5Fclose(file);
          Alert::Exception()
            << "Local DOF count mismatch: file has " << values.size()
            << ", local Vec portion has " << localN << "."
            << Alert::Raise;
        }

        H5Fclose(file);

        std::vector<PetscInt> indices(static_cast<size_t>(localN));
        std::iota(indices.begin(), indices.end(), rb);
        std::vector<PetscScalar> localValues(static_cast<size_t>(localN));
        for (PetscInt i = 0; i < localN; ++i)
          localValues[static_cast<size_t>(i)] = static_cast<PetscScalar>(values[static_cast<size_t>(i)]);

        ierr = VecSetValues(vec, localN, indices.data(), localValues.data(), INSERT_VALUES);
        assert(ierr == PETSC_SUCCESS);

        ierr = VecAssemblyBegin(vec);
        assert(ierr == PETSC_SUCCESS);

        ierr = VecAssemblyEnd(vec);
        assert(ierr == PETSC_SUCCESS);
      }
  };

  /**
   * @brief HDF5 printer for PETSc-backed grid functions.
   *
   * Writes the local portion of a distributed PETSc vector to an HDF5
   * file using the same layout expected by @ref GridFunctionLoader.
   *
   * @tparam FES Finite element space type.
   *
   * @note Stream-based printing is not supported; use the path-based overload.
   *
   * @see GridFunctionLoader
   */
  template <class FES>
  class GridFunctionPrinter<FileFormat::HDF5, FES, ::Vec> final
    : public GridFunctionPrinterBase<FileFormat::HDF5, FES, ::Vec>
  {
    public:
      /// @brief PETSc vector data type (`::Vec`).
      using DataType = ::Vec;
      /// @brief Grid function type being printed.
      using ObjectType = Variational::GridFunction<FES, DataType>;
      /// @brief Parent printer base class.
      using Parent = GridFunctionPrinterBase<FileFormat::HDF5, FES, DataType>;

      /**
       * @brief Constructs an HDF5 printer bound to a PETSc-backed grid function.
       */
      GridFunctionPrinter(const ObjectType& gf)
        : Parent(gf)
      {}

      /**
       * @brief Stream-based overload — not supported for HDF5.
       * @throws Alert::MemberFunctionException Always; use print(path) instead.
       */
      void print(std::ostream&) override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "HDF5 GridFunction printing is file-path based."
          << "Please use the "
          << Alert::Identifier::Function("print(const boost::filesystem::path&)")
          << " overload."
          << Alert::Raise;
      }

      /**
       * @brief Writes grid function data to an HDF5 file.
       * @param filename Path for the output HDF5 file.
       */
      void print(const boost::filesystem::path& filename) override
      {
        const auto& gf = this->getObject();
        const auto& vec = gf.getData();

        PetscInt rb = 0, re = 0;
        auto ierr = VecGetOwnershipRange(vec, &rb, &re);
        assert(ierr == PETSC_SUCCESS);
        const PetscInt localN = re - rb;

        const PetscScalar* raw = PETSC_NULLPTR;
        ierr = VecGetArrayRead(vec, &raw);
        assert(ierr == PETSC_SUCCESS);

        std::vector<double> values(static_cast<size_t>(localN));
        for (PetscInt i = 0; i < localN; ++i)
          values[static_cast<size_t>(i)] = static_cast<double>(PetscRealPart(raw[i]));

        ierr = VecRestoreArrayRead(vec, &raw);
        assert(ierr == PETSC_SUCCESS);

        const auto file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        if (file < 0)
        {
          Alert::Exception()
            << "Failed to create HDF5 file: " << filename
            << Alert::Raise;
        }

        const auto gfGroup = H5Gcreate2(file, "/GridFunction", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        assert(gfGroup >= 0);
        H5Gclose(gfGroup);

        const auto metaGroup = H5Gcreate2(file, "/GridFunction/Meta", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        assert(metaGroup >= 0);
        H5Gclose(metaGroup);

        const auto valuesGroup = H5Gcreate2(file, "/GridFunction/Values", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        assert(valuesGroup >= 0);
        H5Gclose(valuesGroup);

        const hsize_t dims[1] = { static_cast<hsize_t>(values.size()) };
        const auto dataSpace = H5Screate_simple(1, dims, nullptr);
        assert(dataSpace >= 0);

        const auto dataSet = H5Dcreate2(
            file, "/GridFunction/Values/Data", H5T_NATIVE_DOUBLE, dataSpace,
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        assert(dataSet >= 0);

        if (!values.empty())
        {
          const auto writeStatus = H5Dwrite(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                          H5P_DEFAULT, values.data());
          assert(writeStatus >= 0);
          (void) writeStatus;
        }
        H5Dclose(dataSet);
        H5Sclose(dataSpace);

        Internal::writeScalarULL(file, "/GridFunction/Meta/Size",
            static_cast<unsigned long long>(localN));
        Internal::writeScalarULL(file, "/GridFunction/Meta/Dimension",
            static_cast<unsigned long long>(gf.getDimension()));

        H5Fclose(file);
      }

  };
}

#endif
