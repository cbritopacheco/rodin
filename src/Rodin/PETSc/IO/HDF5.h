/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_IO_HDF5_H
#define RODIN_PETSC_IO_HDF5_H

#include <petscvec.h>
#include <petscsys.h>
#include <petscmath.h>
#include <mpi.h>
#include <vector>
#include <numeric>
#include <limits>

#include "Rodin/IO/HDF5.h"
#include "Rodin/Alert.h"

#include <hdf5.h>

namespace Rodin::IO
{
  namespace Internal
  {
    struct GatherResources
    {
      VecScatter scatter = PETSC_NULLPTR;
      Vec gathered = PETSC_NULLPTR;

      ~GatherResources()
      {
        if (scatter != PETSC_NULLPTR)
          VecScatterDestroy(&scatter);
        if (gathered != PETSC_NULLPTR)
          VecDestroy(&gathered);
      }
    };

    inline void assertCondition(bool condition, const char* msg)
    {
      if (!condition)
        Alert::Exception() << msg << Alert::Raise;
    }

    inline std::vector<double> gatherPETScVector(const ::Vec& vec)
    {
      GatherResources resources;
      auto ierr = VecScatterCreateToAll(vec, &resources.scatter, &resources.gathered);
      assertCondition(ierr == PETSC_SUCCESS, "Failed to create PETSc scatter for HDF5 export.");
      ierr = VecScatterBegin(resources.scatter, vec, resources.gathered, INSERT_VALUES, SCATTER_FORWARD);
      assertCondition(ierr == PETSC_SUCCESS, "Failed to begin PETSc scatter for HDF5 export.");
      ierr = VecScatterEnd(resources.scatter, vec, resources.gathered, INSERT_VALUES, SCATTER_FORWARD);
      assertCondition(ierr == PETSC_SUCCESS, "Failed to end PETSc scatter for HDF5 export.");

      PetscInt n = 0;
      ierr = VecGetSize(resources.gathered, &n);
      assertCondition(ierr == PETSC_SUCCESS, "Failed to get PETSc gathered vector size.");

      const PetscScalar* raw = PETSC_NULLPTR;
      ierr = VecGetArrayRead(resources.gathered, &raw);
      assertCondition(ierr == PETSC_SUCCESS, "Failed to get PETSc gathered vector array.");

      std::vector<double> values(static_cast<size_t>(n));
      for (PetscInt i = 0; i < n; ++i)
        values[static_cast<size_t>(i)] = static_cast<double>(PetscRealPart(raw[i]));

      ierr = VecRestoreArrayRead(resources.gathered, &raw);
      assertCondition(ierr == PETSC_SUCCESS, "Failed to restore PETSc gathered vector array.");
      return values;
    }

    inline void writeScalarULL(hid_t file, const char* path, unsigned long long value)
    {
      const auto space = H5Screate(H5S_SCALAR);
      assertCondition(space >= 0, "Failed to create HDF5 scalar dataspace.");
      const auto ds = H5Dcreate2(file, path, H5T_NATIVE_ULLONG, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      assertCondition(ds >= 0, "Failed to create HDF5 scalar dataset.");
      assertCondition(H5Dwrite(ds, H5T_NATIVE_ULLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value) >= 0,
            "Failed to write HDF5 scalar dataset.");
      H5Dclose(ds);
      H5Sclose(space);
    }

    inline std::vector<double> readVectorDouble(hid_t file, const char* path)
    {
      const auto ds = H5Dopen2(file, path, H5P_DEFAULT);
      assertCondition(ds >= 0, "Failed to open HDF5 vector dataset.");
      const auto space = H5Dget_space(ds);
      assertCondition(space >= 0, "Failed to get HDF5 dataspace.");
      hsize_t dims[1] = {0};
      assertCondition(H5Sget_simple_extent_ndims(space) == 1, "Invalid HDF5 vector rank.");
      assertCondition(H5Sget_simple_extent_dims(space, dims, nullptr) >= 0, "Failed to get HDF5 vector dimensions.");
      std::vector<double> v(static_cast<size_t>(dims[0]));
      if (!v.empty())
      {
        assertCondition(H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data()) >= 0,
              "Failed to read HDF5 vector dataset.");
      }
      H5Sclose(space);
      H5Dclose(ds);
      return v;
    }

    inline unsigned long long readScalarULL(hid_t file, const char* path)
    {
      const auto ds = H5Dopen2(file, path, H5P_DEFAULT);
      assertCondition(ds >= 0, "Failed to open HDF5 scalar dataset.");
      unsigned long long value = 0;
      assertCondition(H5Dread(ds, H5T_NATIVE_ULLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value) >= 0,
            "Failed to read HDF5 scalar dataset.");
      H5Dclose(ds);
      return value;
    }

    inline void scatterIntoPETScVector(::Vec& vec, const std::vector<double>& values)
    {
      PetscInt n = 0;
      auto ierr = VecGetSize(vec, &n);
      assertCondition(ierr == PETSC_SUCCESS, "Failed to get PETSc vector size.");
      if (values.size() != static_cast<size_t>(n))
      {
        Alert::Exception()
          << "Invalid PETSc GridFunction data size: expected " << n
          << ", got " << values.size() << "."
          << Alert::Raise;
      }

      PetscInt rb = 0, re = 0;
      ierr = VecGetOwnershipRange(vec, &rb, &re);
      assertCondition(ierr == PETSC_SUCCESS, "Failed to get PETSc ownership range.");
      const PetscInt localN = re - rb;

      std::vector<PetscInt> indices(static_cast<size_t>(localN));
      std::iota(indices.begin(), indices.end(), rb);
      std::vector<PetscScalar> localValues(static_cast<size_t>(localN));
      for (PetscInt i = 0; i < localN; ++i)
        localValues[static_cast<size_t>(i)] = static_cast<PetscScalar>(values[static_cast<size_t>(rb + i)]);

      ierr = VecSetValues(vec, localN, indices.data(), localValues.data(), INSERT_VALUES);
      assertCondition(ierr == PETSC_SUCCESS, "Failed to set PETSc vector values from HDF5.");
      ierr = VecAssemblyBegin(vec);
      assertCondition(ierr == PETSC_SUCCESS, "Failed to begin PETSc vector assembly.");
      ierr = VecAssemblyEnd(vec);
      assertCondition(ierr == PETSC_SUCCESS, "Failed to end PETSc vector assembly.");
    }
  }

  template <class FES>
  class GridFunctionLoader<FileFormat::HDF5, FES, ::Vec>
    : public GridFunctionLoaderBase<FES, ::Vec>
  {
    public:
      using DataType = ::Vec;
      using ObjectType = Variational::GridFunction<FES, DataType>;
      using Parent = GridFunctionLoaderBase<FES, DataType>;

      GridFunctionLoader(ObjectType& gf)
        : Parent(gf)
      {}

      void load(std::istream&) override
      {
        Alert::Exception()
          << "HDF5 GridFunction loading requires file-path based loading."
          << Alert::Raise;
      }

      void load(const boost::filesystem::path& filename) override
      {
        auto& gf = this->getObject();
        auto& vec = gf.getData();

        MPI_Comm comm = MPI_COMM_NULL;
        PetscErrorCode ierr = PetscObjectGetComm(reinterpret_cast<PetscObject>(vec), &comm);
        Internal::assertCondition(ierr == PETSC_SUCCESS, "Failed to get PETSc communicator for HDF5 load.");

        int rank = 0;
        const auto mpiErr = MPI_Comm_rank(comm, &rank);
        Internal::assertCondition(mpiErr == MPI_SUCCESS, "Failed to get MPI rank for HDF5 load.");

        std::vector<double> values;
        if (rank == 0)
        {
          const auto file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
          Internal::assertCondition(file >= 0, "Failed to open HDF5 GridFunction file.");
          values = Internal::readVectorDouble(file, "/GridFunction/Values/Data");
          const auto dofCount = static_cast<size_t>(Internal::readScalarULL(file, "/GridFunction/Meta/Size"));
          const auto vectorDim = static_cast<size_t>(Internal::readScalarULL(file, "/GridFunction/Meta/Dimension"));
          Internal::assertCondition(dofCount == gf.getSize(), "GridFunction size mismatch.");
          Internal::assertCondition(vectorDim == gf.getDimension(), "GridFunction dimension mismatch.");
          H5Fclose(file);
        }

        PetscInt globalN = 0;
        ierr = VecGetSize(vec, &globalN);
        Internal::assertCondition(ierr == PETSC_SUCCESS, "Failed to get PETSc vector size.");
        if (static_cast<size_t>(globalN) > static_cast<size_t>(std::numeric_limits<int>::max()))
        {
          Alert::Exception()
            << "GridFunction vector too large for MPI_Bcast count: size is " << globalN
            << ", maximum is " << std::numeric_limits<int>::max() << "."
            << Alert::Raise;
        }
        if (rank != 0)
          values.resize(static_cast<size_t>(globalN));
        const auto bcastErr = MPI_Bcast(values.data(), static_cast<int>(values.size()), MPI_DOUBLE, 0, comm);
        Internal::assertCondition(bcastErr == MPI_SUCCESS, "Failed to broadcast HDF5 GridFunction data.");

        Internal::scatterIntoPETScVector(vec, values);
      }
  };

  template <class FES>
  class GridFunctionPrinter<FileFormat::HDF5, FES, ::Vec> final
    : public GridFunctionPrinterBase<FileFormat::HDF5, FES, ::Vec>
  {
    public:
      using DataType = ::Vec;
      using ObjectType = Variational::GridFunction<FES, DataType>;
      using Parent = GridFunctionPrinterBase<FileFormat::HDF5, FES, DataType>;

      GridFunctionPrinter(const ObjectType& gf)
        : Parent(gf)
      {}

      void print(std::ostream&) const
      {
        Alert::Exception()
          << "HDF5 GridFunction printing requires file-path based printing."
          << Alert::Raise;
      }

      void print(const boost::filesystem::path& filename) const
      {
        const auto& gf = this->getObject();
        const auto& vec = gf.getData();

        MPI_Comm comm = MPI_COMM_NULL;
        PetscErrorCode ierr = PetscObjectGetComm(reinterpret_cast<PetscObject>(vec), &comm);
        Internal::assertCondition(ierr == PETSC_SUCCESS, "Failed to get PETSc communicator for HDF5 save.");

        int rank = 0;
        const auto mpiErr = MPI_Comm_rank(comm, &rank);
        Internal::assertCondition(mpiErr == MPI_SUCCESS, "Failed to get MPI rank for HDF5 save.");

        const auto values = Internal::gatherPETScVector(vec);
        if (rank == 0)
        {
          const auto file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
          Internal::assertCondition(file >= 0, "Failed to create HDF5 GridFunction file.");

          const auto gfGroup = H5Gcreate2(file, "/GridFunction", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
          Internal::assertCondition(gfGroup >= 0, "Failed to create /GridFunction group.");
          H5Gclose(gfGroup);
          const auto metaGroup = H5Gcreate2(file, "/GridFunction/Meta", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
          Internal::assertCondition(metaGroup >= 0, "Failed to create /GridFunction/Meta group.");
          H5Gclose(metaGroup);
          const auto valuesGroup = H5Gcreate2(file, "/GridFunction/Values", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
          Internal::assertCondition(valuesGroup >= 0, "Failed to create /GridFunction/Values group.");
          H5Gclose(valuesGroup);

          const hsize_t dims[1] = { static_cast<hsize_t>(values.size()) };
          const auto dataSpace = H5Screate_simple(1, dims, nullptr);
          Internal::assertCondition(dataSpace >= 0, "Failed to create HDF5 values dataspace.");
          const auto dataSet = H5Dcreate2(file, "/GridFunction/Values/Data", H5T_NATIVE_DOUBLE, dataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
          Internal::assertCondition(dataSet >= 0, "Failed to create /GridFunction/Values/Data dataset.");
          if (!values.empty())
          {
            Internal::assertCondition(H5Dwrite(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()) >= 0,
                            "Failed to write GridFunction values dataset.");
          }
          H5Dclose(dataSet);
          H5Sclose(dataSpace);

          Internal::writeScalarULL(file, "/GridFunction/Meta/Size",
            static_cast<unsigned long long>(gf.getSize()));
          Internal::writeScalarULL(file, "/GridFunction/Meta/Dimension",
            static_cast<unsigned long long>(gf.getDimension()));

          H5Fclose(file);
        }
      }

  };
}

#endif
