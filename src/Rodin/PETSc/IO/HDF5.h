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

#include "Rodin/Alert/Notation.h"
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

    inline std::vector<double> gatherPETScVector(const ::Vec& vec)
    {
      GatherResources resources;
      auto ierr = VecScatterCreateToAll(vec, &resources.scatter, &resources.gathered);
      assert(ierr == PETSC_SUCCESS);

      ierr = VecScatterBegin(resources.scatter, vec, resources.gathered, INSERT_VALUES, SCATTER_FORWARD);
      assert(ierr == PETSC_SUCCESS);

      ierr = VecScatterEnd(resources.scatter, vec, resources.gathered, INSERT_VALUES, SCATTER_FORWARD);
      assert(ierr == PETSC_SUCCESS);

      PetscInt n = 0;
      ierr = VecGetSize(resources.gathered, &n);
      assert(ierr == PETSC_SUCCESS);

      const PetscScalar* raw = PETSC_NULLPTR;
      ierr = VecGetArrayRead(resources.gathered, &raw);
      assert(ierr == PETSC_SUCCESS);

      std::vector<double> values(static_cast<size_t>(n));
      for (PetscInt i = 0; i < n; ++i)
        values[static_cast<size_t>(i)] = static_cast<double>(PetscRealPart(raw[i]));

      ierr = VecRestoreArrayRead(resources.gathered, &raw);
      assert(ierr == PETSC_SUCCESS);
      return values;
    }

    inline void writeScalarULL(hid_t file, const char* path, unsigned long long value)
    {
      const auto space = H5Screate(H5S_SCALAR);
      assert(space >= 0);
      const auto ds = H5Dcreate2(file, path, H5T_NATIVE_ULLONG, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      assert(ds >= 0);
      assert(H5Dwrite(ds, H5T_NATIVE_ULLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value) >= 0);
      H5Dclose(ds);
      H5Sclose(space);
    }

    inline std::vector<double> readVectorDouble(hid_t file, const char* path)
    {
      const auto ds = H5Dopen2(file, path, H5P_DEFAULT);
      assert(ds >= 0);
      const auto space = H5Dget_space(ds);
      assert(space >= 0);
      hsize_t dims[1] = {0};
      assert(H5Sget_simple_extent_ndims(space) == 1);
      assert(H5Sget_simple_extent_dims(space, dims, nullptr) >= 0);
      std::vector<double> v(static_cast<size_t>(dims[0]));
      if (!v.empty())
      {
        assert(H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data()) >= 0);
      }
      H5Sclose(space);
      H5Dclose(ds);
      return v;
    }

    inline unsigned long long readScalarULL(hid_t file, const char* path)
    {
      const auto ds = H5Dopen2(file, path, H5P_DEFAULT);
      assert(ds >= 0);
      unsigned long long value = 0;
      assert(H5Dread(ds, H5T_NATIVE_ULLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value) >= 0);
      H5Dclose(ds);
      return value;
    }

    inline void scatterIntoPETScVector(::Vec& vec, const std::vector<double>& values)
    {
      PetscInt n = 0;
      auto ierr = VecGetSize(vec, &n);
      assert(ierr == PETSC_SUCCESS);
      if (values.size() != static_cast<size_t>(n))
      {
        Alert::Exception()
          << "Invalid PETSc GridFunction data size: expected " << n
          << ", got " << values.size() << "."
          << Alert::Raise;
      }

      PetscInt rb = 0, re = 0;
      ierr = VecGetOwnershipRange(vec, &rb, &re);
      assert(ierr == PETSC_SUCCESS);
      const PetscInt localN = re - rb;

      std::vector<PetscInt> indices(static_cast<size_t>(localN));
      std::iota(indices.begin(), indices.end(), rb);
      std::vector<PetscScalar> localValues(static_cast<size_t>(localN));
      for (PetscInt i = 0; i < localN; ++i)
        localValues[static_cast<size_t>(i)] = static_cast<PetscScalar>(values[static_cast<size_t>(rb + i)]);

      ierr = VecSetValues(vec, localN, indices.data(), localValues.data(), INSERT_VALUES);
      assert(ierr == PETSC_SUCCESS);

      ierr = VecAssemblyBegin(vec);
      assert(ierr == PETSC_SUCCESS);

      ierr = VecAssemblyEnd(vec);
      assert(ierr == PETSC_SUCCESS);
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
        Alert::MemberFunctionException(*this, __func__)
          << "HDF5 GridFunction loading is file-path based."
          << "Please use the "
          << Alert::Identifier::Function("load(const boost::filesystem::path&)")
          << " overload."
          << Alert::Raise;
      }

      void load(const boost::filesystem::path& filename) override
      {
        // auto& gf = this->getObject();
        // auto& vec = gf.getData();

        // MPI_Comm comm = MPI_COMM_NULL;
        // PetscErrorCode ierr = PetscObjectGetComm(reinterpret_cast<PetscObject>(vec), &comm);
        // assert(ierr == PETSC_SUCCESS);

        // int rank = 0;
        // const auto mpiErr = MPI_Comm_rank(comm, &rank);
        // assert(mpiErr == MPI_SUCCESS);

        // std::vector<double> values;
        // if (rank == 0)
        // {
        //   const auto file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        //   assert(file >= 0);
        //   values = Internal::readVectorDouble(file, "/GridFunction/Values/Data");
        //   const auto dofCount = static_cast<size_t>(Internal::readScalarULL(file, "/GridFunction/Meta/Size"));
        //   const auto vectorDim = static_cast<size_t>(Internal::readScalarULL(file, "/GridFunction/Meta/Dimension"));
        //   assert(dofCount == gf.getSize());
        //   assert(vectorDim == gf.getDimension());
        //   H5Fclose(file);
        // }

        // PetscInt globalN = 0;
        // ierr = VecGetSize(vec, &globalN);
        // assert(ierr == PETSC_SUCCESS);
        // if (static_cast<size_t>(globalN) > static_cast<size_t>(std::numeric_limits<int>::max()))
        // {
        //   Alert::Exception()
        //     << "GridFunction vector too large for MPI_Bcast count: size is " << globalN
        //     << ", maximum is " << std::numeric_limits<int>::max() << "."
        //     << Alert::Raise;
        // }
        // if (rank != 0)
        //   values.resize(static_cast<size_t>(globalN));
        // const auto bcastErr = MPI_Bcast(values.data(), static_cast<int>(values.size()), MPI_DOUBLE, 0, comm);
        // assert(bcastErr == MPI_SUCCESS);

        // Internal::scatterIntoPETScVector(vec, values);
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

      void print(std::ostream&) override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "HDF5 GridFunction printing is file-path based."
          << "Please use the "
          << Alert::Identifier::Function("print(const boost::filesystem::path&)")
          << " overload."
          << Alert::Raise;
      }

      void print(const boost::filesystem::path& filename)
      {
        // const auto& gf = this->getObject();
        // const auto& vec = gf.getData();

        // MPI_Comm comm = MPI_COMM_NULL;
        // PetscErrorCode ierr = PetscObjectGetComm(reinterpret_cast<PetscObject>(vec), &comm);
        // assert(ierr == PETSC_SUCCESS);

        // int rank = 0;
        // const auto mpiErr = MPI_Comm_rank(comm, &rank);
        // assert(mpiErr == MPI_SUCCESS);

        // const auto values = Internal::gatherPETScVector(vec);
        // if (rank == 0)
        // {
        //   const auto file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        //   assert(file >= 0);

        //   const auto gfGroup = H5Gcreate2(file, "/GridFunction", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        //   assert(gfGroup >= 0);

        //   H5Gclose(gfGroup);
        //   const auto metaGroup = H5Gcreate2(file, "/GridFunction/Meta", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        //   assert(metaGroup >= 0);

        //   H5Gclose(metaGroup);
        //   const auto valuesGroup = H5Gcreate2(file, "/GridFunction/Values", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        //   assert(valuesGroup >= 0);

        //   H5Gclose(valuesGroup);

        //   const hsize_t dims[1] = { static_cast<hsize_t>(values.size()) };
        //   const auto dataSpace = H5Screate_simple(1, dims, nullptr);
        //   assert(dataSpace >= 0);

        //   const auto dataSet = H5Dcreate2(file, "/GridFunction/Values/Data", H5T_NATIVE_DOUBLE, dataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        //   assert(dataSet >= 0);

        //   if (!values.empty())
        //   {
        //     assert(H5Dwrite(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()) >= 0);
        //   }
        //   H5Dclose(dataSet);
        //   H5Sclose(dataSpace);

        //   Internal::writeScalarULL(file, "/GridFunction/Meta/Size",
        //     static_cast<unsigned long long>(gf.getSize()));
        //   Internal::writeScalarULL(file, "/GridFunction/Meta/Dimension",
        //     static_cast<unsigned long long>(gf.getDimension()));

        //   H5Fclose(file);
        // }
      }

  };
}

#endif
