/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

#include <gtest/gtest.h>
#include <boost/filesystem.hpp>

#include <petsc.h>

#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/IO.h>
#include <Rodin/PETSc/Variational/GridFunction.h>
#include <Rodin/PETSc/IO.h>

#include <hdf5.h>

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  // Helper to create a mesh for each polytope type
  Geometry::Mesh<Context::Local> makeMesh(Polytope::Type type)
  {
    using LocalMesh = Geometry::Mesh<Context::Local>;
    switch (type)
    {
      case Polytope::Type::Segment:
        return LocalMesh::UniformGrid(type, { 11 });
      case Polytope::Type::Triangle:
      case Polytope::Type::Quadrilateral:
        return LocalMesh::UniformGrid(type, { 4, 4 });
      case Polytope::Type::Tetrahedron:
      case Polytope::Type::Hexahedron:
      case Polytope::Type::Wedge:
        return LocalMesh::UniformGrid(type, { 3, 3, 3 });
      default:
        ADD_FAILURE() << "Unsupported polytope type for makeMesh";
        return Geometry::Mesh<Context::Local>();
    }
  }

  std::string polytopeLabel(Polytope::Type type)
  {
    switch (type)
    {
      case Polytope::Type::Segment:       return "Segment1D";
      case Polytope::Type::Triangle:      return "Triangle2D";
      case Polytope::Type::Quadrilateral: return "Quad2D";
      case Polytope::Type::Tetrahedron:   return "Tet3D";
      case Polytope::Type::Hexahedron:    return "Hex3D";
      case Polytope::Type::Wedge:         return "Wedge3D";
      default:                            return "Unknown";
    }
  }

  // Parameterized test fixture with PETSc initialization
  class PETScHDF5 : public ::testing::TestWithParam<Polytope::Type>
  {
    protected:
      static void SetUpTestSuite()
      {
        PetscInitialize(nullptr, nullptr, nullptr, nullptr);
      }

      static void TearDownTestSuite()
      {
        PetscFinalize();
      }
  };

  // --- PETSc GridFunction HDF5 standalone field file -------------------------

  TEST_P(PETScHDF5, GridFunctionStandalone)
  {
    const auto type = GetParam();
    const std::string gfFile =
        "/tmp/rodin_petsc_hdf5_gf_sa_" + polytopeLabel(type) + ".h5";

    auto mesh = makeMesh(type);
    P1 fes(mesh);

    // Create a PETSc-backed grid function using sequential communicator
    Rodin::PETSc::Variational::GridFunction gf(fes);

    // Set values
    PetscScalar* array = nullptr;
    VecGetArray(gf.getData(), &array);
    PetscInt n = 0;
    VecGetLocalSize(gf.getData(), &n);
    for (PetscInt i = 0; i < n; ++i)
      array[i] = static_cast<PetscScalar>(1.0 + static_cast<double>(i) * 0.1);
    VecRestoreArray(gf.getData(), &array);

    // Save
    gf.save(gfFile, FileFormat::HDF5);

    // Verify HDF5 file structure
    hid_t h5 = H5Fopen(gfFile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    ASSERT_GE(h5, 0);

    hid_t metaSize = H5Dopen2(h5, "/GridFunction/Meta/Size", H5P_DEFAULT);
    ASSERT_GE(metaSize, 0);
    H5Dclose(metaSize);

    hid_t metaDim = H5Dopen2(h5, "/GridFunction/Meta/Dimension", H5P_DEFAULT);
    ASSERT_GE(metaDim, 0);
    H5Dclose(metaDim);

    hid_t values = H5Dopen2(h5, "/GridFunction/Values/Data", H5P_DEFAULT);
    ASSERT_GE(values, 0);
    hid_t dspace = H5Dget_space(values);
    int rank = H5Sget_simple_extent_ndims(dspace);
    ASSERT_EQ(rank, 1);
    hsize_t count = 0;
    H5Sget_simple_extent_dims(dspace, &count, nullptr);
    EXPECT_EQ(static_cast<size_t>(count), gf.getSize());

    H5Sclose(dspace);
    H5Dclose(values);
    H5Fclose(h5);
    std::remove(gfFile.c_str());
  }

  // --- PETSc GridFunction HDF5 round-trip raw DOFs ---------------------------

  TEST_P(PETScHDF5, GridFunctionPersistenceRawDOFs)
  {
    const auto type = GetParam();
    const std::string gfFile =
        "/tmp/rodin_petsc_hdf5_gf_dofs_" + polytopeLabel(type) + ".h5";

    auto mesh = makeMesh(type);
    P1 fes(mesh);

    Rodin::PETSc::Variational::GridFunction gf(fes);

    // Fill with known values
    PetscScalar* array = nullptr;
    VecGetArray(gf.getData(), &array);
    PetscInt n = 0;
    VecGetLocalSize(gf.getData(), &n);
    for (PetscInt i = 0; i < n; ++i)
      array[i] = static_cast<PetscScalar>(std::sin(static_cast<double>(i)));
    VecRestoreArray(gf.getData(), &array);

    gf.save(gfFile, FileFormat::HDF5);

    // Read back raw values and compare
    hid_t h5 = H5Fopen(gfFile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    ASSERT_GE(h5, 0);
    hid_t dset = H5Dopen2(h5, "/GridFunction/Values/Data", H5P_DEFAULT);
    ASSERT_GE(dset, 0);
    hid_t dspace = H5Dget_space(dset);
    hsize_t count = 0;
    H5Sget_simple_extent_dims(dspace, &count, nullptr);
    EXPECT_EQ(static_cast<PetscInt>(count), n);

    std::vector<double> saved(static_cast<size_t>(count));
    H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, saved.data());

    for (PetscInt i = 0; i < n; ++i)
    {
      EXPECT_DOUBLE_EQ(saved[static_cast<size_t>(i)],
                        static_cast<double>(std::sin(static_cast<double>(i))))
          << "DOF mismatch at index " << i;
    }

    H5Sclose(dspace);
    H5Dclose(dset);
    H5Fclose(h5);

    // Round-trip: load into a new grid function
    Rodin::PETSc::Variational::GridFunction loaded(fes);
    loaded.load(gfFile, FileFormat::HDF5);

    // Compare sizes
    PetscInt loadedN = 0;
    VecGetSize(loaded.getData(), &loadedN);
    EXPECT_EQ(loadedN, n);

    // Compare values
    const PetscScalar* origArray = nullptr;
    const PetscScalar* loadedArray = nullptr;
    VecGetArrayRead(gf.getData(), &origArray);
    VecGetArrayRead(loaded.getData(), &loadedArray);
    PetscInt localN = 0;
    VecGetLocalSize(gf.getData(), &localN);
    for (PetscInt i = 0; i < localN; ++i)
    {
      EXPECT_DOUBLE_EQ(
          static_cast<double>(PetscRealPart(loadedArray[i])),
          static_cast<double>(PetscRealPart(origArray[i])))
          << "Round-trip DOF mismatch at index " << i;
    }
    VecRestoreArrayRead(gf.getData(), &origArray);
    VecRestoreArrayRead(loaded.getData(), &loadedArray);

    std::remove(gfFile.c_str());
  }

  struct PETScPolytopeNameGenerator
  {
    std::string operator()(const ::testing::TestParamInfo<Polytope::Type>& info) const
    {
      return polytopeLabel(info.param);
    }
  };

  INSTANTIATE_TEST_SUITE_P(
      AllDimensions,
      PETScHDF5,
      ::testing::Values(
          Polytope::Type::Segment,
          Polytope::Type::Triangle,
          Polytope::Type::Quadrilateral,
          Polytope::Type::Tetrahedron,
          Polytope::Type::Hexahedron,
          Polytope::Type::Wedge
      ),
      PETScPolytopeNameGenerator());
}
