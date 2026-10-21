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
#include <vector>

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
        return LocalMesh::UniformGrid(type, {11});
      case Polytope::Type::Triangle:
      case Polytope::Type::Quadrilateral:
        return LocalMesh::UniformGrid(type, {4, 4});
      case Polytope::Type::Tetrahedron:
      case Polytope::Type::Hexahedron:
      case Polytope::Type::Pyramid:
      case Polytope::Type::Wedge:
        return LocalMesh::UniformGrid(type, {3, 3, 3});
      default:
        ADD_FAILURE() << "Unsupported polytope type for makeMesh";
        return Geometry::Mesh<Context::Local>();
    }
  }

  std::string polytopeLabel(Polytope::Type type)
  {
    switch (type)
    {
      case Polytope::Type::Segment:
        return "Segment1D";
      case Polytope::Type::Triangle:
        return "Triangle2D";
      case Polytope::Type::Quadrilateral:
        return "Quad2D";
      case Polytope::Type::Tetrahedron:
        return "Tet3D";
      case Polytope::Type::Hexahedron:
        return "Hex3D";
      case Polytope::Type::Pyramid:
        return "Pyramid3D";
      case Polytope::Type::Wedge:
        return "Wedge3D";
      default:
        return "Unknown";
    }
  }

  // Parameterized test fixture (PETSc is initialized once in main())
  class PETScHDF5 : public ::testing::TestWithParam<Polytope::Type>
  {};

  // --- PETSc GridFunction HDF5 standalone field file -------------------------

  /// @brief Verifies grid function standalone for PET sc HDF 5 by checking exact expected values.
  TEST_P(PETScHDF5, GridFunctionStandalone)
  {
    const auto type = GetParam();
    const std::string gfFile =
      "/tmp/rodin_petsc_hdf5_gf_sa_" + polytopeLabel(type) + ".h5";

    auto mesh = makeMesh(type);
    P1<PetscScalar> fes(mesh);

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

  /// @brief Verifies grid function persistence raw DO fs for PET sc HDF 5 by checking tolerance-based numerical results, exact expected values.
  TEST_P(PETScHDF5, GridFunctionPersistenceRawDOFs)
  {
    const auto type = GetParam();
    const std::string gfFile =
      "/tmp/rodin_petsc_hdf5_gf_dofs_" + polytopeLabel(type) + ".h5";

    auto mesh = makeMesh(type);
    P1<PetscScalar> fes(mesh);

    Rodin::PETSc::Variational::GridFunction gf(fes);

    // Fill with known values
    PetscScalar* array = nullptr;
    VecGetArray(gf.getData(), &array);
    PetscInt n = 0;
    VecGetLocalSize(gf.getData(), &n);
    for (PetscInt i = 0; i < n; ++i)
    {
      array[i] = static_cast<PetscScalar>(std::sin(static_cast<double>(i)));
#if defined(PETSC_USE_COMPLEX)
      array[i] += PetscCMPLX(0, std::cos(static_cast<double>(i)));
#endif
    }
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
#if defined(PETSC_USE_COMPLEX)
    const hid_t imaginarySet =
      H5Dopen2(h5, "/GridFunction/Values/ImaginaryData", H5P_DEFAULT);
    ASSERT_GE(imaginarySet, 0);
    std::vector<double> imaginary(static_cast<size_t>(count));
    ASSERT_GE(H5Dread(imaginarySet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                imaginary.data()),
      0);
    for (PetscInt i = 0; i < n; ++i)
    {
      EXPECT_DOUBLE_EQ(
        imaginary[static_cast<size_t>(i)], std::cos(static_cast<double>(i)))
        << "Imaginary DOF mismatch at index " << i;
    }
    EXPECT_GE(H5Dclose(imaginarySet), 0);
#endif
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
      EXPECT_EQ(loadedArray[i], origArray[i]) << "Round-trip DOF mismatch at index " << i;
    }
    VecRestoreArrayRead(gf.getData(), &origArray);
    VecRestoreArrayRead(loaded.getData(), &loadedArray);

    std::remove(gfFile.c_str());
  }

  /// @brief Isolates malformed checkpoint files and suppresses expected HDF5 diagnostics.
  class PETScHDF5Contracts : public ::testing::Test
  {
    protected:
      void SetUp() override
      {
        m_path = boost::filesystem::temp_directory_path() /
          boost::filesystem::unique_path("rodin-petsc-hdf5-%%%%-%%%%-%%%%.h5");
        ASSERT_GE(H5Eget_auto2(H5E_DEFAULT, &m_errorHandler, &m_errorData), 0);
        ASSERT_GE(H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr), 0);
      }

      void TearDown() override
      {
        std::remove(m_path.c_str());
        EXPECT_GE(H5Eset_auto2(H5E_DEFAULT, m_errorHandler, m_errorData), 0);
      }

      const boost::filesystem::path& path() const
      {
        return m_path;
      }

      void replaceDataset(const char* name, const std::vector<double>& values,
        std::vector<hsize_t> shape = {})
      {
        const auto file = H5Fopen(path().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        ASSERT_GE(file, 0);
        if (H5Lexists(file, name, H5P_DEFAULT) > 0)
        {
          ASSERT_GE(H5Ldelete(file, name, H5P_DEFAULT), 0);
        }
        if (shape.empty())
          shape = {values.size()};
        const auto space = H5Screate_simple(shape.size(), shape.data(), nullptr);
        ASSERT_GE(space, 0);
        const auto dataset = H5Dcreate2(
          file, name, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        ASSERT_GE(dataset, 0);
        ASSERT_GE(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                    values.data()),
          0);
        EXPECT_GE(H5Dclose(dataset), 0);
        EXPECT_GE(H5Sclose(space), 0);
        EXPECT_GE(H5Fclose(file), 0);
      }

      template <class Field>
      void expectCoefficients(const Field& field, bool constant = false)
      {
        PetscInt size = 0;
        ASSERT_EQ(VecGetSize(field.getData(), &size), PETSC_SUCCESS);
        const PetscScalar* data = nullptr;
        ASSERT_EQ(VecGetArrayRead(field.getData(), &data), PETSC_SUCCESS);
        for (PetscInt i = 0; i < size; ++i)
        {
          PetscScalar expected = constant ? PetscScalar(7) : PetscScalar(1 + 2 * i);
#if defined(PETSC_USE_COMPLEX)
          if (!constant)
            expected += PetscCMPLX(0, -0.5 * i);
#endif
          EXPECT_LT(PetscAbsScalar(data[i] - expected), 1e-12) << i;
        }
        ASSERT_EQ(VecRestoreArrayRead(field.getData(), &data), PETSC_SUCCESS);
      }

      template <class Field>
      void fillCoefficients(Field& field)
      {
        PetscInt size = 0;
        ASSERT_EQ(VecGetSize(field.getData(), &size), PETSC_SUCCESS);
        PetscScalar* data = nullptr;
        ASSERT_EQ(VecGetArray(field.getData(), &data), PETSC_SUCCESS);
        for (PetscInt i = 0; i < size; ++i)
        {
          data[i] = PetscScalar(1 + 2 * i);
#if defined(PETSC_USE_COMPLEX)
          data[i] += PetscCMPLX(0, -0.5 * i);
#endif
        }
        ASSERT_EQ(VecRestoreArray(field.getData(), &data), PETSC_SUCCESS);
      }

    private:
      /// @brief Unique checkpoint path owned by this test.
      boost::filesystem::path m_path;

      /// @brief HDF5 error callback restored after the test.
      H5E_auto2_t m_errorHandler = nullptr;

      /// @brief Original callback context.
      void* m_errorData = nullptr;
  };

  /// @brief Matrix checkpoints retain shape and reject equal-size transposed ranges.
  TEST_F(PETScHDF5Contracts, MatrixPersistencePreservesShapeAndCoefficients)
  {
    for (const auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
           Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
           Polytope::Type::Hexahedron, Polytope::Type::Pyramid, Polytope::Type::Wedge})
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      auto mesh = makeMesh(geometry);
      P1<Math::SpatialMatrix<PetscScalar>> space(mesh, 2, 3);
      Rodin::PETSc::Variational::GridFunction field(space), loaded(space);
      fillCoefficients(field);
      field.save(path(), FileFormat::HDF5);
      loaded.load(path(), FileFormat::HDF5);
      expectCoefficients(loaded);
      const auto file = H5Fopen(path().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      ASSERT_GE(file, 0);
      for (const auto& [name, expected] : {std::pair{"/GridFunction/Meta/Rows", 2ULL},
             std::pair{"/GridFunction/Meta/Columns", 3ULL}})
      {
        const auto dataset = H5Dopen2(file, name, H5P_DEFAULT);
        ASSERT_GE(dataset, 0);
        unsigned long long value = 0;
        ASSERT_GE(
          H5Dread(dataset, H5T_NATIVE_ULLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value), 0);
        EXPECT_EQ(value, expected);
        EXPECT_GE(H5Dclose(dataset), 0);
      }
      EXPECT_GE(H5Fclose(file), 0);
      P1<Math::SpatialMatrix<PetscScalar>> transposedSpace(mesh, 3, 2);
      Rodin::PETSc::Variational::GridFunction transposed(transposedSpace);
      ASSERT_EQ(VecSet(transposed.getData(), 7), PETSC_SUCCESS);
      EXPECT_THROW(transposed.load(path(), FileFormat::HDF5), std::exception);
      expectCoefficients(transposed, true);
      expectCoefficients(field);
    }
  }

  /// @brief Invalid coefficient counts are rejected before changing destination values.
  TEST_F(PETScHDF5Contracts, IncompatiblePayloadPreservesDestination)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
    P1<PetscScalar> space(mesh);
    Rodin::PETSc::Variational::GridFunction field(space), loaded(space);
    fillCoefficients(field);
    for (const bool imaginary : {false, true})
    {
      SCOPED_TRACE(imaginary);
      field.save(path(), FileFormat::HDF5);
      replaceDataset(
        imaginary ? "/GridFunction/Values/ImaginaryData" : "/GridFunction/Values/Data",
        {9, 10});
#if defined(PETSC_USE_COMPLEX)
      if (!imaginary)
        replaceDataset("/GridFunction/Values/ImaginaryData", {0, 0});
#endif
      ASSERT_EQ(VecSet(loaded.getData(), 7), PETSC_SUCCESS);
      EXPECT_THROW(loaded.load(path(), FileFormat::HDF5), std::exception);
      expectCoefficients(loaded, true);
      expectCoefficients(field);
    }
#if !defined(PETSC_USE_COMPLEX)
    field.save(path(), FileFormat::HDF5);
    replaceDataset("/GridFunction/Values/ImaginaryData", {0, 1, 0});
    EXPECT_THROW(loaded.load(path(), FileFormat::HDF5), std::exception);
    expectCoefficients(loaded, true);
#endif
  }

  /// @brief Real-only legacy files load into every range with zero imaginary coefficients.
  TEST_F(PETScHDF5Contracts, RealOnlyLegacyPayloadLoadsWithZeroImaginaryPart)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
    const auto check = [&](auto& space) {
      Rodin::PETSc::Variational::GridFunction field(space), loaded(space);
      ASSERT_EQ(VecSet(field.getData(), 7), PETSC_SUCCESS);
      field.save(path(), FileFormat::HDF5);
      const auto file = H5Fopen(path().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
      ASSERT_GE(file, 0);
      if (H5Lexists(file, "/GridFunction/Values/ImaginaryData", H5P_DEFAULT) > 0)
      {
        ASSERT_GE(H5Ldelete(file, "/GridFunction/Values/ImaginaryData", H5P_DEFAULT), 0);
      }
      EXPECT_GE(H5Fclose(file), 0);
      ASSERT_EQ(VecSet(loaded.getData(), -1), PETSC_SUCCESS);
      loaded.load(path(), FileFormat::HDF5);
      expectCoefficients(loaded, true);
    };
    P1<PetscScalar> scalar(mesh);
    P1<Math::SpatialVector<PetscScalar>> vector(mesh, 2);
    P1<Math::SpatialMatrix<PetscScalar>> matrix(mesh, 2, 3);
    check(scalar);
    check(vector);
    check(matrix);
  }

  /// @brief Unsupported streams and failed opens preserve scalar, vector and matrix fields.
  TEST_F(PETScHDF5Contracts, StreamsAndFailedOpensPreserveFields)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
    const auto check = [&](auto& space) {
      Rodin::PETSc::Variational::GridFunction field(space);
      fillCoefficients(field);
      using Space = std::remove_cvref_t<decltype(space)>;
      GridFunctionLoader<FileFormat::HDF5, Space, ::Vec> loader(field);
      GridFunctionPrinter<FileFormat::HDF5, Space, ::Vec> printer(field);
      std::istringstream input("untouched");
      std::ostringstream output;
      EXPECT_THROW(loader.load(input), std::exception);
      EXPECT_THROW(printer.print(output), std::exception);
      EXPECT_EQ(input.tellg(), std::streampos(0));
      EXPECT_TRUE(output.str().empty());
      EXPECT_FALSE(boost::filesystem::exists(path()));
      EXPECT_THROW(field.load(path(), FileFormat::HDF5), std::exception);
      const auto missing =
        boost::filesystem::path(path().string() + ".missing") / "field.h5";
      EXPECT_THROW(field.save(missing, FileFormat::HDF5), std::exception);
      EXPECT_FALSE(boost::filesystem::exists(missing.parent_path()));
      EXPECT_FALSE(boost::filesystem::exists(path()));
      expectCoefficients(field);
    };
    P1<PetscScalar> scalar(mesh);
    P1<Math::SpatialVector<PetscScalar>> vector(mesh, 2);
    P1<Math::SpatialMatrix<PetscScalar>> matrix(mesh, 2, 3);
    check(scalar);
    check(vector);
    check(matrix);
  }

  /// @brief Invalid ranks and nonscalar metadata reject before touching coefficients.
  TEST_F(PETScHDF5Contracts, MalformedDatasetShapesPreserveDestination)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    P1<Math::SpatialMatrix<PetscScalar>> space(mesh, 2, 3);
    GridFunction<decltype(space), ::Vec> field(space), loaded(space);
    fillCoefficients(field);
    ASSERT_EQ(VecSet(loaded.getData(), PetscScalar(7)), PETSC_SUCCESS);
    const auto handles = H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL);
    for (const char* name : {"/GridFunction/Values/Data", "/GridFunction/Meta/Size",
           "/GridFunction/Meta/Dimension", "/GridFunction/Meta/Rows",
           "/GridFunction/Meta/Columns"})
    {
      SCOPED_TRACE(name);
      field.save(path(), FileFormat::HDF5);
      if (std::string(name) == "/GridFunction/Values/Data")
        replaceDataset(name, {1, 2, 3, 4}, {2, 2});
      else
        replaceDataset(name, {1, 2});
      EXPECT_THROW(loaded.load(path(), FileFormat::HDF5), Alert::Exception);
      expectCoefficients(loaded, true);
      EXPECT_EQ(H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL), handles);
    }
  }

  struct PETScPolytopeNameGenerator
  {
    std::string operator()(const ::testing::TestParamInfo<Polytope::Type>& info) const
    {
      return polytopeLabel(info.param);
    }
  };

  INSTANTIATE_TEST_SUITE_P(AllDimensions, PETScHDF5,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
      Polytope::Type::Hexahedron, Polytope::Type::Pyramid, Polytope::Type::Wedge),
    PETScPolytopeNameGenerator());
}

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  ::testing::InitGoogleTest(&argc, argv);
  int result = RUN_ALL_TESTS();
  PetscFinalize();
  return result;
}
