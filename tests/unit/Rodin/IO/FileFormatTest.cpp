/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <array>
#include <fstream>
#include <set>
#include <tuple>
#include <sstream>
#include <string>
#include "Rodin/Geometry.h"
#include "Rodin/IO.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::IO;
using namespace Rodin::Variational;

/**
 * @brief Topology codes match the external XDMF mixed-stream vocabulary.
 *
 * Expected identifiers come from the XDMF Model and Format topology table,
 * https://www.xdmf.org/index.php/XDMF_Model_and_Format#Topology .
 * The file writers must not validate emitted identifiers against the same
 * mapping helper that produced them.
 */
TEST(FileFormatContracts, XDMFTopologyVocabulary)
{
  using Type = Polytope::Type;
  const std::array cases{std::tuple{Type::Point, 1u, 1u, "Polyvertex"},
    std::tuple{Type::Segment, 2u, 34u, "Edge_3"},
    std::tuple{Type::Triangle, 4u, 36u, "Triangle_6"},
    std::tuple{Type::Quadrilateral, 5u, 35u, "Quadrilateral_9"},
    std::tuple{Type::Tetrahedron, 6u, 38u, "Tetrahedron_10"},
    std::tuple{Type::Pyramid, 7u, 39u, "Pyramid_13"},
    std::tuple{Type::Wedge, 8u, 41u, "Wedge_18"},
    std::tuple{Type::Hexahedron, 9u, 50u, "Hexahedron_27"}};
  for (const auto& [geometry, linear, quadratic, name] : cases)
  {
    SCOPED_TRACE(int(geometry));
    EXPECT_EQ(HDF5::getXDMFMixedTopologyId(geometry), linear);
    EXPECT_EQ(HDF5::getXDMFQuadraticMixedTopologyId(geometry), quadratic);
    EXPECT_STREQ(HDF5::getXDMFQuadraticTopologyName(geometry), name);
  }
}

/// @brief MFEM section names are exact and unknown enumeration values have no spelling.
TEST(FileFormatContracts, MFEMKeywordSpellings)
{
  using Keyword = MFEM::Keyword;
  EXPECT_STREQ(MFEM::toCharString(Keyword::Dimension), "dimension");
  EXPECT_STREQ(MFEM::toCharString(Keyword::Elements), "elements");
  EXPECT_STREQ(MFEM::toCharString(Keyword::Boundary), "boundary");
  EXPECT_STREQ(MFEM::toCharString(Keyword::Vertices), "vertices");
  EXPECT_EQ(MFEM::toCharString(static_cast<Keyword>(-1)), nullptr);
}

/// @brief Known formats retain their names; unknown values retain numeric diagnostics and a usable stream.
TEST(FileFormatContracts, KnownAndUnknownValuesStreamWithoutFailure)
{
  std::ostringstream stream;
  stream << FileFormat::MFEM << ' ' << FileFormat::MEDIT << ' ' << FileFormat::HDF5 << ' '
         << static_cast<FileFormat>(-1) << ' ' << static_cast<FileFormat>(255) << " end";
  EXPECT_TRUE(stream.good());
  EXPECT_EQ(stream.str(), "MFEM MEDIT HDF5 -1 255 end");
  EXPECT_EQ(toCharString(static_cast<FileFormat>(255)), nullptr);
}

/// @brief Unsupported runtime formats preserve the field and filesystem while reporting member context and the rejected value.
TEST(FileFormatContracts, GridFunctionRejectsUnsupportedFormatsWithoutSideEffects)
{
  const auto directory = boost::filesystem::temp_directory_path() /
    boost::filesystem::unique_path("rodin-format-%%%%-%%%%-%%%%");
  const auto filename = directory / "field";
  ASSERT_FALSE(boost::filesystem::exists(directory));
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  const auto check = [&](const auto& space) {
    GridFunction field(space);
    for (Index i = 0; i < field.getSize(); ++i)
      field.getData()(i) = 1 + i;
    const auto original = field.getData().eval();
    for (int value : {-1, 255})
    {
      SCOPED_TRACE(value);
      const auto format = static_cast<FileFormat>(value);
      for (bool loading : {true, false})
      {
        SCOPED_TRACE(loading);
        try
        {
          if (loading)
            field.load(filename, format);
          else
            field.save(filename, format);
          FAIL() << "An unsupported format was accepted";
        }
        catch (const Alert::Exception& error)
        {
          const std::string message = error.what();
          EXPECT_NE(
            message.find(loading ? "In member function load" : "In member function save"),
            std::string::npos);
          EXPECT_NE(message.find(loading ? "Unsupported file format for loading"
                                         : "Unsupported file format for saving"),
            std::string::npos);
          EXPECT_NE(
            message.find("Format: \"" + std::to_string(value) + "\""), std::string::npos);
        }
        ASSERT_EQ(field.getData().size(), original.size());
        for (Eigen::Index i = 0; i < original.size(); ++i)
          EXPECT_EQ(field.getData()(i), original(i));
        EXPECT_FALSE(boost::filesystem::exists(directory));
      }
    }
  };
  check(P1(mesh));
  check(P1(mesh, 2));
}

namespace
{
  /// @brief Owns files used to exercise the common stream/file interface.
  class FileInterfaceContracts : public ::testing::Test
  {
    protected:
      void SetUp() override
      {
        directory = boost::filesystem::temp_directory_path() /
          boost::filesystem::unique_path("rodin-file-interface-%%%%-%%%%");
        boost::filesystem::create_directory(directory);
      }

      void TearDown() override
      {
        boost::filesystem::remove_all(directory);
      }

      boost::filesystem::path directory;
  };
}

/// @brief Generic mesh file overloads dispatch to concrete text formats and truncate old output.
TEST_F(FileInterfaceContracts, MeshFilesDispatchAndTruncate)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  mesh.setAttribute({2, 0}, 17);
  mesh.setAttribute({2, 1}, 23);
  const auto check = [&]<FileFormat Format>() {
    SCOPED_TRACE(static_cast<int>(Format));
    const auto path = directory / "mesh";
    const std::string stale(16384, '!');
    {
      std::ofstream old(path.string());
      old << stale;
    }
    MeshPrinter<Format, Context::Local> concretePrinter(mesh);
    Printer<LocalMesh>& printer = concretePrinter;
    EXPECT_EQ(&printer.getObject(), &mesh);
    printer.print(path);
    std::ifstream input(path.string());
    ASSERT_TRUE(input.good());
    const std::string text((std::istreambuf_iterator<char>(input)), {});
    ASSERT_EQ(text.find('!'), std::string::npos);
    EXPECT_NE(
      text.find(Format == FileFormat::MFEM ? "MFEM mesh" : "MeshVersionFormatted"),
      std::string::npos);
    LocalMesh loaded;
    MeshLoader<Format, Context::Local> concreteLoader(loaded);
    Loader<LocalMesh>& loader = concreteLoader;
    loader.load(path);
    ASSERT_EQ(loaded.getDimension(), 2u);
    ASSERT_EQ(loaded.getSpaceDimension(), 2u);
    ASSERT_EQ(loaded.getVertexCount(), 4u);
    ASSERT_EQ(loaded.getCellCount(), 2u);
    std::set<std::pair<Real, Real>> coordinates;
    for (Index i = 0; i < 4; ++i)
    {
      const auto& point = loaded.getVertexCoordinates(i);
      coordinates.emplace(point[0], point[1]);
    }
    const std::set<std::pair<Real, Real>> expected{{0, 0}, {1, 0}, {0, 1}, {1, 1}};
    EXPECT_EQ(coordinates, expected);
    std::set<Attribute> attributes;
    for (auto cell = loaded.getCell(); cell; ++cell)
    {
      const auto attribute = cell->getAttribute();
      ASSERT_TRUE(attribute);
      attributes.insert(*attribute);
    }
    EXPECT_EQ(attributes, (std::set<Attribute>{17, 23}));
    EXPECT_THROW(loader.load(directory / "absent"), Alert::Exception);
    EXPECT_EQ(loaded.getVertexCount(), 4u);
    EXPECT_EQ(loaded.getCellCount(), 2u);
  };
  check.template operator()<FileFormat::MFEM>();
  check.template operator()<FileFormat::MEDIT>();
}

/// @brief Generic field file overloads preserve each coefficient across scalar/vector/matrix formats.
TEST_F(FileInterfaceContracts, FieldFilesDispatchAndTruncate)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  const auto check = [&]<FileFormat Format>(const auto& space) {
    using FES = std::decay_t<decltype(space)>;
    GridFunction field(space), loaded(space);
    for (Index i = 0; i < field.getSize(); ++i)
      field.getData()(i) = 3 + i;
    loaded.getData().setConstant(-1);
    const auto path = directory / "field";
    {
      std::ofstream old(path.string());
      old << std::string(16384, '!');
    }
    GridFunctionPrinter<Format, FES, Math::Vector<Real>> concretePrinter(field);
    Printer<decltype(field)>& printer = concretePrinter;
    EXPECT_EQ(&printer.getObject(), &field);
    printer.print(path);
    std::ifstream input(path.string());
    ASSERT_TRUE(input.good());
    const std::string text((std::istreambuf_iterator<char>(input)), {});
    ASSERT_EQ(text.find('!'), std::string::npos);
    GridFunctionLoader<Format, FES, Math::Vector<Real>> concreteLoader(loaded);
    Loader<decltype(loaded)>& loader = concreteLoader;
    loader.load(path);
    ASSERT_EQ(loaded.getSize(), field.getSize());
    for (Index i = 0; i < loaded.getSize(); ++i)
      EXPECT_EQ(loaded.getData()(i), 3 + i);
  };
  const auto checkSpaces = [&]<FileFormat Format>() {
    check.template operator()<Format>(P1(mesh));
    check.template operator()<Format>(P1(mesh, 2));
    check.template operator()<Format>(P1(mesh, 2, 3));
  };
  checkSpaces.template operator()<FileFormat::MFEM>();
  checkSpaces.template operator()<FileFormat::MEDIT>();
}

/// @brief Unsupported MFEM collections preserve fields and streams and retain working HDF5 persistence.
TEST_F(FileInterfaceContracts, GlobalConstantsRejectMFEMWithoutChangingFieldsOrStreams)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  const auto check = [&](const auto& space) {
    using FES = std::decay_t<decltype(space)>;
    using Scalar = typename FormLanguage::Traits<FES>::ScalarType;
    GridFunction field(space);
    for (Index i = 0; i < field.getSize(); ++i)
    {
      field.getData()(i) = Scalar(3 + i);
      if constexpr (std::is_same_v<Scalar, Complex>)
        field.getData()(i) += Scalar(0, Real(i) + 0.5);
    }
    const auto original = field.getData().eval();
    std::istringstream input("FiniteElementSpace\nFiniteElementCollection: H1_2D_P1\n"
                             "VDim: 1\nOrdering: 0\n\n1\n");
    std::ostringstream output;
    output << "preserved prefix";
    GridFunctionLoader<FileFormat::MFEM, FES, Math::Vector<Scalar>> loader(field);
    GridFunctionPrinter<FileFormat::MFEM, FES, Math::Vector<Scalar>> printer(field);
    EXPECT_EQ(&printer.getObject(), &field);
    for (const bool loading : {true, false})
    {
      SCOPED_TRACE(loading);
      try
      {
        if (loading)
          loader.load(input);
        else
          printer.print(output);
        FAIL() << "An unsupported MFEM collection was accepted";
      }
      catch (const Alert::Exception& error)
      {
        const std::string message = error.what();
        EXPECT_NE(
          message.find(loading ? "In member function load" : "In member function print"),
          std::string::npos);
        EXPECT_NE(message.find(loading
                      ? "No GridFunctionLoader specialization"
                      : "No grid-function printer for this finite element collection"),
          std::string::npos);
      }
      EXPECT_EQ((field.getData() - original).norm(), 0);
      EXPECT_TRUE(input.good());
      EXPECT_EQ(input.tellg(), std::streampos(0));
      EXPECT_TRUE(output.good());
      EXPECT_EQ(output.str(), "preserved prefix");
    }
    const auto path = directory / "global-constant.h5";
    field.save(path, FileFormat::HDF5);
    field.getData().setZero();
    field.load(path, FileFormat::HDF5);
    EXPECT_EQ((field.getData() - original).norm(), 0);
  };
  check(P0g(mesh));
  check(P0g<Math::SpatialVector<Real>, LocalMesh>(mesh, 2));
  check(P0g<Complex, LocalMesh>(mesh));
}
