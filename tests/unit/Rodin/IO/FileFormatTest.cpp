/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <array>
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
