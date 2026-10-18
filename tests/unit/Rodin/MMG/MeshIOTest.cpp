/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief MMG mesh file-open failures retain topology and required metadata. */
#include <gtest/gtest.h>
#include <boost/filesystem.hpp>
#include <string>
#include "Rodin/MMG/Mesh.h"
#include "Rodin/IO/ForwardDecls.h"
#include "Rodin/Alert/MemberFunctionException.h"

using namespace Rodin;
using namespace Rodin::Geometry;

namespace
{
  class MMGMeshIO : public ::testing::Test
  {
    protected:
      boost::filesystem::path directory;

      void SetUp() override
      {
        directory = boost::filesystem::temp_directory_path() /
          boost::filesystem::unique_path("rodin-mmg-io-%%%%-%%%%");
        boost::filesystem::create_directory(directory);
      }

      void TearDown() override
      {
        boost::filesystem::remove_all(directory);
      }
  };
}

/// @brief Failed MEDIT input/output preserve tags and permit a valid round trip on the same mesh.
TEST_F(MMGMeshIO, FileOpenFailuresPreserveStateAndPermitRecovery)
{
  MMG::Mesh mesh(LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2}));
  mesh.getConnectivity().compute(1, 2);
  mesh.setCorner(0);
  mesh.setRequiredVertex(1);
  const auto cornerCoordinates = mesh.getVertexCoordinates(0);
  const auto requiredCoordinates = mesh.getVertexCoordinates(1);
  const size_t cells = mesh.getCellCount();
  const size_t vertices = mesh.getVertexCount();
  const auto missing = directory / "absent" / "mesh.mesh";
  const auto verifyState = [&] {
    EXPECT_EQ(mesh.getCellCount(), cells);
    EXPECT_EQ(mesh.getVertexCount(), vertices);
    ASSERT_EQ(mesh.getCorners().size(), 1u);
    ASSERT_EQ(mesh.getRequiredVertices().size(), 1u);
    EXPECT_EQ(mesh.getVertexCoordinates(*mesh.getCorners().begin()).getData(),
      cornerCoordinates.getData());
    EXPECT_EQ(mesh.getVertexCoordinates(*mesh.getRequiredVertices().begin()).getData(),
      requiredCoordinates.getData());
  };
  const auto expectOpenFailure = [&](auto&& action, const std::string& mode) {
    try
    {
      action();
      FAIL() << "A path with a nonexistent parent must fail to open.";
    }
    catch (const Alert::Exception& error)
    {
      const std::string message = error.what();
      EXPECT_NE(message.find("In member function"), std::string::npos);
      EXPECT_NE(message.find("Rodin::MMG::Mesh"), std::string::npos);
      EXPECT_NE(message.find("Failed to open"), std::string::npos);
      EXPECT_NE(message.find(missing.string()), std::string::npos);
      EXPECT_NE(message.find(mode), std::string::npos);
    }
    EXPECT_FALSE(boost::filesystem::exists(missing));
    verifyState();
  };
  expectOpenFailure([&] { mesh.save(missing, IO::FileFormat::MEDIT); }, "for writing");
  expectOpenFailure([&] { mesh.load(missing, IO::FileFormat::MEDIT); }, "for reading");
  const auto valid = directory / "valid.mesh";
  ASSERT_NO_THROW(mesh.save(valid, IO::FileFormat::MEDIT));
  ASSERT_TRUE(boost::filesystem::is_regular_file(valid));
  EXPECT_EQ(&mesh.load(valid, IO::FileFormat::MEDIT), &mesh);
  verifyState();
}
