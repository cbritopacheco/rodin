/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <fstream>
#include <gtest/gtest.h>

#include <Rodin/IO.h>
#include <Rodin/IO/MEDIT.h>

#include <Rodin/Geometry.h>
#include <Rodin/Configure.h>

using namespace Rodin::IO;
using namespace Rodin::Geometry;

TEST(Rodin_IO_MeshLoader, SanityTest_MEDIT_2D_Square)
{
  static constexpr const char* filename = "mmg/Square.medit.mesh";
  boost::filesystem::path meshfile;
  meshfile = boost::filesystem::path(RODIN_RESOURCES_DIR);
  meshfile.append(filename);

  boost::filesystem::ifstream in(meshfile);

  Mesh mesh;
  MeshLoader<FileFormat::MEDIT, Rodin::Context::Serial> loader(mesh);
  loader.load(in);

  EXPECT_EQ(mesh.getCount(0), 4);
  EXPECT_EQ(mesh.getCount(0), mesh.getVertexCount());
  EXPECT_EQ(mesh.getAttribute(0, 0), 1);
  EXPECT_EQ(mesh.getAttribute(0, 1), 2);
  EXPECT_EQ(mesh.getAttribute(0, 2), 3);
  EXPECT_EQ(mesh.getAttribute(0, 3), 4);

  EXPECT_EQ(mesh.getCount(1), 5);
  EXPECT_EQ(mesh.getAttribute(1, 0), 1);
  EXPECT_EQ(mesh.getAttribute(1, 1), 2);

  EXPECT_EQ(mesh.getCount(2), 2);
  EXPECT_EQ(mesh.getElementCount(), mesh.getElementCount());
  EXPECT_EQ(mesh.getAttribute(2, 0), 1);
  EXPECT_EQ(mesh.getAttribute(2, 1), 2);
}
