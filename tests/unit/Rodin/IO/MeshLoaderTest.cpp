/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <fstream>
#include <gtest/gtest.h>

#include <Rodin/IO.h>
#include <Rodin/IO/MFEM.h>
#include <Rodin/IO/MEDIT.h>

#include <Rodin/Geometry.h>
#include <Rodin/Configure.h>

using namespace Rodin::IO;
using namespace Rodin::Geometry;

TEST(Rodin_IO_MeshLoader, SanityTest_MEDIT_2D_Square)
{
  static constexpr const char* filename = "mmg/Square.medit.mesh";
  static constexpr size_t sdim = 2;
  boost::filesystem::path meshfile;
  meshfile = boost::filesystem::path(RODIN_RESOURCES_DIR);
  meshfile.append(filename);

  boost::filesystem::ifstream in(meshfile);

  Mesh mesh;
  MeshLoader<FileFormat::MEDIT, Rodin::Context::Serial> loader(mesh);
  loader.load(in);

  EXPECT_EQ(mesh.getSpaceDimension(), sdim);

  size_t d = 0;
  EXPECT_EQ(mesh.getCount(d), 4);
  EXPECT_EQ(mesh.getCount(d), mesh.getVertexCount());
  EXPECT_EQ(mesh.getAttribute(d, 0), RODIN_DEFAULT_POLYTOPE_ATTRIBUTE);
  EXPECT_EQ(mesh.getAttribute(d, 1), RODIN_DEFAULT_POLYTOPE_ATTRIBUTE);
  EXPECT_EQ(mesh.getAttribute(d, 2), RODIN_DEFAULT_POLYTOPE_ATTRIBUTE);
  EXPECT_EQ(mesh.getAttribute(d, 3), RODIN_DEFAULT_POLYTOPE_ATTRIBUTE);

  d = 1;
  EXPECT_EQ(mesh.getCount(d), 5);
  EXPECT_EQ(mesh.getAttribute(d, 0), 1);
  EXPECT_EQ(mesh.getAttribute(d, 1), 2);

  d = 2;
  EXPECT_EQ(mesh.getCount(d), 2);
  EXPECT_EQ(mesh.getCount(d), mesh.getElementCount());
  EXPECT_EQ(mesh.getAttribute(d, 0), 1);
  EXPECT_EQ(mesh.getAttribute(d, 1), 2);
}

TEST(Rodin_IO_MeshLoader, SanityTest_MFEM_2D_Square)
{
  static constexpr const char* filename = "mfem/Square.mfem.mesh";
  boost::filesystem::path meshfile;
  meshfile = boost::filesystem::path(RODIN_RESOURCES_DIR);
  meshfile.append(filename);

  boost::filesystem::ifstream in(meshfile);

  Mesh mesh;
  MeshLoader<FileFormat::MFEM, Rodin::Context::Serial> loader(mesh);
  loader.load(in);

  size_t d = 0;
  EXPECT_EQ(mesh.getCount(d), 4);
  EXPECT_EQ(mesh.getCount(d), mesh.getVertexCount());

  // MFEM format does not support vertex attributes
  EXPECT_EQ(mesh.getAttribute(d, 0), RODIN_DEFAULT_POLYTOPE_ATTRIBUTE);
  EXPECT_EQ(mesh.getAttribute(d, 1), RODIN_DEFAULT_POLYTOPE_ATTRIBUTE);
  EXPECT_EQ(mesh.getAttribute(d, 2), RODIN_DEFAULT_POLYTOPE_ATTRIBUTE);
  EXPECT_EQ(mesh.getAttribute(d, 3), RODIN_DEFAULT_POLYTOPE_ATTRIBUTE);

  d = 1;
  EXPECT_EQ(mesh.getCount(d), 5);
  EXPECT_EQ(mesh.getAttribute(d, 0), 1);
  EXPECT_EQ(mesh.getAttribute(d, 1), 2);

  d = 2;
  EXPECT_EQ(mesh.getCount(d), 2);
  EXPECT_EQ(mesh.getCount(d), mesh.getElementCount());
  EXPECT_EQ(mesh.getAttribute(d, 0), 1);
  EXPECT_EQ(mesh.getAttribute(d, 1), 2);
}

TEST(Rodin_IO_MeshLoader, SanityTest_MFEM_2D_StarSquare)
{
  static constexpr const char* filename = "mfem/StarSquare.mfem.mesh";
  static constexpr size_t vcount = 101;
  static constexpr size_t ecount = 80;
  static constexpr size_t sdim = 2;

  boost::filesystem::path meshfile;
  meshfile = boost::filesystem::path(RODIN_RESOURCES_DIR);
  meshfile.append(filename);

  boost::filesystem::ifstream in(meshfile);

  Mesh mesh;
  MeshLoader<FileFormat::MFEM, Rodin::Context::Serial> loader(mesh);
  loader.load(in);

  EXPECT_EQ(mesh.getSpaceDimension(), sdim);

  size_t d = 0;
  EXPECT_EQ(mesh.getCount(d), vcount);
  EXPECT_EQ(mesh.getCount(d), mesh.getVertexCount());

  for (size_t i = 0; i < vcount; i++)
    EXPECT_EQ(mesh.getAttribute(d, i), RODIN_DEFAULT_POLYTOPE_ATTRIBUTE);

  d = 2;
  for (size_t i = 0; i < ecount; i++)
    EXPECT_EQ(mesh.getAttribute(d, i), 1);

  for (size_t i = 0; i < ecount; i++)
    EXPECT_EQ(mesh.getGeometry(d, i), Polytope::Geometry::Quadrilateral);
}
