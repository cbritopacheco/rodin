/*
 * Generates MFEM-authored H1 GridFunction fixtures used by
 * H1GridFunctionIOTest.cpp.
 *
 * This source is intentionally not part of the Rodin build. Rebuild fixtures
 * manually against a local MFEM checkout, for example:
 *
 *   c++ -std=c++17 -I/path/to/mfem/build -I/path/to/mfem \
 *     tests/unit/Rodin/IO/MFEMH1FixtureGenerator.cpp \
 *     /path/to/mfem/build/libmfem.a -o /tmp/mfem_h1_fixture_generator
 *   /tmp/mfem_h1_fixture_generator resources/tests/mfem/h1
 */

#include "mfem.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

namespace
{
  double scalar_exact(const mfem::Vector& x)
  {
    const double X = x.Size() > 0 ? x(0) : 0.0;
    const double Y = x.Size() > 1 ? x(1) : 0.0;
    const double Z = x.Size() > 2 ? x(2) : 0.0;

    return 1.0 + 2.0 * X - 3.0 * Y + 0.5 * Z
      + X * X + 0.25 * Y * Y - 0.75 * Z * Z
      + X * Y - 0.4 * X * Z + 0.2 * Y * Z;
  }

  void vector_exact(const mfem::Vector& x, mfem::Vector& out)
  {
    const double X = x.Size() > 0 ? x(0) : 0.0;
    const double Y = x.Size() > 1 ? x(1) : 0.0;
    const double Z = x.Size() > 2 ? x(2) : 0.0;

    out(0) = scalar_exact(x);
    out(1) = -2.0 + X + 4.0 * Y - Z + X * Z + 0.5 * Y * Y;
    out(2) = 3.0 - X * Y + Z * Z + 0.25 * X * X;
  }

  void save_mesh(const std::filesystem::path& path, mfem::Mesh& mesh)
  {
    std::ofstream os(path);
    os << std::setprecision(17);
    mesh.Print(os);
  }

  mfem::Mesh make_mixed_2d_mesh()
  {
    mfem::Mesh mesh(2, 6, 3, 6, 2);

    mesh.AddVertex(0.0, 0.0);
    mesh.AddVertex(0.5, 0.0);
    mesh.AddVertex(1.0, 0.0);
    mesh.AddVertex(0.0, 1.0);
    mesh.AddVertex(0.5, 1.0);
    mesh.AddVertex(1.0, 1.0);

    mesh.AddQuad(0, 1, 4, 3, 1);
    mesh.AddTriangle(1, 2, 5, 2);
    mesh.AddTriangle(1, 5, 4, 2);

    mesh.AddBdrSegment(0, 1, 1);
    mesh.AddBdrSegment(1, 2, 1);
    mesh.AddBdrSegment(2, 5, 2);
    mesh.AddBdrSegment(5, 4, 3);
    mesh.AddBdrSegment(4, 3, 3);
    mesh.AddBdrSegment(3, 0, 4);

    mesh.FinalizeMesh(0, true);
    return mesh;
  }

  void save_scalar_case(
      const std::filesystem::path& out,
      const std::string& name,
      mfem::Mesh mesh,
      int order)
  {
    save_mesh(out / (name + ".mesh"), mesh);

    mfem::H1_FECollection fec(order, mesh.Dimension());
    mfem::FiniteElementSpace fes(
      &mesh, &fec, 1, mfem::Ordering::byNODES);

    mfem::FunctionCoefficient coeff(scalar_exact);
    mfem::GridFunction gf(&fes);
    gf.ProjectCoefficient(coeff);

    std::ofstream os(out / (name + ".gf"));
    os << std::setprecision(17);
    gf.Save(os);
  }

  template <class MeshFactory>
  void save_scalar_order_family(
      const std::filesystem::path& out,
      const std::string& geometry,
      MeshFactory make_mesh)
  {
    for (int order = 1; order <= 6; ++order)
    {
      save_scalar_case(
        out,
        "h1_" + geometry + "_p" + std::to_string(order) + "_scalar_mfem_by_nodes",
        make_mesh(),
        order);
    }
  }

  void save_vector_case(
      const std::filesystem::path& out,
      const std::string& name,
      mfem::Mesh mesh,
      int order,
      int ordering)
  {
    save_mesh(out / (name + ".mesh"), mesh);

    mfem::H1_FECollection fec(order, mesh.Dimension());
    mfem::FiniteElementSpace fes(&mesh, &fec, 3, ordering);

    mfem::VectorFunctionCoefficient coeff(3, vector_exact);
    mfem::GridFunction gf(&fes);
    gf.ProjectCoefficient(coeff);

    std::ofstream os(out / (name + ".gf"));
    os << std::setprecision(17);
    gf.Save(os);
  }

  template <class MeshFactory>
  void save_vector_order_family(
      const std::filesystem::path& out,
      const std::string& geometry,
      const std::string& orderingName,
      int ordering,
      MeshFactory make_mesh)
  {
    for (int order = 1; order <= 6; ++order)
    {
      save_vector_case(
        out,
        "h1_" + geometry + "_p" + std::to_string(order) + "_vector_mfem_" + orderingName,
        make_mesh(),
        order,
        ordering);
    }
  }
}

int main(int argc, char** argv)
{
  if (argc != 2)
  {
    std::cerr << "usage: " << argv[0] << " <output-dir>\n";
    return 2;
  }

  const std::filesystem::path out(argv[1]);
  std::filesystem::create_directories(out);

  save_scalar_order_family(out, "triangle", []()
  {
    return mfem::Mesh::MakeCartesian2D(
      2, 2, mfem::Element::TRIANGLE, true, 1.0, 1.0, false);
  });

  save_scalar_order_family(out, "mixed2d", []()
  {
    return make_mixed_2d_mesh();
  });

  save_scalar_order_family(out, "tetrahedron", []()
  {
    return mfem::Mesh::MakeCartesian3D(
      1, 1, 1, mfem::Element::TETRAHEDRON, 1.0, 1.0, 1.0, false);
  });

  save_scalar_order_family(out, "hexahedron", []()
  {
    return mfem::Mesh::MakeCartesian3D(
      1, 1, 1, mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0, false);
  });

  save_scalar_order_family(out, "wedge", []()
  {
    return mfem::Mesh::MakeCartesian3D(
      1, 1, 1, mfem::Element::WEDGE, 1.0, 1.0, 1.0, false);
  });

  save_vector_order_family(out, "triangle", "by_nodes", mfem::Ordering::byNODES, []()
  {
    return mfem::Mesh::MakeCartesian2D(
      2, 2, mfem::Element::TRIANGLE, true, 1.0, 1.0, false);
  });

  save_vector_order_family(out, "mixed2d", "by_nodes", mfem::Ordering::byNODES, []()
  {
    return make_mixed_2d_mesh();
  });

  save_vector_order_family(out, "tetrahedron", "by_nodes", mfem::Ordering::byNODES, []()
  {
    return mfem::Mesh::MakeCartesian3D(
      1, 1, 1, mfem::Element::TETRAHEDRON, 1.0, 1.0, 1.0, false);
  });

  save_vector_order_family(out, "hexahedron", "by_nodes", mfem::Ordering::byNODES, []()
  {
    return mfem::Mesh::MakeCartesian3D(
      1, 1, 1, mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0, false);
  });

  save_vector_order_family(out, "wedge", "by_nodes", mfem::Ordering::byNODES, []()
  {
    return mfem::Mesh::MakeCartesian3D(
      1, 1, 1, mfem::Element::WEDGE, 1.0, 1.0, 1.0, false);
  });

  save_vector_order_family(out, "triangle", "by_vdim", mfem::Ordering::byVDIM, []()
  {
    return mfem::Mesh::MakeCartesian2D(
      2, 2, mfem::Element::TRIANGLE, true, 1.0, 1.0, false);
  });

  save_vector_order_family(out, "mixed2d", "by_vdim", mfem::Ordering::byVDIM, []()
  {
    return make_mixed_2d_mesh();
  });

  save_vector_order_family(out, "tetrahedron", "by_vdim", mfem::Ordering::byVDIM, []()
  {
    return mfem::Mesh::MakeCartesian3D(
      1, 1, 1, mfem::Element::TETRAHEDRON, 1.0, 1.0, 1.0, false);
  });

  save_vector_order_family(out, "hexahedron", "by_vdim", mfem::Ordering::byVDIM, []()
  {
    return mfem::Mesh::MakeCartesian3D(
      1, 1, 1, mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0, false);
  });

  save_vector_order_family(out, "wedge", "by_vdim", mfem::Ordering::byVDIM, []()
  {
    return mfem::Mesh::MakeCartesian3D(
      1, 1, 1, mfem::Element::WEDGE, 1.0, 1.0, 1.0, false);
  });

  return 0;
}
