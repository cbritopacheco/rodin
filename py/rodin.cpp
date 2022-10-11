/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <pybind11/pybind11.h>

#include <Rodin/IO.h>
#include <Rodin/Mesh.h>
#include <Rodin/Alert.h>

#include "type_cast.h"

namespace py = pybind11;

namespace Rodin
{
  class PyMesh
  {
    using PImpl = std::unique_ptr<MeshBase>;

    public:
      PyMesh()
        : m_impl(new Mesh<Context::Serial>)
      {}

      PyMesh(Mesh<Context::Serial>&& impl)
        : m_impl(new Mesh<Context::Serial>(std::move(impl)))
      {}

      PyMesh(Mesh<Context::Parallel>&& impl)
        : m_impl(new Mesh<Context::Parallel>(std::move(impl)))
      {}

      void save(
          const std::string& filename,
          IO::FileFormat fmt,
          int precison) const
      {
        return serial().save(filename, fmt, precison);
      }

      PyMesh& initialize(int dim, int sdim, int nvert)
      {
        serial().initialize(dim, sdim, nvert);
        return *this;
      }

      PyMesh& vertex(const std::vector<double>& x)
      {
        serial().vertex(x);
        return *this;
      }

      PyMesh& element(
            Geometry geom,
            const std::vector<int>& vs,
            std::optional<int> attr = {})
      {
        serial().element(geom, vs, attr);
        return *this;
      }

      PyMesh& finalize()
      {
        serial().finalize();
        return *this;
      }

      bool isParallel() const
      {
        return m_impl->isParallel();
      }

    private:
      Mesh<Context::Serial>& serial()
      {
        if (m_impl->isParallel())
          Alert::Exception() << "Function requires a serial context." << Alert::Raise;
        return static_cast<Mesh<Context::Serial>&>(*m_impl);
      }

      const Mesh<Context::Serial>& serial() const
      {
        if (m_impl->isParallel())
          Alert::Exception() << "Function requires a serial context." << Alert::Raise;
        return static_cast<const Mesh<Context::Serial>&>(*m_impl);
      }

      Mesh<Context::Parallel>& parallel()
      {
        if (!m_impl->isParallel())
          Alert::Exception() << "Function requires a parallel context." << Alert::Raise;
        return static_cast<Mesh<Context::Parallel>&>(*m_impl);
      }

      const Mesh<Context::Parallel>& parallel() const
      {
        if (!m_impl->isParallel())
          Alert::Exception() << "Function requires a parallel context." << Alert::Raise;
        return static_cast<const Mesh<Context::Parallel>&>(*m_impl);
      }

      PImpl m_impl;
  };
}

PYBIND11_MODULE(rodin, m)
{
  py::module alert = m.def_submodule("alert");
  py::register_exception<Rodin::Alert::Exception>(alert, "Exception");

  // Rodin::IO
  py::module io = m.def_submodule("io");

  // Rodin::IO::FileFormat
  py::enum_<Rodin::IO::FileFormat>(io, "FileFormat")
    .value("MFEM", Rodin::IO::FileFormat::MFEM)
    .value("GMSH", Rodin::IO::FileFormat::GMSH)
    .value("MEDIT", Rodin::IO::FileFormat::MEDIT)
    ;

  // Rodin::Geometry
  py::enum_<Rodin::Geometry>(m, "Geometry")
    .value("Triangle",      Rodin::Geometry::Triangle)
    .value("Cube",          Rodin::Geometry::Cube)
    .value("Invalid",       Rodin::Geometry::Invalid)
    .value("Point",         Rodin::Geometry::Point)
    .value("Prism",         Rodin::Geometry::Prism)
    .value("Pyramid",       Rodin::Geometry::Pyramid)
    .value("Segment",       Rodin::Geometry::Segment)
    .value("Square",        Rodin::Geometry::Square)
    .value("Tetrahedron",   Rodin::Geometry::Tetrahedron)
    ;

  // Rodin::Mesh
  py::class_<Rodin::PyMesh>(m, "Mesh")
    .def(py::init<>())
    .def("save", &Rodin::PyMesh::save,
        py::arg("filename"),
        py::arg("fmt") = Rodin::IO::FileFormat::MFEM,
        py::arg("precision") = 16)
    .def("is_parallel", &Rodin::PyMesh::isParallel)
    .def("initialize", &Rodin::PyMesh::initialize)
    .def("vertex", &Rodin::PyMesh::vertex)
    .def("element", &Rodin::PyMesh::element,
        py::arg("geometry"), py::arg("vertices"), py::arg("attribute") = py::none())
    ;
}
