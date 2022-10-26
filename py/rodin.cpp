/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <pybind11/pybind11.h>

#include <Rodin/IO.h>
#include <Rodin/Geometry.h>
#include <Rodin/Alert.h>

#include "type_cast.h"

namespace py = pybind11;

namespace Rodin
{
  namespace Geometry
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

        PyMesh& initialize(int dim, int sdim)
        {
          serial().initialize(dim, sdim);
          return *this;
        }

        PyMesh& vertex(const std::vector<double>& x)
        {
          serial().vertex(x);
          return *this;
        }

        PyMesh& element(
              Type geom,
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
}

PYBIND11_MODULE(rodin, m)
{
  // Rodin::Alert
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
  py::module geometry = m.def_submodule("geometry");

  // Rodin::Geometry::Type
  py::enum_<Rodin::Geometry::Type>(geometry, "Type")
    .value("Triangle",      Rodin::Geometry::Type::Triangle)
    .value("Cube",          Rodin::Geometry::Type::Cube)
    .value("Invalid",       Rodin::Geometry::Type::Invalid)
    .value("Point",         Rodin::Geometry::Type::Point)
    .value("Prism",         Rodin::Geometry::Type::Prism)
    .value("Pyramid",       Rodin::Geometry::Type::Pyramid)
    .value("Segment",       Rodin::Geometry::Type::Segment)
    .value("Square",        Rodin::Geometry::Type::Square)
    .value("Tetrahedron",   Rodin::Geometry::Type::Tetrahedron)
    ;

  // Rodin::Geometry::Mesh
  py::class_<Rodin::Geometry::PyMesh>(geometry, "Mesh")
    .def(py::init<>())
    .def("save", &Rodin::Geometry::PyMesh::save,
        py::arg("filename"),
        py::arg("fmt") = Rodin::IO::FileFormat::MFEM,
        py::arg("precision") = 16)
    .def("is_parallel", &Rodin::Geometry::PyMesh::isParallel)
    .def("initialize", py::overload_cast<int, int>(&Rodin::Geometry::PyMesh::initialize),
        py::arg("dim"), py::arg("sdim"))
    .def("vertex", &Rodin::Geometry::PyMesh::vertex)
    .def("element", &Rodin::Geometry::PyMesh::element,
        py::arg("geometry"), py::arg("vertices"), py::arg("attribute") = py::none())
    ;
}
