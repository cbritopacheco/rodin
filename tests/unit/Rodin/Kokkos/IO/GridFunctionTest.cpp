/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <sstream>
#include <Rodin/IO/MFEM.h>
#include <Rodin/IO/MEDIT.h>
#include <Rodin/IO/HDF5.h>
#include "../Variational/Common.h"

namespace Rodin::Tests::Unit
{
  using namespace Rodin::Variational;

  /// @brief Compares format output with an independent host field and reads it back.
  template <IO::FileFormat Format, class HostFES, class BackendFES>
  void checkFieldIO(const HostFES& hostSpace, const BackendFES& backendSpace)
  {
    GridFunction host(hostSpace);
    Math::Vector<Real> values(hostSpace.getSize());
    for (Index i = 0; i < static_cast<Index>(values.size()); ++i)
      values(i) = Real(i + 1) / Real(values.size() + 1);
    host.setData(values);
    GridFunction backend(backendSpace);
    backend.setData(values);
    std::ostringstream expected, actual;
    IO::GridFunctionPrinter<Format, HostFES, Math::Vector<Real>>(host).print(expected);
    using Data = std::decay_t<decltype(backend.getData())>;
    IO::GridFunctionPrinter<Format, BackendFES, Data> printer(backend);
    printer.print(actual);
    EXPECT_EQ(actual.str(), expected.str());

    GridFunction loaded(backendSpace);
    IO::GridFunctionLoader<Format, BackendFES, Data> loader(loaded);
    std::istringstream input(expected.str());
    loader.load(input);
    expectSameVector(values, loaded.getData());

    const auto path = boost::filesystem::temp_directory_path() /
      boost::filesystem::unique_path("rodin-field-%%%%-%%%%");
    printer.print(path);
    loaded.getData().setZero();
    loader.load(path);
    boost::filesystem::remove(path);
    expectSameVector(values, loaded.getData());
  }

  /// @brief P1 text formats select their local parser through both overloads.
  TEST(Kokkos_IO_GridFunction, P1TextFormats)
  {
    auto source = HostMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    source.getConnectivity().compute(2, 1);
    source.getConnectivity().compute(1, 0);
    const DeviceMesh mesh(TestContext{}, source);
    P1 hostSpace(source);
    P1 backendSpace(mesh);
    checkFieldIO<IO::FileFormat::MFEM>(hostSpace, backendSpace);
    checkFieldIO<IO::FileFormat::MEDIT>(hostSpace, backendSpace);
  }

  /// @brief Higher-order MFEM output retains the host node conversion and ordering.
  TEST(Kokkos_IO_GridFunction, CubicMFEM)
  {
    auto source = HostMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    source.getConnectivity().compute(2, 1);
    source.getConnectivity().compute(1, 0);
    const DeviceMesh mesh(TestContext{}, source);
    constexpr std::integral_constant<size_t, 3> degree{};
    H1 hostSpace(degree, source);
    H1 backendSpace(degree, mesh);
    checkFieldIO<IO::FileFormat::MFEM>(hostSpace, backendSpace);
  }

  /// @brief HDF5 dispatch keeps its specialized path implementation.
  TEST(Kokkos_IO_GridFunction, HDF5Path)
  {
    auto source = HostMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    source.getConnectivity().compute(2, 1);
    source.getConnectivity().compute(1, 0);
    const DeviceMesh mesh(TestContext{}, source);
    P1 space(mesh);
    GridFunction field(space);
    Math::Vector<Real> values = Math::Vector<Real>::LinSpaced(space.getSize(), 1, 2);
    field.setData(values);
    using FES = decltype(space);
    using Data = std::decay_t<decltype(field.getData())>;
    const auto path = boost::filesystem::temp_directory_path() /
      boost::filesystem::unique_path("rodin-field-%%%%-%%%%.h5");
    IO::GridFunctionPrinter<IO::FileFormat::HDF5, FES, Data>(field).print(path);
    field.getData().setZero();
    IO::GridFunctionLoader<IO::FileFormat::HDF5, FES, Data>(field).load(path);
    boost::filesystem::remove(path);
    expectSameVector(values, field.getData());
  }
}
