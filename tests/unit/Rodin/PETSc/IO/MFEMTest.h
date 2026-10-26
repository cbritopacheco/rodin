/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_TEST_PETSC_IO_MFEMTEST_H
#define RODIN_TEST_PETSC_IO_MFEMTEST_H
/** @file @brief Shared local and distributed MFEM coefficient and analytic oracles. */
#include <gtest/gtest.h>
#include <petsc.h>
#include <memory>
#include <sstream>
#include <type_traits>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/PETSc.h>
#include <Rodin/PETSc/IO.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::IO;

namespace MFEMTest
{
  template <size_t Family, size_t Rank, class MeshType, class Scalar = PetscScalar>
  void checkExport(MeshType& mesh)
  {
    using Range = std::conditional_t<Rank == 0, Scalar,
      std::conditional_t<Rank == 1, Math::SpatialVector<Scalar>,
        Math::SpatialMatrix<Scalar>>>;
    const size_t dimension = mesh.getDimension();
    auto space = [&] {
      if constexpr (Family == 0)
      {
        if constexpr (Rank == 0)
          return P0<Range, MeshType>(mesh);
        else if constexpr (Rank == 1)
          return P0<Range, MeshType>(mesh, 2);
        else
          return P0<Range, MeshType>(mesh, 2, 2);
      }
      else if constexpr (Family == 1)
      {
        if constexpr (Rank == 0)
          return P1<Range, MeshType>(mesh);
        else if constexpr (Rank == 1)
          return P1<Range, MeshType>(mesh, 2);
        else
          return P1<Range, MeshType>(mesh, 2, 2);
      }
      else
      {
        if constexpr (Rank == 0)
          return H1<3, Range, MeshType>(std::integral_constant<size_t, 3>{}, mesh);
        else if constexpr (Rank == 1)
          return H1<3, Range, MeshType>(std::integral_constant<size_t, 3>{}, mesh, 2);
        else
          return H1<3, Range, MeshType>(std::integral_constant<size_t, 3>{}, mesh, 2, 2);
      }
    }();
    const Scalar factor = [] {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return Scalar(1, -2);
      else
        return Scalar(1);
    }();
    PETSc::Variational::GridFunction field(space);
    const auto& localSpace = [&]() -> const auto& {
      if constexpr (requires { space.getShard(); })
        return space.getShard();
      else
        return space;
    }();
    using LocalFES = std::remove_cvref_t<decltype(localSpace)>;
    GridFunction expected(localSpace);
    using FES = std::decay_t<decltype(space)>;
    using Base = GridFunctionPrinterBase<FileFormat::MFEM, FES, Vec>;
    std::unique_ptr<Base> printer(
      new GridFunctionPrinter<FileFormat::MFEM, FES, Vec>(field));
    EXPECT_EQ(&printer->getObject(), &field);
    // Cubic node conversions and 17-digit formatting contribute rounding error.
    constexpr Real PersistenceTolerance = 1e-10;
    for (size_t repeat = 0; repeat < 2; ++repeat)
    {
      const auto value = [=](const Point& point) -> Range {
        auto component = [&](size_t row) {
          Real polynomial = Real(1 + row + repeat);
          for (size_t axis = 0; axis < dimension; ++axis)
            polynomial += Real((row + 1) * (axis + 2)) * point(axis);
          if constexpr (Family != 1)
            polynomial += Real(row + 2) * point(0) * point(0) / 8;
          return Scalar(factor * polynomial);
        };
        if constexpr (Rank == 1)
        {
          Range result(2);
          for (size_t row = 0; row < 2; ++row)
            result[row] = component(row);
          return result;
        }
        else if constexpr (Rank == 2)
        {
          Range result(2, 2);
          for (size_t row = 0; row < 2; ++row)
          {
            for (size_t column = 0; column < 2; ++column)
              result(row, column) = component(2 * row + column);
          }
          return result;
        }
        else
          return component(0);
      };
      field = value;
      expected = value;
      std::stringstream stream;
      printer->print(stream);
      GridFunction restored(localSpace);
      restored.getData().setConstant(Scalar(-123));
      GridFunctionLoader<FileFormat::MFEM, LocalFES, Math::Vector<Scalar>>(restored).load(
        stream);
      if (expected.getData().size() != 0)
      {
        EXPECT_LT(
          (restored.getData() - expected.getData()).template lpNorm<Eigen::Infinity>(),
          PersistenceTolerance);
      }
      stream >> std::ws;
      EXPECT_TRUE(stream.eof());
      for (auto cell = localSpace.getMesh().getCell(); cell; ++cell)
      {
        const Point point(*cell, Polytope::Traits(cell->getGeometry()).getCentroid());
        if constexpr (Rank != 0)
        {
          EXPECT_LT(
            (restored.getValue(point) - value(point)).norm(), PersistenceTolerance);
          EXPECT_LT((field.getValue(point) - value(point)).norm(), PersistenceTolerance);
        }
        else
        {
          EXPECT_LT(
            std::abs(restored.getValue(point) - value(point)), PersistenceTolerance);
          EXPECT_LT(std::abs(field.getValue(point) - value(point)), PersistenceTolerance);
        }
      }
    }
  }

}
#endif
