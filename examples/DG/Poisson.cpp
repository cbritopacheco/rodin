/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
//! [dg-source]
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

#include <Eigen/SparseCholesky>
#include <Rodin/Geometry.h>
#include <Rodin/QF.h>
#include <Rodin/Variational.h>
#include <Rodin/IO/XDMF.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main()
{
  constexpr Real pi = 3.14159265358979323846;
  constexpr Real penalty = 20;
  constexpr Real algebraTolerance = 1e-10;
  std::ofstream history("DGConvergence.csv");
  history.precision(17);
  history << "intervals,dofs,l2,h1,jump,residual,symmetry,min_pivot,patch\n";
  bool valid = true;
  Real previousL2 = 0, previousH1 = 0, previousJump = 0;

  //! [dg-cases]
  // The affine patch tests nonzero boundary data; the smooth case tests rates.
  for (const bool patch : {true, false})
  {
    for (const size_t intervals : {size_t{4}, size_t{8}, size_t{16}, size_t{32}})
    {
      if (patch && intervals != 4)
        continue;
      Mesh mesh =
        Mesh().UniformGrid(Polytope::Type::Triangle, {intervals + 1, intervals + 1});
      mesh.scale(1.0 / intervals);
      mesh.getConnectivity().compute(1, 2);
      const size_t count = mesh.getCellCount();
      const size_t dofs = 3 * count;
      const auto exact = [&](const Point& p) {
        return 1 + p.x() + 2 * p.y() +
          (patch ? 0 : std::sin(pi * p.x()) * std::sin(pi * p.y()));
      };
      const auto exactGradient = [&](const Point& p) {
        return Math::SpatialVector<Real>{
          1 + (patch ? 0 : pi * std::cos(pi * p.x()) * std::sin(pi * p.y())),
          2 + (patch ? 0 : pi * std::sin(pi * p.x()) * std::cos(pi * p.y()))};
      };
      //! [dg-cases]

      //! [dg-local-space]
      // Affine triangle basis: (1-r-s, r, s). Coefficients belong to cells.
      const std::array<Math::SpatialVector<Real>, 3> referenceGradients = {
        Math::SpatialVector<Real>{-1, -1}, {1, 0}, {0, 1}};
      std::vector<std::array<Math::SpatialVector<Real>, 3>> gradients(count);
      for (auto cell = mesh.getCell(); !cell.end(); ++cell)
      {
        const Point center(*cell, Math::SpatialPoint{1.0 / 3, 1.0 / 3});
        for (size_t i = 0; i < 3; ++i)
        {
          gradients[cell->getIndex()][i] =
            center.getJacobianInverse().transpose() * referenceGradients[i];
        }
      }
      Math::SparseMatrix<Real> matrix(dofs, dofs);
      Math::Vector<Real> rhs = Math::Vector<Real>::Zero(dofs);
      std::vector<Math::SparseTriplet<Real>> entries;
      //! [dg-local-space]

      //! [dg-cells]
      const auto& cellFormula =
        QF::QuadratureFormula::get(size_t{10}, Polytope::Type::Triangle);
      for (auto cell = mesh.getCell(); !cell.end(); ++cell)
      {
        const Index k = cell->getIndex();
        for (size_t i = 0; i < 3; ++i)
        {
          for (size_t j = 0; j < 3; ++j)
          {
            entries.emplace_back(3 * k + i, 3 * k + j,
              cell->getMeasure() * gradients[k][i].dot(gradients[k][j]));
          }
        }
        PolytopeQuadrature points(*cell, cellFormula);
        for (size_t q = 0; q < points.getSize(); ++q)
        {
          const auto& p = points.getPoint(q);
          const auto& r = cellFormula.getPoint(q);
          const std::array<Real, 3> basis = {1 - r.x() - r.y(), r.x(), r.y()};
          const Real f =
            patch ? 0 : 2 * pi * pi * std::sin(pi * p.x()) * std::sin(pi * p.y());
          const Real weight = cellFormula.getWeight(q) * p.getJacobianDeterminant();
          for (size_t i = 0; i < 3; ++i)
            rhs[3 * k + i] += weight * f * basis[i];
        }
      }
      //! [dg-cells]

      //! [dg-facets]
      const auto& faceFormula =
        QF::QuadratureFormula::get(size_t{10}, Polytope::Type::Segment);
      for (auto face = mesh.getFace(); !face.end(); ++face)
      {
        const auto& incident =
          mesh.getConnectivity().getIncidence({1, 2}, face->getIndex());
        std::vector<Index> cells(incident.begin(), incident.end());
        const bool boundary = cells.size() == 1;
        const auto firstCell = mesh.getCell(cells[0]);
        const Point center(*firstCell, Math::SpatialPoint{1.0 / 3, 1.0 / 3});
        const Point midpoint(*face, Math::SpatialPoint{0.5});
        const auto vertices = face->getVertex();
        auto next = face->getVertex();
        ++next;
        const Math::SpatialVector<Real> tangent =
          next->getCoordinates() - vertices->getCoordinates();
        Math::SpatialVector<Real> normal{tangent.y(), -tangent.x()};
        const Real length = normal.norm();
        normal /= length;
        if (normal.dot(
              midpoint.getPhysicalCoordinates() - center.getPhysicalCoordinates()) < 0)
          normal *= -1;
        const Real sigma = penalty / length;
        PolytopeQuadrature points(*face, faceFormula);
        for (size_t q = 0; q < points.getSize(); ++q)
        {
          const auto& p = points.getPoint(q);
          const Real weight = faceFormula.getWeight(q) * p.getJacobianDeterminant();
          std::array<Real, 6> jump{};
          std::array<Real, 6> normalDerivative{};
          for (size_t side = 0; side < cells.size(); ++side)
          {
            const auto cell = mesh.getCell(cells[side]);
            Math::SpatialPoint r;
            cell->getTransformation().inverse(r, p.getPhysicalCoordinates());
            const std::array<Real, 3> basis = {1 - r.x() - r.y(), r.x(), r.y()};
            for (size_t i = 0; i < 3; ++i)
            {
              const size_t local = 3 * side + i;
              jump[local] = (side == 0 ? 1 : -1) * basis[i];
              normalDerivative[local] =
                (boundary ? 1 : 0.5) * gradients[cells[side]][i].dot(normal);
            }
          }
          const size_t localSize = 3 * cells.size();
          for (size_t i = 0; i < localSize; ++i)
          {
            const Index row = 3 * cells[i / 3] + i % 3;
            for (size_t j = 0; j < localSize; ++j)
            {
              const Index column = 3 * cells[j / 3] + j % 3;
              const Real value = -normalDerivative[j] * jump[i] -
                normalDerivative[i] * jump[j] + sigma * jump[i] * jump[j];
              entries.emplace_back(row, column, weight * value);
            }
            if (boundary)
              rhs[row] += weight * exact(p) * (-normalDerivative[i] + sigma * jump[i]);
          }
        }
      }
      //! [dg-facets]

      //! [dg-solve]
      matrix.setFromTriplets(entries.begin(), entries.end());
      const Math::SparseMatrix<Real> transpose = matrix.transpose();
      const Real symmetry = (matrix - transpose).norm() / matrix.norm();
      Eigen::SimplicialLDLT<Math::SparseMatrix<Real>> solver;
      solver.compute(matrix);
      if (solver.info() != Eigen::Success || solver.vectorD().minCoeff() <= 0)
        return EXIT_FAILURE;
      const Math::Vector<Real> solution = solver.solve(rhs);
      if (solver.info() != Eigen::Success)
        return EXIT_FAILURE;
      const Real residual = (matrix * solution - rhs).norm() / rhs.norm();
      valid = valid && symmetry < algebraTolerance && residual < algebraTolerance;
      //! [dg-solve]

      //! [dg-errors]
      Real l2 = 0, h1 = 0, jumpNorm = 0;
      for (auto cell = mesh.getCell(); !cell.end(); ++cell)
      {
        const Index k = cell->getIndex();
        Math::SpatialVector<Real> gradient{0, 0};
        for (size_t i = 0; i < 3; ++i)
          gradient += solution[3 * k + i] * gradients[k][i];
        PolytopeQuadrature points(*cell, cellFormula);
        for (size_t q = 0; q < points.getSize(); ++q)
        {
          const auto& p = points.getPoint(q);
          const auto& r = cellFormula.getPoint(q);
          const Real value = solution[3 * k] * (1 - r.x() - r.y()) +
            solution[3 * k + 1] * r.x() + solution[3 * k + 2] * r.y();
          const Real weight = cellFormula.getWeight(q) * p.getJacobianDeterminant();
          l2 += weight * std::pow(value - exact(p), 2);
          h1 += weight * (gradient - exactGradient(p)).squaredNorm();
        }
      }
      for (auto face = mesh.getInterface(); !face.end(); ++face)
      {
        const auto& incident =
          mesh.getConnectivity().getIncidence({1, 2}, face->getIndex());
        PolytopeQuadrature points(*face, faceFormula);
        for (size_t q = 0; q < points.getSize(); ++q)
        {
          const auto& p = points.getPoint(q);
          Real difference = 0;
          Real sign = 1;
          for (const Index k : incident)
          {
            const auto cell = mesh.getCell(k);
            Math::SpatialPoint r;
            cell->getTransformation().inverse(r, p.getPhysicalCoordinates());
            difference += sign *
              (solution[3 * k] * (1 - r.x() - r.y()) + solution[3 * k + 1] * r.x() +
                solution[3 * k + 2] * r.y());
            sign *= -1;
          }
          jumpNorm += faceFormula.getWeight(q) * p.getJacobianDeterminant() * difference *
            difference;
        }
      }
      l2 = std::sqrt(l2);
      h1 = std::sqrt(h1);
      jumpNorm = std::sqrt(jumpNorm);
      history << intervals << ',' << dofs << ',' << l2 << ',' << h1 << ',' << jumpNorm
              << ',' << residual << ',' << symmetry << ',' << solver.vectorD().minCoeff()
              << ',' << patch << '\n';
      std::cout << (patch ? "Patch" : "Smooth") << ": n=" << intervals
                << ", DOFs=" << dofs << ", L2=" << l2 << ", H1=" << h1
                << ", jump=" << jumpNorm << ", residual=" << residual << std::endl;
      if (patch)
        valid = valid && l2 < algebraTolerance && h1 < algebraTolerance;
      else if (intervals > 4)
      {
        // Asymptotic P1 SIPG rates: L2 order two and broken H1 order one.
        valid = valid && previousL2 / l2 > 3.5 && previousH1 / h1 > 1.8 &&
          jumpNorm < previousJump;
      }
      if (!patch)
      {
        previousL2 = l2;
        previousH1 = h1;
        previousJump = jumpNorm;
      }
      //! [dg-errors]

      //! [dg-output]
      if (!patch && (intervals == 4 || intervals == 32))
      {
        // Separate geometric vertices preserve the two traces in visualization.
        auto builder = Mesh<Context::Local>::Builder();
        builder.initialize(2).nodes(dofs);
        for (auto cell = mesh.getCell(); !cell.end(); ++cell)
        {
          for (auto vertex = cell->getVertex(); !vertex.end(); ++vertex)
            builder.vertex(vertex->getCoordinates());
        }
        for (size_t k = 0; k < count; ++k)
          builder.polytope(Polytope::Type::Triangle, {3 * k, 3 * k + 1, 3 * k + 2});
        Mesh output = builder.finalize();
        P1 outputSpace(output);
        GridFunction field(outputSpace);
        for (size_t k = 0; k < count; ++k)
        {
          const auto& indices = outputSpace.getDOFs(2, k);
          for (size_t i = 0; i < 3; ++i)
            field[indices[i]] = solution[3 * k + i];
        }
        GridFunction reference(outputSpace);
        reference = RealFunction(exact);
        GridFunction error(outputSpace);
        error = field - reference;
        IO::XDMF xdmf(intervals == 4 ? "DGCoarse" : "DGPoisson");
        xdmf.grid().setMesh(output).add("u", field).add("error", error);
        xdmf.write().flush();
      }
      //! [dg-output]
    }
  }
  history.flush();
  return valid && history ? EXIT_SUCCESS : EXIT_FAILURE;
}
//! [dg-source]
