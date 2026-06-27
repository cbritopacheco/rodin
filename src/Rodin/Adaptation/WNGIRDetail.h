/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_WNGIRDETAIL_H
#define RODIN_ADAPTATION_WNGIRDETAIL_H

#include "WNGIRCommon.h"

namespace Rodin::Adaptation
{
  namespace Detail
  {
    inline Real referenceCellMargin(
        Geometry::Polytope::Type geometry,
        const Math::SpatialPoint& rc,
        std::size_t& mostViolatedFace)
    {
      const Geometry::Polytope::Traits traits(geometry);
      const auto& hs = traits.getHalfSpace();
      Real margin = std::numeric_limits<Real>::infinity();
      mostViolatedFace = 0;
      for (std::size_t j = 0; j < static_cast<std::size_t>(hs.vector.size()); ++j)
      {
        const Real phi = hs.vector[j] - rc.dot(hs.matrix.row(j).transpose());
        if (phi < margin)
        {
          margin = phi;
          mostViolatedFace = j;
        }
      }
      return margin;
    }

    inline Geometry::Point makeTranslatedPoint(
        const Geometry::Point& source,
        const Math::SpatialVector<Real>& yPhysical,
        Real tolerance)
    {
      const auto& sourcePolytope = source.getPolytope();
      const auto& mesh = sourcePolytope.getMesh();
      const std::size_t cd = mesh.getDimension();
      const auto& conn = mesh.getConnectivity();
      const Real tol = tolerance > Real(0)
        ? tolerance
        : Real(64) * std::numeric_limits<Real>::epsilon();

      Index cell = Index(-1);
      Math::SpatialPoint rc;

      auto cellMargin =
        [&](Index candidate, Math::SpatialPoint& out) -> Real
        {
          mesh.getPolytopeTransformation(cd, candidate).inverse(
              out, yPhysical);
          std::size_t face = 0;
          return referenceCellMargin(mesh.getGeometry(cd, candidate), out, face);
        };

      if (sourcePolytope.getDimension() == cd)
      {
        cell = sourcePolytope.getIndex();
        cellMargin(cell, rc);
      }
      else if (sourcePolytope.getDimension() + 1 == cd)
      {
        const Index face = sourcePolytope.getIndex();
        const auto& adjacent = conn.getIncidence(cd - 1, cd).at(face);
        if (adjacent.empty())
          throw std::runtime_error(
              "WNGIR translated-point evaluation failed: face has no adjacent cell.");

        Real bestMargin = -std::numeric_limits<Real>::infinity();
        Math::SpatialPoint bestRc;
        for (const Index candidate : adjacent)
        {
          Math::SpatialPoint candidateRc;
          const Real margin = cellMargin(candidate, candidateRc);
          if (margin > bestMargin)
          {
            bestMargin = margin;
            bestRc = candidateRc;
            cell = candidate;
          }
        }
        rc = bestRc;
      }
      else
      {
        throw std::runtime_error(
            "WNGIR translated-point evaluation requires a cell or face source point.");
      }

      for (std::size_t hop = 0; hop < 64; ++hop)
      {
        mesh.getPolytopeTransformation(cd, cell).inverse(rc, yPhysical);

        std::size_t mostViolatedFace = 0;
        const Real margin = referenceCellMargin(
            mesh.getGeometry(cd, cell), rc, mostViolatedFace);
        if (margin >= -tol)
        {
          const auto it = mesh.getPolytope(cd, cell);
          return Geometry::Point(*it, rc, yPhysical);
        }

        const auto& faces = conn.getIncidence(cd, cd - 1).at(cell);
        if (mostViolatedFace >= faces.size())
          break;

        const Index face = faces[mostViolatedFace];
        if (mesh.isBoundary(face))
          break;

        const auto& adjacent = conn.getIncidence(cd - 1, cd).at(face);
        if (adjacent.size() != 2)
          break;

        const Index next = (adjacent[0] == cell) ? adjacent[1] : adjacent[0];
        if (next == cell)
          break;
        cell = next;
      }

      const auto it = mesh.getPolytope(cd, cell);
      mesh.getPolytopeTransformation(cd, cell).inverse(rc, yPhysical);
      return Geometry::Point(*it, rc, yPhysical);
    }

    template <class Function, class Vector>
    decltype(auto) evaluateTranslatedPoint(
        const Function& f,
        const Geometry::Point& source,
        const Vector& displacement,
        Real tolerance)
    {
      Math::SpatialVector<Real> y(source.getPolytope().getMesh().getDimension());
      const auto& x = source.getPhysicalCoordinates();
      for (std::size_t r = 0; r < static_cast<std::size_t>(y.size()); ++r)
        y(static_cast<Eigen::Index>(r)) =
          x(static_cast<Eigen::Index>(r))
          + displacement(static_cast<Eigen::Index>(r));
      const auto p = makeTranslatedPoint(source, y, tolerance);
      if constexpr (requires { f.getValue(p); })
        return f.getValue(p);
      else
        return f(p);
    }

    inline Math::SpatialMatrix<Real> makeZeroMatrix(std::size_t dim)
    {
      Math::SpatialMatrix<Real> out(dim, dim);
      out.setZero();
      return out;
    }

    template <class FE, class ReferencePoint, class JacobianInverse>
    Math::SpatialMatrix<Real> physicalJacobian(
        const FE& fe,
        std::size_t local,
        const ReferencePoint& rc,
        const JacobianInverse& Jinv,
        std::size_t dim)
    {
      const auto jref = fe.getBasis(local).getJacobian()(rc);
      auto jp = makeZeroMatrix(dim);
      for (std::size_t r = 0; r < dim; ++r)
        for (std::size_t c = 0; c < dim; ++c)
          for (std::size_t a = 0; a < dim; ++a)
            jp(static_cast<Eigen::Index>(r), static_cast<Eigen::Index>(c))
              += jref(r, a) * Jinv(a, c);
      return jp;
    }

    template <class Displacement>
    Math::SpatialMatrix<Real> deformationGradient(
        const Displacement& u,
        const Geometry::Polytope& polytope,
        const Variational::IntegrationPoint& ip,
        std::size_t dim)
    {
      const auto& fes = u.getFiniteElementSpace();
      const auto& fe = fes.getFiniteElement(
          polytope.getDimension(), polytope.getIndex());
      const auto& dofs = fes.getDOFs(
          polytope.getDimension(), polytope.getIndex());
      const auto Jinv = ip.getPoint().getJacobianInverse();
      const auto& rc = ip.getPoint().getReferenceCoordinates();

      auto F = Math::SpatialMatrix<Real>::Identity(dim, dim);
      for (std::size_t l = 0; l < fe.getCount(); ++l)
        F += u.getData()(dofs[l])
          * physicalJacobian(fe, l, rc, Jinv, dim);
      return F;
    }
  }
}

#endif
