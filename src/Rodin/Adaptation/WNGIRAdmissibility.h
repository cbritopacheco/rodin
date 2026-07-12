/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file WNGIRAdmissibility.h
 * @brief Admissibility constraints and safeguards for WNGIR displacements.
 */
#ifndef RODIN_ADAPTATION_WNGIRADMISSIBILITY_H
#define RODIN_ADAPTATION_WNGIRADMISSIBILITY_H

#include <algorithm>
#include <cstddef>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "Rodin/Math/SpatialMatrix.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/QF/PolytopeQuadratureFormula.h"
#include "Rodin/Types.h"
#include "Rodin/Variational/IntegrationPoint.h"
#include "Rodin/Variational/Jacobian.h"

namespace Rodin::Adaptation
{
  struct WNGIRAdmissibilityReport
  {
    Real minJ = std::numeric_limits<Real>::infinity();
    std::size_t inadmissibleCount = 0;
    Real maxQRel = Real(0);
  };

  inline std::size_t wngirAdmissibilityQuadratureOrder(
      std::size_t feOrder)
  {
    return std::max<std::size_t>(2, 2 * feOrder);
  }

  template <class Displacement>
  WNGIRAdmissibilityReport evaluateWNGIRAdmissibilitySampled(
      Displacement& u,
      const Math::Vector<Real>& uData,
      Real jMin,
      std::size_t quadratureOrder = 0)
  {
    using Variational::IntegrationPoint;
    using Variational::Jacobian;

    auto& uMutable = u;
    uMutable.getData() = uData;

    WNGIRAdmissibilityReport rep;
    const auto& fes = uMutable.getFiniteElementSpace();
    const auto& mesh = fes.getMesh();
    const std::size_t dim = mesh.getDimension();
    const std::size_t vdim = fes.getVectorDimension();
    if (dim != vdim)
      throw std::runtime_error(
          "evaluateWNGIRAdmissibilitySampled: displacement dimension "
          "must equal mesh dimension.");

    auto gradU = Jacobian(uMutable);
    for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
    {
      const auto& cell = *cellIt;
      const auto& fe =
        fes.getFiniteElement(cell.getDimension(), cell.getIndex());
      const auto& qf =
        QF::PolytopeQuadratureFormula::get(
            quadratureOrder > 0
              ? quadratureOrder
              : wngirAdmissibilityQuadratureOrder(fe.getOrder()),
            cell.getGeometry());
      const auto& quadrature = cell.getQuadrature(qf);
      for (std::size_t q = 0; q < quadrature.getSize(); ++q)
      {
        const auto& pt = quadrature.getPoint(q);
        const IntegrationPoint ip(pt, &qf, q);
        const Math::SpatialMatrix<Real> F =
          Math::SpatialMatrix<Real>::Identity(dim, dim)
          + gradU.getValue(ip);
        const Real j = F.determinant();

        rep.minJ = std::min(rep.minJ, j);
        if (j <= jMin)
          ++rep.inadmissibleCount;

        if (j > Real(0))
        {
          const Real qRel =
            F.squaredNorm()
            / (static_cast<Real>(dim)
               * std::pow(j, Real(2) / static_cast<Real>(dim)));
          rep.maxQRel = std::max(rep.maxQRel, qRel);
        }
      }
    }
    return rep;
  }
}

#endif
