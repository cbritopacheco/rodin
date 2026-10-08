/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_WNGIR_ADMISSIBILITY_H
#define RODIN_ADAPTATION_WNGIR_ADMISSIBILITY_H

#include <algorithm>
#include <cstddef>
#include <cmath>
#include <limits>
#include "Rodin/Alert.h"

#include "Rodin/QF/PolytopeQuadratureFormula.h"
#include "Rodin/Types.h"
#include "Rodin/Variational/IntegrationPoint.h"
#include "Rodin/Variational/Jacobian.h"

#include "../CellDeformation.h"
#include "Parameters.h"

namespace Rodin::Adaptation::WNGIR
{
  /// @brief Sampled geometric admissibility diagnostics.
  struct AdmissibilityReport
  {
      /// @brief Minimum sampled Jacobian determinant.
      Real minJ = std::numeric_limits<Real>::infinity();
      /// @brief Number of failed Jacobian or distortion validity checks.
      std::size_t inadmissibleCount = 0;
      /// @brief Maximum finite sampled relative distortion.
      Real maxQRel = Real(0);
  };

  template <class Displacement>
  /**
   * @brief Evaluates sampled admissibility without changing the displacement.
   * @param u Displacement field to sample.
   * @param jacobian Lower admissible relative Jacobian bound.
   * @param quadratureOrder Sampling order; zero selects the finite-element rule.
   * @returns Sampled Jacobian, relative distortion and invalid-sample count.
   */
  AdmissibilityReport evaluateAdmissibility(const Displacement& u,
    Real jacobian, std::size_t quadratureOrder = 0)
  {
    using Variational::IntegrationPoint;
    using Variational::Jacobian;

    AdmissibilityReport rep;
    const auto& fes = u.getFiniteElementSpace();
    const auto& mesh = fes.getMesh();
    const std::size_t dim = mesh.getDimension();
    const std::size_t vdim = fes.getVectorDimension();
    if (dim != vdim)
      Alert::Exception() << "WNGIR displacement dimension must equal mesh dimension."
                         << Alert::Raise;

    auto gradU = Jacobian(u);
    CellDeformation deformation(dim);
    for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
    {
      const auto& cell = *cellIt;
      const auto& fe = fes.getFiniteElement(cell.getDimension(), cell.getIndex());
      const auto& qf = QF::PolytopeQuadratureFormula::get(quadratureOrder > 0
          ? quadratureOrder
          : Parameters::Quadrature::getCellOrder(fe.getOrder()),
        cell.getGeometry());
      const auto& quadrature = cell.getQuadrature(qf);
      for (std::size_t q = 0; q < quadrature.getSize(); ++q)
      {
        const auto& pt = quadrature.getPoint(q);
        const IntegrationPoint ip(pt, &qf, q);
        deformation.setDisplacementGradient(gradU.getValue(ip));
        const Real j = deformation.getJacobian();

        rep.minJ = std::min(rep.minJ, j);
        bool invalid = !std::isfinite(j) || j <= jacobian || !deformation.isAdmissible();

        if (deformation.isAdmissible())
        {
          const Real qRel = deformation.getRelativeDistortion();
          if (!std::isfinite(qRel))
            invalid = true;
          else
            rep.maxQRel = std::max(rep.maxQRel, qRel);
        }
        if (invalid)
          ++rep.inadmissibleCount;
      }
    }
    return rep;
  }
}

#endif
