/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_INTEGRATIONPOINT_H
#define RODIN_VARIATIONAL_INTEGRATIONPOINT_H

#include <cassert>

#include "Rodin/Types.h"
#include "Rodin/Geometry/Point.h"
#include "Rodin/Geometry/PolytopeQuadrature.h"
#include "Rodin/QF/QuadratureFormula.h"

namespace Rodin::Variational
{
  /**
   * @brief Evaluation point used by variational operators.
   *
   * An IntegrationPoint borrows a geometric point and optionally records the
   * formula/index of its originating integration loop. This sampling metadata
   * does not certify that the point still has that formula's reference
   * coordinates: mapping a face sample into a cell changes the chart.
   *
   * ## Reference association
   *
   * Construction from @ref Geometry::PolytopeQuadrature and a logical index
   * establishes reference association through owned samples, without coordinate
   * comparisons. @ref getReferenceQuadratureFormula returns that formula only
   * while its lifetime/assignment identity matches the quadrature's snapshot.
   * Direct point/formula/index construction records provenance only. Reference
   * tabulation and coordinate-sensitive reuse must use the reference accessor;
   * an absent association selects the supplied point's reference coordinates.
   * Coefficient snapshots may still use the originating sampling metadata.
   *
   * @note The quadrature formula pointer is non-owning. The referenced
   * quadrature object must outlive all IntegrationPoint uses that dereference
   * this pointer.
   */
  class IntegrationPoint
  {
    public:
      /**
       * @brief Constructs an integration point from a geometric point.
       * @param p Point at which the operation is evaluated.
       */
    IntegrationPoint(const Geometry::Point& p)
      : m_p(p),
        m_qf(nullptr),
        m_qp(0),
        m_referenceIdentity(std::nullopt)
    {}

      /**
       * @brief Records sampling provenance without certifying reference coordinates.
       * @param[in] p Geometric point, possibly mapped from the quadrature sample into another polytope
       * @param[in] qf Non-owning pointer to the quadrature formula that owns
       *                the sample indexing, or @c nullptr for direct
       *                pointwise evaluation
       * @param[in] qp Quadrature sample index in @p qf
       */
    IntegrationPoint(
      const Geometry::Point& p, const QF::QuadratureFormulaBase* qf, size_t qp)
      : m_p(p),
        m_qf(qf),
        m_qp(qp),
        m_referenceIdentity(std::nullopt)
    {}

      /**
       * @brief Borrows an owned mapped sample with logical reference association.
       *
       * No reference to the quadrature container is retained. The point and
       * formula must outlive this binding, including after a container move.
       * Formula reassignment disables reference reuse without changing the
       * owned point or its originating metadata. Geometry changes still require
       * normal mapped-point reconstruction and rebinding.
       *
       * @param quadrature Owner of the mapped point and its reference association.
       * @param qp Logical sample index in the mapped quadrature.
       */
    IntegrationPoint(const Geometry::PolytopeQuadrature& quadrature, size_t qp)
      : m_p(quadrature.getPoint(qp)),
        m_qf(&quadrature.getQuadratureFormula()),
        m_qp(qp),
        m_referenceIdentity(quadrature.getQuadratureFormulaIdentity())
    {}

      /// @returns Underlying geometric point context.
    const Geometry::Point& getPoint() const
    {
      return m_p;
    }

      /**
       * @returns Non-owning pointer to the quadrature formula, or @c nullptr
       * when this integration point is not in a quadrature context.
       */
    const QF::QuadratureFormulaBase* getQuadratureFormula() const
    {
      return m_qf;
    }

      /**
       * @brief Returns the logically certified reference formula, when available.
       *
       * The comparison uses formula lifetime/assignment identity, never point
       * coordinates. Reads are local and non-collective and create no identity.
       *
       * @returns Formula authorizing reference tabulation, or @c nullptr.
       */
    const QF::QuadratureFormulaBase* getReferenceQuadratureFormula() const noexcept
    {
      return m_referenceIdentity && *m_referenceIdentity == m_qf->getCacheIdentity()
        ? m_qf
        : nullptr;
    }

      /**
       * @returns Quadrature index associated with this point.
       * @note Meaningful only when getQuadratureFormula() is not @c nullptr.
       */
    size_t getIndex() const
    {
      return m_qp;
    }

    private:
    std::reference_wrapper<const Geometry::Point> m_p;
    const QF::QuadratureFormulaBase* m_qf;
    size_t m_qp;
    Optional<size_t> m_referenceIdentity;
  };
}

#endif
