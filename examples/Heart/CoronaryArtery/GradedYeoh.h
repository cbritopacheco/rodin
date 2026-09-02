/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef EXAMPLES_HEART_CORONARYARTERY_GRADEDYEOH_H
#define EXAMPLES_HEART_CORONARYARTERY_GRADEDYEOH_H

#include <functional>
#include <utility>

#include <Rodin/Solid.h>
#include <Rodin/Types.h>

namespace Rodin::Examples::Heart
{
  /**
   * @brief Position-graded compressible Yeoh hyperelastic law.
   *
   * Wraps a base @ref Rodin::Solid::Yeoh law with a spatial multiplier
   * @f$ m(x) @f$ evaluated at the quadrature point:
   *
   * @f[
   *   W_{\mathrm{graded}}(x, \mathbf{F}) = m(x) \, W_{\mathrm{Yeoh}}(\mathbf{F}).
   * @f]
   *
   * The Yeoh strain energy is JOINTLY LINEAR in its material constants
   * @f$ (c_1, c_2, c_3, \kappa) @f$, so scaling the energy by @f$ m(x) @f$
   * scales the first Piola-Kirchhoff stress and the material tangent by
   * exactly the same factor.  The wrapper therefore delegates every
   * evaluation to the base law and multiplies the OUTPUT -- no constitutive
   * re-derivation, and the Newton tangent remains exactly consistent with
   * the residual.
   *
   * The multiplier receives the Geometry::Point carried by the
   * ConstitutivePoint (the integrators construct it with full geometric
   * context), so any position-dependent grading works: through-thickness
   * intima/media/adventitia profiles, axial tapering, plaque patches, etc.
   * If the ConstitutivePoint carries no geometric context, or no scale
   * function was given, the factor defaults to 1 (plain Yeoh).
   */
  class GradedYeoh final : public Solid::HyperElasticLaw<GradedYeoh>
  {
    public:
      /// Same precomputed cache as the base law.
      using Cache = Solid::Yeoh::Cache;

      /// Spatial multiplier evaluated at the quadrature point.
      using ScaleFunction = std::function<Real(const Geometry::Point&)>;

      /**
       * @brief Constructs a graded Yeoh law.
       * @param base  Base Yeoh law (media-level constants).
       * @param scale Spatial multiplier m(x); empty means m = 1 everywhere.
       */
      GradedYeoh(Solid::Yeoh base, ScaleFunction scale)
        : m_base(std::move(base)),
          m_scale(std::move(scale))
      {}

      /// @brief Copy constructor.
      GradedYeoh(const GradedYeoh&) = default;
      /// @brief Move constructor.
      GradedYeoh(GradedYeoh&&) = default;

      /// @brief Populates the constitutive cache (delegates to the base law).
      void setCache(Cache& cache, const Solid::ConstitutivePoint& cp) const
      {
        m_base.setCache(cache, cp);
      }

      /// @brief Graded strain-energy density m(x) W(F).
      Real getStrainEnergyDensity(
        const Cache& cache, const Solid::ConstitutivePoint& cp) const
      {
        return factor(cp) * m_base.getStrainEnergyDensity(cache, cp);
      }

      /// @brief Graded first Piola-Kirchhoff stress m(x) P(F).
      void getFirstPiolaKirchhoffStress(Math::SpatialMatrix<Real>& P, const Cache& cache,
        const Solid::ConstitutivePoint& cp) const
      {
        m_base.getFirstPiolaKirchhoffStress(P, cache, cp);
        P *= factor(cp);
      }

      /// @brief Graded material tangent m(x) dP(F)[dF].
      void getMaterialTangent(Math::SpatialMatrix<Real>& dP, const Cache& cache,
        const Solid::ConstitutivePoint& cp, const Math::SpatialMatrix<Real>& dF) const
      {
        m_base.getMaterialTangent(dP, cache, cp, dF);
        dP *= factor(cp);
      }

      /// @brief Gets the base (ungraded) Yeoh law.
      const Solid::Yeoh& getBaseLaw() const
      {
        return m_base;
      }

    private:
      Real factor(const Solid::ConstitutivePoint& cp) const
      {
        if (!m_scale)
          return 1.0;
        const auto& pt = cp.getPoint();
        if (!pt)
          return 1.0;
        return m_scale(pt->get());
      }

      Solid::Yeoh m_base;
      ScaleFunction m_scale;
  };
}

#endif
