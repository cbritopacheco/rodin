/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_CORONARYARTERY_GRADEDYEOH_H
#define RODIN_CORONARYARTERY_GRADEDYEOH_H

/**
 * @file GradedYeoh.h
 * @brief Position-graded compressible Yeoh hyperelastic law.
 */

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
   * The multiplier is sampled on every binding and remains fixed for all
   * subsequent const evaluations until the next binding.
   */
  class GradedYeoh final : public Solid::HyperElasticLaw<GradedYeoh>
  {
    public:
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

      /**
       * @brief Captures kinematics and samples the spatial multiplier on every bind.
       * @param cp Incoming constitutive point.
       * @returns This prepared law.
       */
      GradedYeoh& setConstitutivePoint(const Solid::ConstitutivePoint& cp)
      {
        m_set = false;
        const Real scale = factor(cp);
        m_base.setConstitutivePoint(cp);
        m_factor = scale;
        m_set = true;
        return *this;
      }

      /// @brief Returns owned kinematics of the prepared base law.
      const Solid::KinematicState& getKinematicState() const
      {
        this->require(m_set, __func__);
        return m_base.getKinematicState();
      }

      /// @brief Returns the graded strain-energy density.
      Real getStrainEnergyDensity() const
      {
        this->require(m_set, __func__);
        return m_factor * m_base.getStrainEnergyDensity();
      }

      /// @brief Computes the graded first Piola-Kirchhoff stress.
      void getFirstPiolaKirchhoffStress(Math::SpatialMatrix<Real>& P) const
      {
        this->require(m_set, __func__);
        m_base.getFirstPiolaKirchhoffStress(P);
        P *= m_factor;
      }

      /// @brief Computes the graded material tangent action.
      void getMaterialTangent(
        Math::SpatialMatrix<Real>& dP, const Math::SpatialMatrix<Real>& dF) const
      {
        this->require(m_set, __func__);
        m_base.getMaterialTangent(dP, dF);
        dP *= m_factor;
      }

      /// @brief Gets the base (ungraded) Yeoh law.
      const Solid::Yeoh& getBaseLaw() const
      {
        return m_base;
      }

    private:
      /// @brief Captured spatial multiplier for the current binding.
      Real m_factor = 1;

      /// @brief Whether base preparation and multiplier evaluation succeeded.
      bool m_set = false;

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
