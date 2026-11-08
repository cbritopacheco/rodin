/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_VARIATIONAL_HYPERELASTICLAW_H
#define RODIN_KOKKOS_VARIATIONAL_HYPERELASTICLAW_H

/**
 * @file HyperElasticLaw.h
 * @brief Device hyperelastic laws of constant parameters.
 *
 * Each device law evaluates, at a displacement gradient @f$ H @f$ and
 * @f$ F = I + H @f$, the first Piola-Kirchhoff stress @f$ P(F) @f$ and the
 * elasticity tensor
 * @f[
 *   \mathbb{A}_{iJkL} = \frac{\partial P_{iJ}}{\partial F_{kL}},
 * @f]
 * whose action on @f$ dF @f$ is the host
 * HyperElasticLaw::getMaterialTangent(). The device law of a host law is
 * given by DeviceLaw.
 */

#include <concepts>

#include <Kokkos_Core.hpp>

#include "Rodin/Solid/Constitutive/MooneyRivlin.h"
#include "Rodin/Solid/Constitutive/SaintVenantKirchhoff.h"
#include "Rodin/Solid/Constitutive/Yeoh.h"

#include "Rodin/Kokkos/Geometry/MeshView.h"

namespace Rodin::Kokkos::Variational
{
  /// @brief Stride of the indices of an elasticity tensor.
  constexpr Index ElasticityStride = RODIN_MAXIMAL_SPACE_DIMENSION;

  /**
   * @brief Number of entries of a row @f$ \mathbb{A}_{i \cdot \cdot \cdot} @f$
   * of an elasticity tensor.
   *
   * Entry @f$ \mathbb{A}_{iJkL} @f$ of row i is at
   * @f$ (J s + k) s + L @f$ for the stride s of ElasticityStride.
   */
  constexpr Index ElasticityRow = ElasticityStride * ElasticityStride * ElasticityStride;

  /**
   * @brief Position of an entry in a row of an elasticity tensor.
   * @param J Column of the stress.
   * @param k Row of the deformation gradient.
   * @param L Column of the deformation gradient.
   * @return Offset of @f$ \mathbb{A}_{iJkL} @f$ in row i.
   */
  KOKKOS_INLINE_FUNCTION
  constexpr Index getElasticityIndex(Index J, Index k, Index L)
  {
    return (J * ElasticityStride + k) * ElasticityStride + L;
  }

  /**
   * @brief Kinematics of the isochoric invariants.
   *
   * Holds @f$ F @f$, @f$ F^{-T} @f$, @f$ J @f$,
   * @f$ I_1 = F : F @f$ and @f$ J^{-2/d} @f$, as the host Mooney-Rivlin and
   * Yeoh laws.
   */
  struct IsochoricPoint
  {
    /// @param h Displacement gradient @f$ H @f$, a square matrix.
    KOKKOS_INLINE_FUNCTION
    explicit IsochoricPoint(const Geometry::Matrix& h)
      : F(h),
        FinvT(h.rows(), h.cols())
    {
      KOKKOS_ASSERT(h.rows() == h.cols());
      const Index d = F.rows();
      for (Index i = 0; i < d; ++i)
        F(i, i) += 1;
      J = F.determinant();
      KOKKOS_ASSERT(J > 0);
      const Geometry::Matrix Finv = F.inverse();
      I1 = 0;
      for (Index i = 0; i < d; ++i)
      {
        for (Index j = 0; j < d; ++j)
        {
          FinvT(i, j) = Finv(j, i);
          I1 += F(i, j) * F(i, j);
        }
      }
      Jm2d = ::Kokkos::pow(J, -Real(2) / Real(d));
    }

    /**
     * @brief Gets @f$ \partial \bar I_1 / \partial F @f$.
     * @param i Row.
     * @param j Column.
     * @return @f$ J^{-2/d} (2 F - \frac{2}{d} I_1 F^{-T})_{ij} @f$.
     */
    KOKKOS_INLINE_FUNCTION
    Real getIsochoricGradient(Index i, Index j) const
    {
      const Real d = Real(F.rows());
      return Jm2d * (2 * F(i, j) - (2 / d) * I1 * FinvT(i, j));
    }

    /**
     * @brief Gets the derivative of @f$ \partial \bar I_1 / \partial F @f$.
     * @param i Row of the first deformation gradient entry.
     * @param J Column of the first deformation gradient entry.
     * @param k Row of the second deformation gradient entry.
     * @param L Column of the second deformation gradient entry.
     * @return @f$ \partial^2 \bar I_1 / \partial F_{iJ} \partial F_{kL} @f$.
     */
    KOKKOS_INLINE_FUNCTION
    Real getIsochoricHessian(Index i, Index J, Index k, Index L) const
    {
      const Real d = Real(F.rows());
      Real v =
        -(4 / d) * F(k, L) * FinvT(i, J) + (2 / d) * I1 * FinvT(i, L) * FinvT(k, J);
      if (i == k && J == L)
        v += 2;
      return -(2 / d) * FinvT(k, L) * getIsochoricGradient(i, J) + Jm2d * v;
    }

    /**
     * @brief Gets the elasticity of @f$ \frac{\kappa}{2} (J - 1)^2 @f$.
     * @param kappa Bulk modulus.
     * @param i Row of the first deformation gradient entry.
     * @param j Column of the first deformation gradient entry.
     * @param k Row of the second deformation gradient entry.
     * @param l Column of the second deformation gradient entry.
     * @return Entry @f$ iJkL @f$ of the volumetric elasticity.
     */
    KOKKOS_INLINE_FUNCTION
    Real getVolumetricElasticity(Real kappa, Index i, Index j, Index k, Index l) const
    {
      return kappa *
        ((2 * J - 1) * J * FinvT(k, l) * FinvT(i, j) -
          (J - 1) * J * FinvT(i, l) * FinvT(k, j));
    }

    /// @brief Deformation gradient @f$ F @f$.
    Geometry::Matrix F;

    /// @brief Inverse transpose @f$ F^{-T} @f$.
    Geometry::Matrix FinvT;

    /// @brief Jacobian @f$ J = \det F @f$.
    Real J = 1;

    /// @brief First invariant @f$ I_1 = \operatorname{tr} C @f$.
    Real I1 = 0;

    /// @brief Isochoric scaling @f$ J^{-2/d} @f$.
    Real Jm2d = 1;
  };

  /**
   * @brief Device Saint Venant-Kirchhoff law.
   *
   * With @f$ E = \frac{1}{2} (F^T F - I) @f$ and
   * @f$ S = \lambda \operatorname{tr} E \ I + 2 \mu E @f$, the stress is
   * @f$ P = F S @f$ and
   * @f[
   *   \mathbb{A}_{iJkL} = \delta_{ik} S_{LJ} + \lambda F_{iJ} F_{kL}
   *   + \mu (F_{iL} F_{kJ} + (F F^T)_{ik} \delta_{JL}),
   * @f]
   * as Solid::SaintVenantKirchhoff.
   */
  struct DeviceSaintVenantKirchhoff
  {
    /// @brief Kinematics at a quadrature point.
    struct Point
    {
      /// @brief Deformation gradient @f$ F @f$.
      Geometry::Matrix F;

      /// @brief Second Piola-Kirchhoff stress @f$ S @f$.
      Geometry::Matrix S;

      /// @brief Left Cauchy-Green tensor @f$ F F^T @f$.
      Geometry::Matrix B;
    };

    /**
     * @param h Displacement gradient @f$ H @f$, a square matrix.
     * @return Kinematics at @f$ F = I + H @f$.
     */
    KOKKOS_INLINE_FUNCTION
    Point getPoint(const Geometry::Matrix& h) const
    {
      KOKKOS_ASSERT(h.rows() == h.cols());
      const Index d = h.rows();
      Point p{h, Geometry::Matrix(d, d), Geometry::Matrix(d, d)};
      for (Index i = 0; i < d; ++i)
        p.F(i, i) += 1;
      Real trE = 0;
      for (Index m = 0; m < d; ++m)
      {
        for (Index n = 0; n < d; ++n)
        {
          Real c = 0;
          Real b = 0;
          for (Index a = 0; a < d; ++a)
          {
            c += p.F(a, m) * p.F(a, n);
            b += p.F(m, a) * p.F(n, a);
          }
          p.S(m, n) = Real(0.5) * (c - (m == n ? 1 : 0));
          p.B(m, n) = b;
        }
        trE += p.S(m, m);
      }
      for (Index m = 0; m < d; ++m)
      {
        for (Index n = 0; n < d; ++n)
          p.S(m, n) *= 2 * mu;
        p.S(m, m) += lambda * trE;
      }
      return p;
    }

    /**
     * @param p Kinematics.
     * @param[out] P First Piola-Kirchhoff stress @f$ F S @f$.
     */
    KOKKOS_INLINE_FUNCTION
    void getStress(const Point& p, Geometry::Matrix& P) const
    {
      const Index d = p.F.rows();
      P = Geometry::Matrix(d, d);
      for (Index i = 0; i < d; ++i)
      {
        for (Index j = 0; j < d; ++j)
        {
          Real v = 0;
          for (Index m = 0; m < d; ++m)
            v += p.F(i, m) * p.S(m, j);
          P(i, j) = v;
        }
      }
    }

    /**
     * @param p Kinematics.
     * @param i Row of the stress.
     * @param scale Factor of the entries.
     * @param[out] row Entries of @f$ \text{scale} \ \mathbb{A}_{i \cdot \cdot \cdot} @f$
     * at getElasticityIndex().
     */
    KOKKOS_INLINE_FUNCTION
    void getElasticity(const Point& p, Index i, Real scale, Real* row) const
    {
      const Index d = p.F.rows();
      for (Index J = 0; J < d; ++J)
      {
        for (Index k = 0; k < d; ++k)
        {
          for (Index L = 0; L < d; ++L)
          {
            Real v = lambda * p.F(i, J) * p.F(k, L) + mu * p.F(i, L) * p.F(k, J);
            if (i == k)
              v += p.S(L, J);
            if (J == L)
              v += mu * p.B(i, k);
            row[getElasticityIndex(J, k, L)] = scale * v;
          }
        }
      }
    }

    /// @brief First Lamé parameter @f$ \lambda @f$.
    Real lambda;

    /// @brief Shear modulus @f$ \mu @f$.
    Real mu;
  };

  /**
   * @brief Device Yeoh law.
   *
   * With @f$ x = \bar I_1 - d @f$, the energy
   * @f$ c_1 x + c_2 x^2 + c_3 x^3 + \frac{\kappa}{2} (J - 1)^2 @f$ gives
   * @f$ P = W_1 \partial \bar I_1 / \partial F + \kappa (J - 1) J F^{-T} @f$,
   * @f$ W_1 = c_1 + 2 c_2 x + 3 c_3 x^2 @f$, as Solid::Yeoh.
   */
  struct DeviceYeoh
  {
    /// @brief Kinematics at a quadrature point.
    using Point = IsochoricPoint;

    /**
     * @param h Displacement gradient @f$ H @f$, a square matrix.
     * @return Kinematics at @f$ F = I + H @f$.
     */
    KOKKOS_INLINE_FUNCTION
    Point getPoint(const Geometry::Matrix& h) const
    {
      return Point(h);
    }

    /**
     * @param p Kinematics.
     * @param[out] P First Piola-Kirchhoff stress.
     */
    KOKKOS_INLINE_FUNCTION
    void getStress(const Point& p, Geometry::Matrix& P) const
    {
      const Index d = p.F.rows();
      const Real x = p.Jm2d * p.I1 - Real(d);
      const Real w1 = c1 + 2 * c2 * x + 3 * c3 * x * x;
      const Real volume = kappa * (p.J - 1) * p.J;
      P = Geometry::Matrix(d, d);
      for (Index i = 0; i < d; ++i)
      {
        for (Index j = 0; j < d; ++j)
          P(i, j) = w1 * p.getIsochoricGradient(i, j) + volume * p.FinvT(i, j);
      }
    }

    /**
     * @param p Kinematics.
     * @param i Row of the stress.
     * @param scale Factor of the entries.
     * @param[out] row Entries of @f$ \text{scale} \ \mathbb{A}_{i \cdot \cdot \cdot} @f$
     * at getElasticityIndex().
     */
    KOKKOS_INLINE_FUNCTION
    void getElasticity(const Point& p, Index i, Real scale, Real* row) const
    {
      const Index d = p.F.rows();
      const Real x = p.Jm2d * p.I1 - Real(d);
      const Real w1 = c1 + 2 * c2 * x + 3 * c3 * x * x;
      const Real w2 = 2 * c2 + 6 * c3 * x;
      for (Index J = 0; J < d; ++J)
      {
        const Real gi = p.getIsochoricGradient(i, J);
        for (Index k = 0; k < d; ++k)
        {
          for (Index L = 0; L < d; ++L)
          {
            const Real v = w2 * gi * p.getIsochoricGradient(k, L) +
              w1 * p.getIsochoricHessian(i, J, k, L) +
              p.getVolumetricElasticity(kappa, i, J, k, L);
            row[getElasticityIndex(J, k, L)] = scale * v;
          }
        }
      }
    }

    /// @brief Material constant @f$ c_1 @f$.
    Real c1;

    /// @brief Material constant @f$ c_2 @f$.
    Real c2;

    /// @brief Material constant @f$ c_3 @f$.
    Real c3;

    /// @brief Bulk modulus @f$ \kappa @f$.
    Real kappa;
  };

  /**
   * @brief Device Mooney-Rivlin law.
   *
   * The energy
   * @f$ c_1 (\bar I_1 - d) + c_2 (\bar I_2 - d) + \frac{\kappa}{2} (J - 1)^2 @f$,
   * @f$ \bar I_2 = J^{-4/d} I_2 @f$, gives
   * @f[
   *   P = c_1 \frac{\partial \bar I_1}{\partial F}
   *   + c_2 J^{-4/d} (2 I_1 F - 2 F C - \tfrac{4}{d} I_2 F^{-T})
   *   + \kappa (J - 1) J F^{-T},
   * @f]
   * as Solid::MooneyRivlin.
   */
  struct DeviceMooneyRivlin
  {
    /// @brief Kinematics at a quadrature point.
    struct Point : IsochoricPoint
    {
      /// @param h Displacement gradient @f$ H @f$, a square matrix.
      KOKKOS_INLINE_FUNCTION
      explicit Point(const Geometry::Matrix& h)
        : IsochoricPoint(h),
          C(h.rows(), h.cols()),
          FC(h.rows(), h.cols()),
          B(h.rows(), h.cols())
      {
        const Index d = F.rows();
        for (Index m = 0; m < d; ++m)
        {
          for (Index n = 0; n < d; ++n)
          {
            Real c = 0;
            Real b = 0;
            for (Index a = 0; a < d; ++a)
            {
              c += F(a, m) * F(a, n);
              b += F(m, a) * F(n, a);
            }
            C(m, n) = c;
            B(m, n) = b;
          }
        }
        Real trC2 = 0;
        for (Index m = 0; m < d; ++m)
        {
          for (Index n = 0; n < d; ++n)
          {
            Real v = 0;
            for (Index a = 0; a < d; ++a)
              v += F(m, a) * C(a, n);
            FC(m, n) = v;
            trC2 += C(m, n) * C(n, m);
          }
        }
        I2 = Real(0.5) * (I1 * I1 - trC2);
        Jm4d = ::Kokkos::pow(J, -Real(4) / Real(d));
      }

      /**
       * @brief Gets @f$ \partial \bar I_2 / \partial F @f$.
       * @param i Row.
       * @param j Column.
       * @return @f$ J^{-4/d} (2 I_1 F - 2 F C - \frac{4}{d} I_2 F^{-T})_{ij} @f$.
       */
      KOKKOS_INLINE_FUNCTION
      Real getSecondGradient(Index i, Index j) const
      {
        const Real d = Real(F.rows());
        return Jm4d * (2 * I1 * F(i, j) - 2 * FC(i, j) - (4 / d) * I2 * FinvT(i, j));
      }

      /**
       * @brief Gets the derivative of @f$ \partial \bar I_2 / \partial F @f$.
       * @param i Row of the first deformation gradient entry.
       * @param J Column of the first deformation gradient entry.
       * @param k Row of the second deformation gradient entry.
       * @param L Column of the second deformation gradient entry.
       * @return @f$ \partial^2 \bar I_2 / \partial F_{iJ} \partial F_{kL} @f$.
       */
      KOKKOS_INLINE_FUNCTION
      Real getSecondHessian(Index i, Index J, Index k, Index L) const
      {
        const Real d = Real(F.rows());
        const Real dI2 = 2 * I1 * F(k, L) - 2 * FC(k, L);
        Real v = 4 * F(k, L) * F(i, J) - 2 * F(i, L) * F(k, J) -
          (4 / d) * dI2 * FinvT(i, J) + (4 / d) * I2 * FinvT(i, L) * FinvT(k, J);
        if (i == k)
          v -= 2 * C(L, J);
        if (J == L)
          v -= 2 * B(i, k);
        if (i == k && J == L)
          v += 2 * I1;
        return -(4 / d) * FinvT(k, L) * getSecondGradient(i, J) + Jm4d * v;
      }

      /// @brief Right Cauchy-Green tensor @f$ C = F^T F @f$.
      Geometry::Matrix C;

      /// @brief Product @f$ F C @f$.
      Geometry::Matrix FC;

      /// @brief Left Cauchy-Green tensor @f$ F F^T @f$.
      Geometry::Matrix B;

      /// @brief Second invariant @f$ I_2 @f$.
      Real I2 = 0;

      /// @brief Isochoric scaling @f$ J^{-4/d} @f$.
      Real Jm4d = 1;
    };

    /**
     * @param h Displacement gradient @f$ H @f$, a square matrix.
     * @return Kinematics at @f$ F = I + H @f$.
     */
    KOKKOS_INLINE_FUNCTION
    Point getPoint(const Geometry::Matrix& h) const
    {
      return Point(h);
    }

    /**
     * @param p Kinematics.
     * @param[out] P First Piola-Kirchhoff stress.
     */
    KOKKOS_INLINE_FUNCTION
    void getStress(const Point& p, Geometry::Matrix& P) const
    {
      const Index d = p.F.rows();
      const Real volume = kappa * (p.J - 1) * p.J;
      P = Geometry::Matrix(d, d);
      for (Index i = 0; i < d; ++i)
      {
        for (Index j = 0; j < d; ++j)
        {
          P(i, j) = c1 * p.getIsochoricGradient(i, j) + c2 * p.getSecondGradient(i, j) +
            volume * p.FinvT(i, j);
        }
      }
    }

    /**
     * @param p Kinematics.
     * @param i Row of the stress.
     * @param scale Factor of the entries.
     * @param[out] row Entries of @f$ \text{scale} \ \mathbb{A}_{i \cdot \cdot \cdot} @f$
     * at getElasticityIndex().
     */
    KOKKOS_INLINE_FUNCTION
    void getElasticity(const Point& p, Index i, Real scale, Real* row) const
    {
      const Index d = p.F.rows();
      for (Index J = 0; J < d; ++J)
      {
        for (Index k = 0; k < d; ++k)
        {
          for (Index L = 0; L < d; ++L)
          {
            const Real v = c1 * p.getIsochoricHessian(i, J, k, L) +
              c2 * p.getSecondHessian(i, J, k, L) +
              p.getVolumetricElasticity(kappa, i, J, k, L);
            row[getElasticityIndex(J, k, L)] = scale * v;
          }
        }
      }
    }

    /// @brief Material constant @f$ c_1 @f$.
    Real c1;

    /// @brief Material constant @f$ c_2 @f$.
    Real c2;

    /// @brief Bulk modulus @f$ \kappa @f$.
    Real kappa;
  };

  /**
   * @brief Device law of a host hyperelastic law.
   *
   * Specializations define Type, a trivially copyable device law, and
   * get(), building it from the host law. Laws without a specialization,
   * among them the laws reading constitutive inputs, stay on the host.
   * Solid::NeoHookean has its own device integrators.
   */
  template <class Law>
  struct DeviceLaw
  {};

  /// @brief Device law of Solid::SaintVenantKirchhoff.
  template <>
  struct DeviceLaw<Solid::SaintVenantKirchhoff>
  {
    /// @brief Device law type.
    using Type = DeviceSaintVenantKirchhoff;

    /// @param law Host law. @return Device law.
    static Type get(const Solid::SaintVenantKirchhoff& law)
    {
      return {law.getLameFirstParameter(), law.getShearModulus()};
    }
  };

  /// @brief Device law of Solid::Yeoh.
  template <>
  struct DeviceLaw<Solid::Yeoh>
  {
    /// @brief Device law type.
    using Type = DeviceYeoh;

    /// @param law Host law. @return Device law.
    static Type get(const Solid::Yeoh& law)
    {
      return {law.getMaterialConstantC1(), law.getMaterialConstantC2(),
        law.getMaterialConstantC3(), law.getBulkModulus()};
    }
  };

  /// @brief Device law of Solid::MooneyRivlin.
  template <>
  struct DeviceLaw<Solid::MooneyRivlin>
  {
    /// @brief Device law type.
    using Type = DeviceMooneyRivlin;

    /// @param law Host law. @return Device law.
    static Type get(const Solid::MooneyRivlin& law)
    {
      return {
        law.getMaterialConstantC1(), law.getMaterialConstantC2(), law.getBulkModulus()};
    }
  };

  /// @brief Host hyperelastic laws with a device law in DeviceLaw.
  template <class Law>
  concept DeviceHyperElasticLaw = requires(const Law& law) {
    { DeviceLaw<Law>::get(law) } -> std::same_as<typename DeviceLaw<Law>::Type>;
  };
}

#endif
