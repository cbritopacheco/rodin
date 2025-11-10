/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file FaceIntegral.h
 * @brief Face integral classes for discontinuous Galerkin methods.
 *
 * This file defines the FaceIntegral classes which represent integrals over
 * mesh faces (element interfaces) in discontinuous Galerkin (DG) formulations.
 * Face integrals are crucial for coupling elements in DG methods where
 * solutions are discontinuous across element boundaries.
 *
 * ## Mathematical Foundation
 * Face integrals compute:
 * @f[
 *   \int_{\mathcal{F}_h} f(x) \, ds
 * @f]
 * where @f$ \mathcal{F}_h @f$ is the set of all mesh faces (both interior
 * and boundary faces).
 *
 * ## Discontinuous Galerkin Context
 * In DG methods, face integrals handle:
 * - **Numerical fluxes**: Coupling between elements
 * - **Jump terms**: @f$ \int_{\mathcal{F}_h} [\![u]\!] \cdot [\![v]\!] \, ds @f$
 * - **Average terms**: @f$ \int_{\mathcal{F}_h} \{\!\{u\}\!\} \cdot [\![v]\!] \, ds @f$
 * - **Penalty terms**: Interior penalty for stability
 *
 * ## Usage Example
 * ```cpp
 * // Interior penalty DG method
 * auto penalty = FaceIntegral(Jump(u), Jump(v));
 * auto consistency = FaceIntegral(Average(Grad(u)), Jump(v));
 * auto symmetry = FaceIntegral(Jump(u), Average(Grad(v)));
 * ```
 *
 * @see InterfaceIntegral, BoundaryIntegral, Jump, Average
 */
#ifndef RODIN_VARIATIONAL_FACEINTEGRAL_H
#define RODIN_VARIATIONAL_FACEINTEGRAL_H

#include <memory>
#include <utility>

#include "Rodin/Geometry/Region.h"

#include "ForwardDecls.h"

namespace Rodin::Variational
{
  /**
   * @defgroup FaceIntegralSpecializations FaceIntegral Template Specializations
   * @brief Template specializations of the FaceIntegral class.
   *
   * @see FaceIntegral
   */

  /**
   * @ingroup FaceIntegralSpecializations
   * @brief Integration of the dot product of a trial and test operators.
   *
   * Given two operators defined over trial and test spaces @f$ U_h
   * @f$ and @f$ V_h @f$,
   * @f[
   *   A : U_h \rightarrow \mathbb{R}^{p \times q}, \quad B : V_h \rightarrow \mathbb{R}^{p \times q},
   * @f]
   * this class represents the integral of their dot product:
   * @f[
   *   \int_{\mathcal{F}_h} A(u) : B(v) \ d\sigma(x) \ .
   * @f]
   */
  template <class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  class FaceIntegral<Dot<
          ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>,
          ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>> final
    : public QuadratureRule<Dot<
          ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>,
          ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>
  {
    public:
      using LHSType = ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>;

      using RHSType = ShapeFunctionBase<RHSDerived, TestFES, TestSpace>;

      using IntegrandType = Dot<LHSType, RHSType>;

      using Parent = QuadratureRule<IntegrandType>;

      /**
       * @brief Constructs face integral from trial and test operators.
       * @param[in] lhs Trial operator @f$ A(u) @f$
       * @param[in] rhs Test operator @f$ B(v) @f$
       *
       * Represents the face integral over all mesh faces:
       * @f[
       *   \int_{\mathcal{F}_h} A(u) : B(v) \, d\sigma
       * @f]
       * where @f$ \mathcal{F}_h @f$ denotes all faces in the mesh.
       */
      FaceIntegral(const LHSType& lhs, const RHSType& rhs)
        : FaceIntegral(Dot(lhs, rhs))
      {}

      /**
       * @brief Constructs face integral from dot product.
       * @param[in] prod Dot product of trial and test operators
       */
      FaceIntegral(const IntegrandType& prod)
        : Parent(prod)
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Face integral to copy
       */
      FaceIntegral(const FaceIntegral& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Face integral to move from
       */
      FaceIntegral(FaceIntegral&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Gets the integration region.
       * @return Region::Faces indicating integration over all faces
       *
       * Face integrals are computed over all faces @f$ F \in \mathcal{F}_h @f$,
       * including both interior faces and boundary faces. Used extensively in
       * discontinuous Galerkin (DG) methods.
       */
      Geometry::Region getRegion() const override
      {
        return Geometry::Region::Faces;
      }

      /**
       * @brief Creates a polymorphic copy of this face integral.
       * @return Pointer to a new copy
       */
      FaceIntegral* copy() const noexcept override
      {
        return new FaceIntegral(*this);
      }

    private:
      std::unique_ptr<IntegrandType> m_prod;
  };

  template <class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  FaceIntegral(const Dot<ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>, ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>&)
    -> FaceIntegral<Dot<ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>, ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>;

  template <class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  FaceIntegral(const ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>&, const ShapeFunctionBase<RHSDerived, TestFES, TestSpace>&)
    -> FaceIntegral<Dot<ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>, ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>;

  /**
   * @ingroup FaceIntegralSpecializations
   * @brief Integration of a test operator.
   *
   * Given an operator defined over a test space @f$ V_h @f$
   * @f[
   *   A : V_h \rightarrow \mathbb{R},
   * @f]
   * this class will represent its integral
   * @f[
   *   \int_{\mathcal{F}_h} A(v) \ d\sigma(x) \ .
   * @f]
   */
  template <class NestedDerived, class FES>
  class FaceIntegral<ShapeFunctionBase<NestedDerived, FES, TestSpace>> final
    : public QuadratureRule<ShapeFunctionBase<NestedDerived, FES, TestSpace>>
  {
    public:
      using IntegrandType = ShapeFunctionBase<NestedDerived, FES, TestSpace>;

      using Parent = QuadratureRule<IntegrandType>;

      template <class LHSDerived, class RHSDerived>
      FaceIntegral(
          const FunctionBase<LHSDerived>& lhs,
          const ShapeFunctionBase<RHSDerived, FES, TestSpace>& rhs)
        : FaceIntegral(Dot(lhs, rhs))
      {}

      FaceIntegral(const IntegrandType& integrand)
        : Parent(integrand)
      {}

      FaceIntegral(const FaceIntegral& other)
        : Parent(other)
      {}

      FaceIntegral(FaceIntegral&& other)
        : Parent(std::move(other))
      {}

      Geometry::Region getRegion() const override
      {
        return Geometry::Region::Faces;
      }

      FaceIntegral* copy() const noexcept override
      {
        return new FaceIntegral(*this);
      }

    private:
      std::unique_ptr<IntegrandType> m_integrand;
  };

  template <class NestedDerived, class FES>
  FaceIntegral(const ShapeFunctionBase<NestedDerived, FES, TestSpace>&)
    -> FaceIntegral<ShapeFunctionBase<NestedDerived, FES, TestSpace>>;

  template <class LHSDerived, class RHSDerived, class FES>
  FaceIntegral(const FunctionBase<LHSDerived>&, const ShapeFunctionBase<RHSDerived, FES, TestSpace>&)
    -> FaceIntegral<ShapeFunctionBase<Dot<FunctionBase<LHSDerived>, ShapeFunctionBase<RHSDerived, FES, TestSpace>>, FES, TestSpace>>;
}

#endif
