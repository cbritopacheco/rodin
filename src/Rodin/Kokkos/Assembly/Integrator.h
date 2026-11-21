/*
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_ASSEMBLY_INTEGRATOR_H
#define RODIN_KOKKOS_ASSEMBLY_INTEGRATOR_H

/**
 * @file
 * @brief Host-side dispatch of polymorphic Rodin integrators to Kokkos.
 *
 * Assembly receives type-erased core integrator lists. Runtime casts unwrap
 * scaling and negation and select the backend interface once while preparing
 * each assembly input. They do not execute in element or quadrature kernels.
 * Unsupported integrators raise instead of silently selecting host evaluation.
 */

#include <utility>

#include "Rodin/Alert/Exception.h"
#include "Rodin/Variational/BilinearFormIntegrator.h"
#include "Rodin/Variational/LinearFormIntegrator.h"
#include "Rodin/Variational/Mult.h"
#include "Rodin/Variational/UnaryMinus.h"
#include "Rodin/Kokkos/Variational/BilinearFormIntegrator.h"
#include "Rodin/Kokkos/Variational/LinearFormIntegrator.h"

namespace Rodin::Assembly::Internal
{
  /**
   * @brief Finds the device implementation of a local integrator.
   * @tparam Backend Kokkos backend of the mesh context.
   * @param lbfi Local integrator, possibly negated or scaled by a real.
   * @param sign Factor of the contributions.
   * @return Device integrator and the factor of its contributions.
   */
  template <class Backend>
  std::pair<::Rodin::Kokkos::Variational::BilinearFormIntegrator<Backend>*, Real>
  getIntegrator(
    Variational::LocalBilinearFormIntegratorBase<Real>& lbfi, Real sign = Real(1))
  {
    using Negated =
      Variational::UnaryMinus<Variational::LocalBilinearFormIntegratorBase<Real>>;
    using Scaled =
      Variational::Mult<Real, Variational::LocalBilinearFormIntegratorBase<Real>>;
    if (auto* negated = dynamic_cast<Negated*>(&lbfi))
      return getIntegrator<Backend>(negated->getOperand(), -sign);
    if (auto* scaled = dynamic_cast<Scaled*>(&lbfi))
      return getIntegrator<Backend>(scaled->getRHS(), scaled->getLHS() * sign);
    auto* integrator =
      dynamic_cast<::Rodin::Kokkos::Variational::BilinearFormIntegrator<Backend>*>(&lbfi);
    if (!integrator)
    {
      Alert::Exception() << "Local integrator has no device implementation."
                         << Alert::Raise;
    }
    return {integrator, sign};
  }

  /**
   * @brief Finds the device implementation of a global integrator.
   * @tparam Backend Kokkos backend of the mesh context.
   * @param gbfi Global integrator, possibly negated.
   * @param sign Factor of the contributions.
   * @return Device integrator and the factor of its contributions.
   */
  template <class Backend>
  std::pair<::Rodin::Kokkos::Variational::BilinearFormIntegrator<Backend>*, Real>
  getIntegrator(
    Variational::GlobalBilinearFormIntegratorBase<Real>& gbfi, Real sign = Real(1))
  {
    using Negated =
      Variational::UnaryMinus<Variational::GlobalBilinearFormIntegratorBase<Real>>;
    if (auto* negated = dynamic_cast<Negated*>(&gbfi))
      return getIntegrator<Backend>(negated->getOperand(), -sign);
    auto* integrator =
      dynamic_cast<::Rodin::Kokkos::Variational::BilinearFormIntegrator<Backend>*>(&gbfi);
    if (!integrator)
    {
      Alert::Exception() << "Global integrator has no device implementation."
                         << Alert::Raise;
    }
    return {integrator, sign};
  }

  /**
   * @brief Finds the device implementation of a linear integrator.
   * @tparam Backend Kokkos backend of the mesh context.
   * @param lfi Integrator, possibly negated or scaled by a real.
   * @param sign Factor of the contributions.
   * @return Device integrator and the factor of its contributions.
   */
  template <class Backend>
  std::pair<::Rodin::Kokkos::Variational::LinearFormIntegrator<Backend>*, Real>
  getIntegrator(Variational::LinearFormIntegratorBase<Real>& lfi, Real sign = Real(1))
  {
    using Negated = Variational::UnaryMinus<Variational::LinearFormIntegratorBase<Real>>;
    using Scaled = Variational::Mult<Real, Variational::LinearFormIntegratorBase<Real>>;
    if (auto* negated = dynamic_cast<Negated*>(&lfi))
      return getIntegrator<Backend>(negated->getOperand(), -sign);
    if (auto* scaled = dynamic_cast<Scaled*>(&lfi))
      return getIntegrator<Backend>(scaled->getRHS(), scaled->getLHS() * sign);
    auto* integrator =
      dynamic_cast<::Rodin::Kokkos::Variational::LinearFormIntegrator<Backend>*>(&lfi);
    if (!integrator)
      Alert::Exception() << "Integrator has no device implementation." << Alert::Raise;
    return {integrator, sign};
  }
}

#endif
