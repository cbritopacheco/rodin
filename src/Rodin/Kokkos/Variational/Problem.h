/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_VARIATIONAL_PROBLEM_H
#define RODIN_KOKKOS_VARIATIONAL_PROBLEM_H

/**
 * @file
 * @brief Backend linear-system deduction for variational problems.
 */

#include "Rodin/Variational/Problem.h"
#include "Rodin/Kokkos/Math/LinearSystem.h"
#include "DOFs.h"
#include "GridFunctionData.h"

namespace Rodin::Kokkos::Variational
{
  /// @brief Space of a trial or test function.
  template <class Function>
  using FunctionSpace = typename ::Rodin::FormLanguage::Traits<Function>::FESType;

  /// @brief Trial function carrying backend coefficients.
  template <class U>
  concept BackendTrial = Space<FunctionSpace<U>> && requires(const U& u) {
    requires std::is_same_v<std::decay_t<decltype(u.getSolution().getData())>,
      Data<FunctionSpace<U>>>;
  };

  /// @brief Algebraic system matching a backend trial function.
  template <class U>
  using System = ::Rodin::Math::LinearSystem<
    ::Rodin::Kokkos::Math::SparseMatrix<
      typename ::Rodin::FormLanguage::Traits<FunctionSpace<U>>::ScalarType,
      typename FunctionSpace<U>::MeshType::ContextType::BackendType>,
    Data<FunctionSpace<U>>>;

}

namespace Rodin::Variational
{
  /// @param u Backend trial function. @param v Test function. @return Backend problem.
  template <class U, class V>
    requires ::Rodin::Kokkos::Variational::BackendTrial<U>
  Problem(U& u, V& v) -> Problem<::Rodin::Kokkos::Variational::System<U>, U, V>;

  /**
   * @param u1 First backend trial function. @param u2 Second function.
   * @param u3 Third function. @param us Remaining functions.
   * @return Backend problem containing the supplied fields.
   */
  template <class U1, class U2, class U3, class... Us>
    requires ::Rodin::Kokkos::Variational::BackendTrial<U1>
  Problem(U1& u1, U2& u2, U3& u3,
    Us&... us) -> Problem<::Rodin::Kokkos::Variational::System<U1>, U1, U2, U3, Us...>;
}

#endif
