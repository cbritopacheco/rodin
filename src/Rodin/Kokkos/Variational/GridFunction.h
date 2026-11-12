/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_VARIATIONAL_GRIDFUNCTION_H
#define RODIN_KOKKOS_VARIATIONAL_GRIDFUNCTION_H
/**
 * @file
 * @brief Backend coefficient storage selected by Kokkos finite element spaces.
 */
#include "Rodin/Variational/GridFunction.h"
#include "Rodin/Variational/TrialFunction.h"
#include "Rodin/Variational/Problem.h"
#include "Rodin/Kokkos/Math.h"
#include "Integral.h"

namespace Rodin::Kokkos::Variational
{
  /// @brief Coefficient vector matching a Kokkos finite element space.
  template <class FES>
  using Coefficients =
    ::Rodin::Kokkos::Math::Vector<typename ::Rodin::FormLanguage::Traits<FES>::ScalarType,
      Internal::BackendOf<FES>>;

  /// @brief Space of a trial or test function.
  template <class Function>
  using FunctionSpace = typename ::Rodin::FormLanguage::Traits<Function>::FESType;

  /// @brief Trial function carrying backend coefficients.
  template <class U>
  concept BackendTrial = Space<FunctionSpace<U>> && requires(const U& u) {
    requires std::is_same_v<std::decay_t<decltype(u.getSolution().getData())>,
      Coefficients<FunctionSpace<U>>>;
  };

  /// @brief Algebraic system matching a backend trial function.
  template <class U>
  using System = ::Rodin::Math::LinearSystem<
    ::Rodin::Kokkos::Math::SparseMatrix<
      typename ::Rodin::FormLanguage::Traits<FunctionSpace<U>>::ScalarType,
      Internal::BackendOf<FunctionSpace<U>>>,
    Coefficients<FunctionSpace<U>>>;

  /// @brief Borrows coefficients already stored in the mesh memory space.
  template <class Scalar, class Backend>
  struct GridFunctionCoefficients<::Rodin::Kokkos::Math::Vector<Scalar, Backend>>
  {
    /// @param fes Owning space. @param data Backend coefficients. @return Borrowed coefficients, without a transfer.
    template <class FES>
    static auto getView([[maybe_unused]] const FES& fes,
      const ::Rodin::Kokkos::Math::Vector<Scalar, Backend>& data)
    {
      static_assert(std::is_same_v<typename FES::MeshType::MemorySpace,
        typename ::Rodin::Kokkos::Math::Vector<Scalar, Backend>::MemorySpace>);
      return data.getView();
    }
  };
}

namespace Rodin::Variational
{
  /// @param u Trial function. @param v Test function on a Kokkos space.
  template <class Solution, ::Rodin::Kokkos::Variational::Space TrialFES,
    ::Rodin::Kokkos::Variational::Space TestFES>
  BilinearForm(const TrialFunction<Solution, TrialFES>& u, const TestFunction<TestFES>& v)
    -> BilinearForm<Solution, TrialFES, TestFES,
      ::Rodin::Kokkos::Math::SparseMatrix<
        typename FormLanguage::Mult<typename FormLanguage::Traits<TrialFES>::ScalarType,
          typename FormLanguage::Traits<TestFES>::ScalarType>::Type,
        ::Rodin::Kokkos::Variational::Internal::BackendOf<TrialFES>>>;

  /// @param v Test function on a Kokkos space.
  template <::Rodin::Kokkos::Variational::Space FES>
  LinearForm(const TestFunction<FES>& v)
    -> LinearForm<FES, ::Rodin::Kokkos::Variational::Coefficients<FES>>;

  /// @param fes Kokkos space. @return Grid function with backend coefficients.
  template <::Rodin::Kokkos::Variational::Space FES>
  GridFunction(const FES& fes)
    -> GridFunction<FES, ::Rodin::Kokkos::Variational::Coefficients<FES>>;

  /// @param fes Kokkos space. @return Trial function with backend coefficients.
  template <::Rodin::Kokkos::Variational::Space FES>
  TrialFunction(const FES& fes)
    -> TrialFunction<GridFunction<FES, ::Rodin::Kokkos::Variational::Coefficients<FES>>,
      FES>;

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
  Problem(U1& u1, U2& u2, U3& u3, Us&... us)
    -> Problem<::Rodin::Kokkos::Variational::System<U1>, U1, U2, U3, Us...>;
}
#endif
