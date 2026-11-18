/*
 *          Copyright Carlos BRITO PACHECO 2026.
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_VARIATIONAL_INTEGRAL_TRAITS_H
#define RODIN_KOKKOS_VARIATIONAL_INTEGRAL_TRAITS_H

/**
 * @file
 * @brief Compile-time space and differential-form classification for integration.
 */

#include <type_traits>
#include "Rodin/Variational/H1/ForwardDecls.h"
#include "Rodin/Kokkos/Variational/View.h"

namespace Rodin::Kokkos::Variational::Internal
{
  /// @brief Kokkos backend of a device space.
  template <class Mesh>
  struct MeshBackend;

  /// @brief Backend of a resident mesh.
  template <class Backend>
  struct MeshBackend<
    ::Rodin::Geometry::Mesh<::Rodin::Context::Kokkos<::Rodin::Context::Local, Backend>>>
  {
    /// @brief Kokkos backend of the resident mesh.
    using Type = Backend;
  };

  /// @brief Backend of a distributed resident mesh.
  template <class Backend>
  struct MeshBackend<
    ::Rodin::Geometry::Mesh<::Rodin::Context::Kokkos<::Rodin::Context::MPI, Backend>>>
  {
    /// @brief Kokkos backend of the distributed resident mesh.
    using Type = Backend;
  };

  /// @brief Kokkos backend of the mesh of a device space.
  template <class FES>
  using BackendOf = typename MeshBackend<typename FES::MeshType>::Type;

  /// @brief Whether a space is a P1 space.
  template <class FES>
  struct IsP1Space
  {
    /// @brief Whether the space is a P1 space, false in general.
    static constexpr bool Value = false;
  };

  /// @brief P1 space.
  template <class Range, class Mesh>
  struct IsP1Space<::Rodin::Variational::P1<Range, Mesh>>
  {
    /// @brief Whether the space is a P1 space, true here.
    static constexpr bool Value = true;
  };

  /// @brief Whether a space is an H1 space.
  template <class FES>
  struct IsH1Space
  {
    /// @brief Whether the space is an H1 space, false in general.
    static constexpr bool Value = false;
  };

  /// @brief H1 space.
  template <size_t K, class Range, class Mesh>
  struct IsH1Space<::Rodin::Variational::H1<K, Range, Mesh>>
  {
    /// @brief Whether the space is an H1 space, true here.
    static constexpr bool Value = true;
  };

  /// @brief Derivatives of a trial/test pair handled by a dedicated host handler.
  enum class DerivativeKind
  {
    /// @brief No dedicated handler.
    None,

    /// @brief Gradients or Jacobians on both sides.
    Gradients,

    /// @brief Divergence of the trial function against a scalar test function.
    TrialDivergence,

    /// @brief Scalar trial function against the divergence of the test function.
    TestDivergence,

    /// @brief Jacobian of the trial function times a function against a test function.
    JacobianValue
  };

  /**
   * @brief Derivatives of a trial/test pair, selecting the fallback
   * quadrature order of the host handlers.
   *
   * When both spaces are P1, the host handlers of forms with derivatives,
   * such as @f$ \nabla u \cdot \nabla v @f$, @f$ \nabla u : \nabla v @f$
   * or @f$ (\nabla \cdot u) q @f$, fall back to the quadrature order
   * @f$ k_u + k_v - 2 @f$. When both are H1, the order loses one degree per
   * derivative. The device reproduces both.
   */
  template <class LHSDerived, class RHSDerived>
  struct DerivativeForm
  {
    /// @brief Derivatives of the pair, none in general.
    static constexpr DerivativeKind Kind = DerivativeKind::None;
  };

  /// @brief Gradient trial and test operators.
  template <class LHSDerived, class RHSDerived, class TrialFES, class TestFES>
  struct DerivativeForm<
    ::Rodin::Variational::Grad<::Rodin::Variational::ShapeFunction<LHSDerived, TrialFES,
      ::Rodin::Variational::ShapeFunctionSpaceType::Trial>>,
    ::Rodin::Variational::Grad<::Rodin::Variational::ShapeFunction<RHSDerived, TestFES,
      ::Rodin::Variational::ShapeFunctionSpaceType::Test>>>
  {
    /// @brief Gradients on both sides.
    static constexpr DerivativeKind Kind = DerivativeKind::Gradients;
  };

  /// @brief Scaled gradient trial operator and gradient test operator.
  template <class CoefficientDerived, class LHSDerived, class RHSDerived, class TrialFES,
    class TestFES>
  struct DerivativeForm<
    ::Rodin::Variational::Mult<::Rodin::Variational::FunctionBase<CoefficientDerived>,
      ::Rodin::Variational::ShapeFunctionBase<
        ::Rodin::Variational::Grad<::Rodin::Variational::ShapeFunction<LHSDerived,
          TrialFES, ::Rodin::Variational::ShapeFunctionSpaceType::Trial>>,
        TrialFES, ::Rodin::Variational::ShapeFunctionSpaceType::Trial>>,
    ::Rodin::Variational::Grad<::Rodin::Variational::ShapeFunction<RHSDerived, TestFES,
      ::Rodin::Variational::ShapeFunctionSpaceType::Test>>>
  {
    /// @brief Gradients on both sides.
    static constexpr DerivativeKind Kind = DerivativeKind::Gradients;
  };

  /// @brief Divergence trial operator and scalar test function.
  template <class LHSDerived, class RHSDerived, class TrialFES, class TestFES>
  struct DerivativeForm<
    ::Rodin::Variational::Div<::Rodin::Variational::ShapeFunction<LHSDerived, TrialFES,
      ::Rodin::Variational::ShapeFunctionSpaceType::Trial>>,
    ::Rodin::Variational::ShapeFunction<RHSDerived, TestFES,
      ::Rodin::Variational::ShapeFunctionSpaceType::Test>>
  {
    /**
     * @brief Divergence of the trial function when the test space is scalar,
     * none otherwise.
     */
    static constexpr DerivativeKind Kind =
      std::is_same_v<typename TestFES::RangeType, Real>
      ? DerivativeKind::TrialDivergence
      : DerivativeKind::None;
  };

  /// @brief Scalar trial function and divergence test operator.
  template <class LHSDerived, class RHSDerived, class TrialFES, class TestFES>
  struct DerivativeForm<::Rodin::Variational::ShapeFunction<LHSDerived, TrialFES,
                          ::Rodin::Variational::ShapeFunctionSpaceType::Trial>,
    ::Rodin::Variational::Div<::Rodin::Variational::ShapeFunction<RHSDerived, TestFES,
      ::Rodin::Variational::ShapeFunctionSpaceType::Test>>>
  {
    /**
     * @brief Divergence of the test function when the trial space is scalar,
     * none otherwise.
     */
    static constexpr DerivativeKind Kind =
      std::is_same_v<typename TrialFES::RangeType, Real>
      ? DerivativeKind::TestDivergence
      : DerivativeKind::None;
  };

  /// @brief Jacobian trial and test operators.
  template <class LHSDerived, class RHSDerived, class TrialFES, class TestFES>
  struct DerivativeForm<
    ::Rodin::Variational::Jacobian<::Rodin::Variational::ShapeFunction<LHSDerived,
      TrialFES, ::Rodin::Variational::ShapeFunctionSpaceType::Trial>>,
    ::Rodin::Variational::Jacobian<::Rodin::Variational::ShapeFunction<RHSDerived,
      TestFES, ::Rodin::Variational::ShapeFunctionSpaceType::Test>>>
  {
    /// @brief Jacobians on both sides.
    static constexpr DerivativeKind Kind = DerivativeKind::Gradients;
  };

  /// @brief Scaled Jacobian trial operator and Jacobian test operator.
  template <class CoefficientDerived, class LHSDerived, class RHSDerived, class TrialFES,
    class TestFES>
  struct DerivativeForm<
    ::Rodin::Variational::Mult<::Rodin::Variational::FunctionBase<CoefficientDerived>,
      ::Rodin::Variational::ShapeFunctionBase<
        ::Rodin::Variational::Jacobian<::Rodin::Variational::ShapeFunction<LHSDerived,
          TrialFES, ::Rodin::Variational::ShapeFunctionSpaceType::Trial>>,
        TrialFES, ::Rodin::Variational::ShapeFunctionSpaceType::Trial>>,
    ::Rodin::Variational::Jacobian<::Rodin::Variational::ShapeFunction<RHSDerived,
      TestFES, ::Rodin::Variational::ShapeFunctionSpaceType::Test>>>
  {
    /// @brief Jacobians on both sides.
    static constexpr DerivativeKind Kind = DerivativeKind::Gradients;
  };

  /// @brief Jacobian trial operator times a function, against a test function.
  template <class CoefficientDerived, class LHSDerived, class RHSDerived, class TrialFES,
    class TestFES>
  struct DerivativeForm<
    ::Rodin::Variational::Mult<
      ::Rodin::Variational::ShapeFunctionBase<
        ::Rodin::Variational::Jacobian<::Rodin::Variational::ShapeFunction<LHSDerived,
          TrialFES, ::Rodin::Variational::ShapeFunctionSpaceType::Trial>>,
        TrialFES, ::Rodin::Variational::ShapeFunctionSpaceType::Trial>,
      ::Rodin::Variational::FunctionBase<CoefficientDerived>>,
    ::Rodin::Variational::ShapeFunction<RHSDerived, TestFES,
      ::Rodin::Variational::ShapeFunctionSpaceType::Test>>
  {
    /// @brief Jacobian of the trial function times a function.
    static constexpr DerivativeKind Kind = DerivativeKind::JacobianValue;
  };

}

#endif
