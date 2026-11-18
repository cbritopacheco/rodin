/*
 *          Copyright Carlos BRITO PACHECO 2026.
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_VARIATIONAL_INTEGRAL_H
#define RODIN_KOKKOS_VARIATIONAL_INTEGRAL_H

/**
 * @file
 * @brief Kokkos specializations of Rodin integral expressions.
 */

#include "BilinearIntegral.h"
#include "LinearIntegral.h"

namespace Rodin::Kokkos::Variational
{
  using ::Rodin::Kokkos::Geometry::isSelected;
  using ::Rodin::Kokkos::Geometry::isIntegrated;
}

namespace Rodin::Variational
{
  /**
   * @ingroup IntegralSpecializations
   * @brief Device integration of a trial/test dot product over the cells of a
   * Kokkos mesh.
   *
   * Represents
   * @f[
   *   \int_{\mathcal{T}_h} A(u) \cdot B(v) \ dx
   * @f]
   * for device trial and test spaces on a resident Kokkos mesh, as described in
   * Rodin::Kokkos::Variational::BilinearIntegral.
   */
  template <class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
    requires Rodin::Kokkos::Variational::Space<TrialFES> &&
    Rodin::Kokkos::Variational::Space<TestFES>
  class Integral<Dot<ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>,
    ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>
    final : public Rodin::Kokkos::Variational::BilinearIntegral<LHSDerived, RHSDerived,
              TrialFES, TestFES>
  {
    public:
      /// @brief Device integral type.
      using Parent = Rodin::Kokkos::Variational::BilinearIntegral<LHSDerived, RHSDerived,
        TrialFES, TestFES>;

      /// @brief Type of the left operand of the dot product.
      using LHSType = typename Parent::LHSType;

      /// @brief Type of the right operand of the dot product.
      using RHSType = typename Parent::RHSType;

      /// @brief Type of the integrand.
      using IntegrandType = typename Parent::IntegrandType;

      /**
       * @brief Constructs the integral of the dot product of two operators.
       * @param[in] lhs Trial operator @f$ A(u) @f$.
       * @param[in] rhs Test operator @f$ B(v) @f$.
       */
      Integral(const LHSType& lhs, const RHSType& rhs)
        : Integral(Dot(lhs, rhs))
      {}

      /**
       * @brief Constructs the integral of a dot product.
       * @param[in] integrand Dot product of a trial and a test operator.
       */
      Integral(const IntegrandType& integrand)
        : Parent(integrand)
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Integral to copy.
       */
      Integral(const Integral& other) = default;

      /**
       * @brief Move constructor.
       * @param[in] other Integral to move.
       */
      Integral(Integral&& other) = default;

      /**
       * @brief Gets the integration region.
       * @return Region::Cells, the integral running over mesh cells.
       */
      Geometry::Region getRegion() const override
      {
        return Geometry::Region::Cells;
      }

      /**
       * @brief Restricts the integration to mesh attributes, as the parent
       * over().
       * @param[in] attrs Attribute, attributes or set of attributes of the
       * integrated polytopes.
       * @return Reference to this integrator, which keeps its device type.
       */
      template <class... Attributes>
      Integral& over(const Attributes&... attrs)
      {
        Parent::over(attrs...);
        return *this;
      }

      /**
       * @brief Creates a polymorphic copy of this integral.
       * @return Pointer to a new copy.
       */
      Integral* copy() const noexcept override
      {
        return new Integral(*this);
      }
  };

  /**
   * @ingroup IntegralSpecializations
   * @brief Device integration of a test operator over the cells of a Kokkos
   * mesh.
   *
   * Represents
   * @f[
   *   \int_{\mathcal{T}_h} A(v) \ dx
   * @f]
   * for a scalar valued operator of a device test function on a resident Kokkos
   * mesh, as described in Rodin::Kokkos::Variational::LinearIntegral.
   */
  template <class NestedDerived, class FES>
    requires ScalarValuedIntegrand<ShapeFunctionBase<NestedDerived, FES, TestSpace>> &&
    Rodin::Kokkos::Variational::Space<FES>
  class Integral<ShapeFunctionBase<NestedDerived, FES, TestSpace>> final
    : public Rodin::Kokkos::Variational::LinearIntegral<NestedDerived, FES>
  {
    public:
      /// @brief Device integral type.
      using Parent = Rodin::Kokkos::Variational::LinearIntegral<NestedDerived, FES>;

      /// @brief Finite element space of the test functions.
      using FESType = typename Parent::FESType;

      /// @brief Type of the integrand.
      using IntegrandType = typename Parent::IntegrandType;

      /**
       * @brief Constructs the integral of the dot product of a function and a
       * test operator.
       * @param[in] lhs Function @f$ f @f$.
       * @param[in] rhs Test operator @f$ B(v) @f$.
       */
      template <class LHSDerived, class RHSDerived>
      Integral(const FunctionBase<LHSDerived>& lhs,
        const ShapeFunctionBase<RHSDerived, FESType, TestSpace>& rhs)
        : Integral(Dot(lhs, rhs))
      {}

      /**
       * @brief Constructs the integral of a test operator.
       * @param[in] integrand Scalar valued test operator.
       */
      Integral(const IntegrandType& integrand)
        : Parent(integrand)
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Integral to copy.
       */
      Integral(const Integral& other) = default;

      /**
       * @brief Move constructor.
       * @param[in] other Integral to move.
       */
      Integral(Integral&& other) = default;

      /**
       * @brief Gets the integration region.
       * @return Region::Cells, the integral running over mesh cells.
       */
      Geometry::Region getRegion() const override
      {
        return Geometry::Region::Cells;
      }

      /**
       * @brief Restricts the integration to mesh attributes, as the parent
       * over().
       * @param[in] attrs Attribute, attributes or set of attributes of the
       * integrated polytopes.
       * @return Reference to this integrator, which keeps its device type.
       */
      template <class... Attributes>
      Integral& over(const Attributes&... attrs)
      {
        Parent::over(attrs...);
        return *this;
      }

      /**
       * @brief Creates a polymorphic copy of this integral.
       * @return Pointer to a new copy.
       */
      Integral* copy() const noexcept override
      {
        return new Integral(*this);
      }
  };

  /**
   * @ingroup IntegralSpecializations
   * @brief Device integration of a trial/test dot product over the boundary faces of a
   * Kokkos mesh.
   *
   * Represents
   * @f[
   *   \int_{\mathcal{B}_h} A(u) \cdot B(v) \ d\sigma(x)
   * @f]
   * for device trial and test spaces on a resident Kokkos mesh, as described in
   * Rodin::Kokkos::Variational::BilinearIntegral.
   */
  template <class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
    requires Rodin::Kokkos::Variational::Space<TrialFES> &&
    Rodin::Kokkos::Variational::Space<TestFES>
  class BoundaryIntegral<Dot<ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>,
    ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>
    final : public Rodin::Kokkos::Variational::BilinearIntegral<LHSDerived, RHSDerived,
              TrialFES, TestFES>
  {
    public:
      /// @brief Device integral type.
      using Parent = Rodin::Kokkos::Variational::BilinearIntegral<LHSDerived, RHSDerived,
        TrialFES, TestFES>;

      /// @brief Type of the left operand of the dot product.
      using LHSType = typename Parent::LHSType;

      /// @brief Type of the right operand of the dot product.
      using RHSType = typename Parent::RHSType;

      /// @brief Type of the integrand.
      using IntegrandType = typename Parent::IntegrandType;

      /**
       * @brief Constructs the integral of the dot product of two operators.
       * @param[in] lhs Trial operator @f$ A(u) @f$.
       * @param[in] rhs Test operator @f$ B(v) @f$.
       */
      BoundaryIntegral(const LHSType& lhs, const RHSType& rhs)
        : BoundaryIntegral(Dot(lhs, rhs))
      {}

      /**
       * @brief Constructs the integral of a dot product.
       * @param[in] integrand Dot product of a trial and a test operator.
       */
      BoundaryIntegral(const IntegrandType& integrand)
        : Parent(integrand)
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Integral to copy.
       */
      BoundaryIntegral(const BoundaryIntegral& other) = default;

      /**
       * @brief Move constructor.
       * @param[in] other Integral to move.
       */
      BoundaryIntegral(BoundaryIntegral&& other) = default;

      /**
       * @brief Gets the integration region.
       * @return Region::Boundary, the integral running over boundary faces.
       */
      Geometry::Region getRegion() const override
      {
        return Geometry::Region::Boundary;
      }

      /**
       * @brief Restricts the integration to mesh attributes, as the parent
       * over().
       * @param[in] attrs Attribute, attributes or set of attributes of the
       * integrated polytopes.
       * @return Reference to this integrator, which keeps its device type.
       */
      template <class... Attributes>
      BoundaryIntegral& over(const Attributes&... attrs)
      {
        Parent::over(attrs...);
        return *this;
      }

      /**
       * @brief Creates a polymorphic copy of this integral.
       * @return Pointer to a new copy.
       */
      BoundaryIntegral* copy() const noexcept override
      {
        return new BoundaryIntegral(*this);
      }
  };

  /**
   * @ingroup IntegralSpecializations
   * @brief Device integration of a test operator over the boundary faces of a Kokkos
   * mesh.
   *
   * Represents
   * @f[
   *   \int_{\mathcal{B}_h} A(v) \ d\sigma(x)
   * @f]
   * for a scalar valued operator of a device test function on a resident Kokkos
   * mesh, as described in Rodin::Kokkos::Variational::LinearIntegral.
   */
  template <class NestedDerived, class FES>
    requires ScalarValuedIntegrand<ShapeFunctionBase<NestedDerived, FES, TestSpace>> &&
    Rodin::Kokkos::Variational::Space<FES>
  class BoundaryIntegral<ShapeFunctionBase<NestedDerived, FES, TestSpace>> final
    : public Rodin::Kokkos::Variational::LinearIntegral<NestedDerived, FES>
  {
    public:
      /// @brief Device integral type.
      using Parent = Rodin::Kokkos::Variational::LinearIntegral<NestedDerived, FES>;

      /// @brief Finite element space of the test functions.
      using FESType = typename Parent::FESType;

      /// @brief Type of the integrand.
      using IntegrandType = typename Parent::IntegrandType;

      /**
       * @brief Constructs the integral of the dot product of a function and a
       * test operator.
       * @param[in] lhs Function @f$ f @f$.
       * @param[in] rhs Test operator @f$ B(v) @f$.
       */
      template <class LHSDerived, class RHSDerived>
      BoundaryIntegral(const FunctionBase<LHSDerived>& lhs,
        const ShapeFunctionBase<RHSDerived, FESType, TestSpace>& rhs)
        : BoundaryIntegral(Dot(lhs, rhs))
      {}

      /**
       * @brief Constructs the integral of a test operator.
       * @param[in] integrand Scalar valued test operator.
       */
      BoundaryIntegral(const IntegrandType& integrand)
        : Parent(integrand)
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Integral to copy.
       */
      BoundaryIntegral(const BoundaryIntegral& other) = default;

      /**
       * @brief Move constructor.
       * @param[in] other Integral to move.
       */
      BoundaryIntegral(BoundaryIntegral&& other) = default;

      /**
       * @brief Gets the integration region.
       * @return Region::Boundary, the integral running over boundary faces.
       */
      Geometry::Region getRegion() const override
      {
        return Geometry::Region::Boundary;
      }

      /**
       * @brief Restricts the integration to mesh attributes, as the parent
       * over().
       * @param[in] attrs Attribute, attributes or set of attributes of the
       * integrated polytopes.
       * @return Reference to this integrator, which keeps its device type.
       */
      template <class... Attributes>
      BoundaryIntegral& over(const Attributes&... attrs)
      {
        Parent::over(attrs...);
        return *this;
      }

      /**
       * @brief Creates a polymorphic copy of this integral.
       * @return Pointer to a new copy.
       */
      BoundaryIntegral* copy() const noexcept override
      {
        return new BoundaryIntegral(*this);
      }
  };

  /**
   * @ingroup IntegralSpecializations
   * @brief Device integration of a trial/test dot product over every face of a
   * Kokkos mesh.
   *
   * Represents
   * @f[
   *   \int_{\mathcal{F}_h} A(u) \cdot B(v) \ d\sigma(x)
   * @f]
   * for device trial and test spaces on a resident Kokkos mesh, as described in
   * Rodin::Kokkos::Variational::BilinearIntegral.
   */
  template <class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
    requires Rodin::Kokkos::Variational::Space<TrialFES> &&
    Rodin::Kokkos::Variational::Space<TestFES>
  class FaceIntegral<Dot<ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>,
    ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>
    final : public Rodin::Kokkos::Variational::BilinearIntegral<LHSDerived, RHSDerived,
              TrialFES, TestFES>
  {
    public:
      /// @brief Device integral type.
      using Parent = Rodin::Kokkos::Variational::BilinearIntegral<LHSDerived, RHSDerived,
        TrialFES, TestFES>;

      /// @brief Type of the left operand of the dot product.
      using LHSType = typename Parent::LHSType;

      /// @brief Type of the right operand of the dot product.
      using RHSType = typename Parent::RHSType;

      /// @brief Type of the integrand.
      using IntegrandType = typename Parent::IntegrandType;

      /**
       * @brief Constructs the integral of the dot product of two operators.
       * @param[in] lhs Trial operator @f$ A(u) @f$.
       * @param[in] rhs Test operator @f$ B(v) @f$.
       */
      FaceIntegral(const LHSType& lhs, const RHSType& rhs)
        : FaceIntegral(Dot(lhs, rhs))
      {}

      /**
       * @brief Constructs the integral of a dot product.
       * @param[in] integrand Dot product of a trial and a test operator.
       */
      FaceIntegral(const IntegrandType& integrand)
        : Parent(integrand)
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Integral to copy.
       */
      FaceIntegral(const FaceIntegral& other) = default;

      /**
       * @brief Move constructor.
       * @param[in] other Integral to move.
       */
      FaceIntegral(FaceIntegral&& other) = default;

      /**
       * @brief Gets the integration region.
       * @return Region::Faces, the integral running over every face.
       */
      Geometry::Region getRegion() const override
      {
        return Geometry::Region::Faces;
      }

      /**
       * @brief Restricts the integration to mesh attributes, as the parent
       * over().
       * @param[in] attrs Attribute, attributes or set of attributes of the
       * integrated polytopes.
       * @return Reference to this integrator, which keeps its device type.
       */
      template <class... Attributes>
      FaceIntegral& over(const Attributes&... attrs)
      {
        Parent::over(attrs...);
        return *this;
      }

      /**
       * @brief Creates a polymorphic copy of this integral.
       * @return Pointer to a new copy.
       */
      FaceIntegral* copy() const noexcept override
      {
        return new FaceIntegral(*this);
      }
  };

  /**
   * @ingroup IntegralSpecializations
   * @brief Device integration of a test operator over every face of a Kokkos
   * mesh.
   *
   * Represents
   * @f[
   *   \int_{\mathcal{F}_h} A(v) \ d\sigma(x)
   * @f]
   * for a scalar valued operator of a device test function on a resident Kokkos
   * mesh, as described in Rodin::Kokkos::Variational::LinearIntegral.
   */
  template <class NestedDerived, class FES>
    requires ScalarValuedIntegrand<ShapeFunctionBase<NestedDerived, FES, TestSpace>> &&
    Rodin::Kokkos::Variational::Space<FES>
  class FaceIntegral<ShapeFunctionBase<NestedDerived, FES, TestSpace>> final
    : public Rodin::Kokkos::Variational::LinearIntegral<NestedDerived, FES>
  {
    public:
      /// @brief Device integral type.
      using Parent = Rodin::Kokkos::Variational::LinearIntegral<NestedDerived, FES>;

      /// @brief Finite element space of the test functions.
      using FESType = typename Parent::FESType;

      /// @brief Type of the integrand.
      using IntegrandType = typename Parent::IntegrandType;

      /**
       * @brief Constructs the integral of the dot product of a function and a
       * test operator.
       * @param[in] lhs Function @f$ f @f$.
       * @param[in] rhs Test operator @f$ B(v) @f$.
       */
      template <class LHSDerived, class RHSDerived>
      FaceIntegral(const FunctionBase<LHSDerived>& lhs,
        const ShapeFunctionBase<RHSDerived, FESType, TestSpace>& rhs)
        : FaceIntegral(Dot(lhs, rhs))
      {}

      /**
       * @brief Constructs the integral of a test operator.
       * @param[in] integrand Scalar valued test operator.
       */
      FaceIntegral(const IntegrandType& integrand)
        : Parent(integrand)
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Integral to copy.
       */
      FaceIntegral(const FaceIntegral& other) = default;

      /**
       * @brief Move constructor.
       * @param[in] other Integral to move.
       */
      FaceIntegral(FaceIntegral&& other) = default;

      /**
       * @brief Gets the integration region.
       * @return Region::Faces, the integral running over every face.
       */
      Geometry::Region getRegion() const override
      {
        return Geometry::Region::Faces;
      }

      /**
       * @brief Restricts the integration to mesh attributes, as the parent
       * over().
       * @param[in] attrs Attribute, attributes or set of attributes of the
       * integrated polytopes.
       * @return Reference to this integrator, which keeps its device type.
       */
      template <class... Attributes>
      FaceIntegral& over(const Attributes&... attrs)
      {
        Parent::over(attrs...);
        return *this;
      }

      /**
       * @brief Creates a polymorphic copy of this integral.
       * @return Pointer to a new copy.
       */
      FaceIntegral* copy() const noexcept override
      {
        return new FaceIntegral(*this);
      }
  };

  /**
   * @ingroup IntegralSpecializations
   * @brief Device integration of a trial/test dot product over the interface faces of a
   * Kokkos mesh.
   *
   * Represents
   * @f[
   *   \int_{\mathcal{I}_h} A(u) \cdot B(v) \ d\sigma(x)
   * @f]
   * for device trial and test spaces on a resident Kokkos mesh, as described in
   * Rodin::Kokkos::Variational::BilinearIntegral.
   */
  template <class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
    requires Rodin::Kokkos::Variational::Space<TrialFES> &&
    Rodin::Kokkos::Variational::Space<TestFES>
  class InterfaceIntegral<Dot<ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>,
    ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>
    final : public Rodin::Kokkos::Variational::BilinearIntegral<LHSDerived, RHSDerived,
              TrialFES, TestFES>
  {
    public:
      /// @brief Device integral type.
      using Parent = Rodin::Kokkos::Variational::BilinearIntegral<LHSDerived, RHSDerived,
        TrialFES, TestFES>;

      /// @brief Type of the left operand of the dot product.
      using LHSType = typename Parent::LHSType;

      /// @brief Type of the right operand of the dot product.
      using RHSType = typename Parent::RHSType;

      /// @brief Type of the integrand.
      using IntegrandType = typename Parent::IntegrandType;

      /**
       * @brief Constructs the integral of the dot product of two operators.
       * @param[in] lhs Trial operator @f$ A(u) @f$.
       * @param[in] rhs Test operator @f$ B(v) @f$.
       */
      InterfaceIntegral(const LHSType& lhs, const RHSType& rhs)
        : InterfaceIntegral(Dot(lhs, rhs))
      {}

      /**
       * @brief Constructs the integral of a dot product.
       * @param[in] integrand Dot product of a trial and a test operator.
       */
      InterfaceIntegral(const IntegrandType& integrand)
        : Parent(integrand)
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Integral to copy.
       */
      InterfaceIntegral(const InterfaceIntegral& other) = default;

      /**
       * @brief Move constructor.
       * @param[in] other Integral to move.
       */
      InterfaceIntegral(InterfaceIntegral&& other) = default;

      /**
       * @brief Gets the integration region.
       * @return Region::Interface, the integral running over interface faces.
       */
      Geometry::Region getRegion() const override
      {
        return Geometry::Region::Interface;
      }

      /**
       * @brief Restricts the integration to mesh attributes, as the parent
       * over().
       * @param[in] attrs Attribute, attributes or set of attributes of the
       * integrated polytopes.
       * @return Reference to this integrator, which keeps its device type.
       */
      template <class... Attributes>
      InterfaceIntegral& over(const Attributes&... attrs)
      {
        Parent::over(attrs...);
        return *this;
      }

      /**
       * @brief Creates a polymorphic copy of this integral.
       * @return Pointer to a new copy.
       */
      InterfaceIntegral* copy() const noexcept override
      {
        return new InterfaceIntegral(*this);
      }
  };

  /**
   * @ingroup IntegralSpecializations
   * @brief Device integration of a test operator over the interface faces of a Kokkos
   * mesh.
   *
   * Represents
   * @f[
   *   \int_{\mathcal{I}_h} A(v) \ d\sigma(x)
   * @f]
   * for a scalar valued operator of a device test function on a resident Kokkos
   * mesh, as described in Rodin::Kokkos::Variational::LinearIntegral.
   */
  template <class NestedDerived, class FES>
    requires ScalarValuedIntegrand<ShapeFunctionBase<NestedDerived, FES, TestSpace>> &&
    Rodin::Kokkos::Variational::Space<FES>
  class InterfaceIntegral<ShapeFunctionBase<NestedDerived, FES, TestSpace>> final
    : public Rodin::Kokkos::Variational::LinearIntegral<NestedDerived, FES>
  {
    public:
      /// @brief Device integral type.
      using Parent = Rodin::Kokkos::Variational::LinearIntegral<NestedDerived, FES>;

      /// @brief Finite element space of the test functions.
      using FESType = typename Parent::FESType;

      /// @brief Type of the integrand.
      using IntegrandType = typename Parent::IntegrandType;

      /**
       * @brief Constructs the integral of the dot product of a function and a
       * test operator.
       * @param[in] lhs Function @f$ f @f$.
       * @param[in] rhs Test operator @f$ B(v) @f$.
       */
      template <class LHSDerived, class RHSDerived>
      InterfaceIntegral(const FunctionBase<LHSDerived>& lhs,
        const ShapeFunctionBase<RHSDerived, FESType, TestSpace>& rhs)
        : InterfaceIntegral(Dot(lhs, rhs))
      {}

      /**
       * @brief Constructs the integral of a test operator.
       * @param[in] integrand Scalar valued test operator.
       */
      InterfaceIntegral(const IntegrandType& integrand)
        : Parent(integrand)
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Integral to copy.
       */
      InterfaceIntegral(const InterfaceIntegral& other) = default;

      /**
       * @brief Move constructor.
       * @param[in] other Integral to move.
       */
      InterfaceIntegral(InterfaceIntegral&& other) = default;

      /**
       * @brief Gets the integration region.
       * @return Region::Interface, the integral running over interface faces.
       */
      Geometry::Region getRegion() const override
      {
        return Geometry::Region::Interface;
      }

      /**
       * @brief Restricts the integration to mesh attributes, as the parent
       * over().
       * @param[in] attrs Attribute, attributes or set of attributes of the
       * integrated polytopes.
       * @return Reference to this integrator, which keeps its device type.
       */
      template <class... Attributes>
      InterfaceIntegral& over(const Attributes&... attrs)
      {
        Parent::over(attrs...);
        return *this;
      }

      /**
       * @brief Creates a polymorphic copy of this integral.
       * @return Pointer to a new copy.
       */
      InterfaceIntegral* copy() const noexcept override
      {
        return new InterfaceIntegral(*this);
      }
  };
}

#endif
