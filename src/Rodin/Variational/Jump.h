/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Jump.h
 * @brief Jump operator for interface discontinuities in DG methods.
 *
 * This file defines the Jump class, which computes the jump (difference) in
 * function values across element interfaces. This operator is fundamental for
 * Discontinuous Galerkin (DG) formulations where functions are allowed to be
 * discontinuous at element boundaries.
 */
#ifndef RODIN_VARIATIONAL_JUMP_H
#define RODIN_VARIATIONAL_JUMP_H

#include <algorithm>

#include "Rodin/Types.h"
#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/Point.h"

#include "ForwardDecls.h"
#include "IntegrationPoint.h"
#include "ShapeFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup JumpSpecializations Jump Template Specializations
   * @brief Template specializations of the Jump class.
   * @see Jump
   */

  /**
   * @ingroup JumpSpecializations
   * @brief Jump of a FunctionBase instance across an interface.
   *
   * The Jump operator computes the jump (discontinuity) in a function's value
   * across an interface (face) between two adjacent elements:
   * @f[
   *   [\![u]\!] = u^+ - u^-
   * @f]
   * where @f$ u^+ @f$ and @f$ u^- @f$ denote the values from the two sides
   * of the interface (with a chosen orientation).
   *
   * ## Mathematical Foundation
   * In Discontinuous Galerkin (DG) methods, functions are allowed to be
   * discontinuous across element interfaces. The jump operator captures this
   * discontinuity and is used in:
   * - **Interior Penalty terms**: @f$ \int_{\Gamma_h} [\![u]\!] \cdot [\![v]\!] \, ds @f$
   * - **Consistency terms**: @f$ \int_{\Gamma_h} \{\!\{\nabla u\}\!\} \cdot [\![v]\!] \, ds @f$
   * - **Flux terms**: @f$ \int_{\Gamma_h} [\![u]\!] \cdot \{\!\{\nabla v\}\!\} \, ds @f$
   *
   * where @f$ \Gamma_h @f$ denotes the set of interior faces.
   *
   * ## Usage Example
   * ```cpp
   * // Evaluate jump of a function on an interface
   * auto jump_u = Jump(u);
   * ```
   *
   * @tparam FunctionDerived Type of the function being jumped
   *
   * @see Average
   */
  template <class FunctionDerived>
  class Jump<FunctionBase<FunctionDerived>> final
    : public FunctionBase<Jump<FunctionBase<FunctionDerived>>>
  {
    public:
      /// @brief Type of the operand function
      using OperandType = FunctionBase<FunctionDerived>;

      /// @brief Parent class type
      using Parent = FunctionBase<Jump<FunctionBase<FunctionDerived>>>;

      using Parent::traceOf;

      using Parent::operator();

      /**
       * @brief Constructs a jump operator for the given function.
       * @param[in] op Function to compute the jump of across interfaces
       */
      constexpr
      Jump(const OperandType& op)
        : m_operand(op.copy())
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Jump operator to copy
       */
      constexpr
      Jump(const Jump& other)
        : Parent(other),
          m_operand(other.m_operand->copy())
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Jump operator to move
       */
      constexpr
      Jump(Jump&& other)
        : Parent(std::move(other)),
          m_operand(std::move(other.m_operand))
      {}

      /**
       * @brief Gets the operand function.
       * @returns Const reference to the operand
       */
      constexpr
      const auto& getOperand() const
      {
        assert(m_operand);
        return *m_operand;
      }

      /**
       * @brief Evaluates the jump at a point on an interface.
       * @param[in] p Point on an interior face
       * @returns Jump value @f$ u^+(p) - u^-(p) @f$
       *
       * The point must be on an interior face with exactly two adjacent elements.
       */
      auto getValue(const Geometry::Point& p) const
      {
        static thread_local Math::SpatialPoint s_rc1;
        static thread_local Math::SpatialPoint s_rc2;

        assert(p.getPolytope().isFace());
        const auto& face = p.getPolytope();
        const size_t d = face.getDimension();
        const auto& mesh = face.getMesh();
        const auto& inc = mesh.getConnectivity().getIncidence({ d, d + 1 }, face.getIndex());
        assert(inc.size() == 2);
        const Index idx1 = *inc.begin();
        const Index idx2 = *std::next(inc.begin());
        const auto it1 = mesh.getPolytope(d + 1, idx1);
        const auto it2 = mesh.getPolytope(d + 1, idx2);
        const auto& pc = p.getPhysicalCoordinates();
        it1->getTransformation().inverse(s_rc1, pc);
        it2->getTransformation().inverse(s_rc2, pc);
        const Geometry::Point p1(std::cref(*it1), std::cref(s_rc1), pc);
        const Geometry::Point p2(std::cref(*it2), std::cref(s_rc2), pc);
        return this->object(getOperand().getValue(p1)) - this->object(getOperand().getValue(p2));
      }

      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope& p) const noexcept
      {
        return this->getOperand().getOrder(p);
      }

      /**
       * @brief Creates a copy of this jump operator.
       * @returns Pointer to newly allocated copy
       */
      Jump* copy() const noexcept override
      {
        return new Jump(*this);
      }

    private:
      std::unique_ptr<OperandType> m_operand;
  };

  /**
   * @brief Deduction guide for Jump with FunctionBase.
   */
  template <class Derived>
  Jump(const FunctionBase<Derived>&) -> Jump<FunctionBase<Derived>>;

  /**
   * @ingroup JumpSpecializations
   * @brief Jump of a ShapeFunctionBase instance across an interface.
   *
   * Computes the jump of shape (trial or test) functions across element
   * interfaces in finite element formulations. This is used in DG bilinear
   * and linear form assembly.
   *
   * ## Usage Example
   * ```cpp
   * // DG interior penalty term
   * auto penalty = InterfaceIntegral(Jump(u), Jump(v));
   * ```
   *
   * @tparam NestedDerived Type of the shape function
   * @tparam FES Finite element space type
   * @tparam Space Trial or test function space
   */
  template <class NestedDerived, class FES, ShapeFunctionSpaceType Space>
  class Jump<ShapeFunctionBase<NestedDerived, FES, Space>> final
    : public ShapeFunctionBase<Jump<ShapeFunctionBase<NestedDerived, FES, Space>>, FES, Space>
  {
    public:
      using FESType = FES;

      using OperandType = ShapeFunctionBase<NestedDerived, FES, Space>;

      using Parent = ShapeFunctionBase<Jump<OperandType>, FES, Space>;

      /**
       * @brief Constructs the jump of a shape function.
       * @param[in] op Shape function to compute the jump of
       */
      constexpr
      Jump(const OperandType& op)
        : Parent(op.getFiniteElementSpace()),
          m_operand(op.copy()),
          m_ip(nullptr)
      {}

      constexpr
      Jump(const Jump& other)
        : Parent(other),
          m_operand(other.m_operand->copy()),
          m_ip(nullptr)
      {}

      constexpr
      Jump(Jump&& other)
        : Parent(std::move(other)),
          m_operand(std::move(other.m_operand)),
          m_ip(std::exchange(other.m_ip, nullptr))
      {}

      /**
       * @brief Gets the underlying shape function.
       * @returns Reference to the operand
       */
      constexpr
      const OperandType& getOperand() const
      {
        assert(m_operand);
        return *m_operand;
      }

      /**
       * @brief Gets the leaf (underlying trial/test function).
       * @returns Reference to the leaf function
       */
      constexpr
      const auto& getLeaf() const
      {
        return getOperand().getLeaf();
      }

      /**
       * @brief Gets number of degrees of freedom on a polytope.
       * @param[in] element Mesh polytope
       * @returns Number of DOFs
       */
      constexpr
      size_t getDOFs(const Geometry::Polytope& element) const
      {
        return getOperand().getDOFs(element);
      }

      /**
       * @brief Gets the current integration point.
       * @returns Reference to the integration point
       */
      const IntegrationPoint& getIntegrationPoint() const
      {
        assert(m_ip);
        return *m_ip;
      }

      /**
       * @brief Sets the integration point for evaluation.
       * @param[in] ip Integration point on an interface face
       * @returns Reference to this jump operator
       */
      Jump& setIntegrationPoint(const IntegrationPoint& ip)
      {
        m_ip = &ip;
        return *this;
      }

      /**
       * @brief Gets the jump of the basis function for local DOF.
       * @param[in] local Local DOF index
       * @returns Jump of the basis function @f$ \phi^+(p) - \phi^-(p) @f$
       *
       * Evaluates the operand's basis function at the same physical point
       * mapped to both adjacent elements, and returns the difference.
       */
      auto getBasis(size_t local) const
      {
        static thread_local Math::SpatialPoint s_rc1;
        static thread_local Math::SpatialPoint s_rc2;

        assert(m_ip);
        const auto& p = m_ip->getPoint();
        assert(p.getPolytope().isFace());
        const auto& face = p.getPolytope();
        const size_t d = face.getDimension();
        const auto& mesh = face.getMesh();
        const auto& inc = mesh.getConnectivity().getIncidence({ d, d + 1 }, face.getIndex());
        assert(inc.size() == 2);
        const Index idx1 = *inc.begin();
        const Index idx2 = *std::next(inc.begin());
        const auto it1 = mesh.getPolytope(d + 1, idx1);
        const auto it2 = mesh.getPolytope(d + 1, idx2);
        const auto& pc = p.getPhysicalCoordinates();
        it1->getTransformation().inverse(s_rc1, pc);
        it2->getTransformation().inverse(s_rc2, pc);
        const Geometry::Point p1(std::cref(*it1), std::cref(s_rc1), pc);
        const Geometry::Point p2(std::cref(*it2), std::cref(s_rc2), pc);
        const IntegrationPoint ip1(p1, m_ip->getQuadratureFormula(), m_ip->getIndex());
        const IntegrationPoint ip2(p2, m_ip->getQuadratureFormula(), m_ip->getIndex());
        m_operand->setIntegrationPoint(ip1);
        const auto& val1 = this->object(m_operand->getBasis(local));
        m_operand->setIntegrationPoint(ip2);
        const auto& val2 = this->object(m_operand->getBasis(local));
        return val1 - val2;
      }

      /**
       * @brief Gets the finite element space.
       * @returns Reference to the FE space
       */
      constexpr
      const FES& getFiniteElementSpace() const
      {
        return getOperand().getFiniteElementSpace();
      }

      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope& p) const noexcept
      {
        return this->getOperand().getOrder(p);
      }

      Jump* copy() const noexcept override
      {
        return new Jump(*this);
      }

    private:
      std::unique_ptr<OperandType> m_operand;
      const IntegrationPoint* m_ip;
  };

  /**
   * @brief Deduction guide for Jump with ShapeFunctionBase.
   */
  template <class NestedDerived, class FES, ShapeFunctionSpaceType Space>
  Jump(const ShapeFunctionBase<NestedDerived, FES, Space>&)
    -> Jump<ShapeFunctionBase<NestedDerived, FES, Space>>;
}

#endif
