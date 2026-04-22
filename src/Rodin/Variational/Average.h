/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Average.h
 * @brief Average operator for interface jump terms in DG methods.
 *
 * This file defines the Average class, which computes the average value of
 * a function across element interfaces. This operator is essential for
 * Discontinuous Galerkin (DG) formulations and interior penalty methods.
 */
#ifndef RODIN_VARIATIONAL_AVERAGE_H
#define RODIN_VARIATIONAL_AVERAGE_H

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
   * @defgroup AverageSpecializations Average Template Specializations
   * @brief Template specializations of the Average class.
   * @see Average
   */

  /**
   * @ingroup AverageSpecializations
   * @brief Average of a FunctionBase instance across an interface.
   *
   * The Average operator computes the average value of a function across an
   * interface (face) between two adjacent elements:
   * @f[
   *   \{\!\{u\}\!\} = \frac{1}{2}(u^+ + u^-)
   * @f]
   * where @f$ u^+ @f$ and @f$ u^- @f$ denote the values from the two sides
   * of the interface.
   *
   * ## Mathematical Foundation
   * In DG methods, functions are allowed to be discontinuous across element
   * interfaces. The average operator is used to define consistent numerical
   * fluxes and penalization terms:
   * - **Interior Penalty**: @f$ \int_{\Gamma_h} \{\!\{\nabla u\}\!\} \cdot [\![v]\!] \, ds @f$
   * - **Consistency**: @f$ \int_{\Gamma_h} \{\!\{\nabla v\}\!\} \cdot [\![u]\!] \, ds @f$
   *
   * where @f$ \Gamma_h @f$ denotes the set of interior faces.
   *
   * ## Usage Example
   * ```cpp
   * // DG interior penalty term
   * auto penalty = InterfaceIntegral(Average(Grad(u)), Jump(v));
   * ```
   *
   * @tparam FunctionDerived Type of the function being averaged
   *
   * @see Jump
   */
  template <class FunctionDerived>
  class Average<FunctionBase<FunctionDerived>> final
    : public FunctionBase<Average<FunctionBase<FunctionDerived>>>
  {
    public:
      /// @brief Type of the operand function
      using OperandType = FunctionBase<FunctionDerived>;

      /// @brief Parent class type
      using Parent = FunctionBase<Average<FunctionBase<FunctionDerived>>>;

      using Parent::traceOf;

      using Parent::operator();

      /**
       * @brief Constructs an average operator for the given function.
       * @param[in] op Function to average across interfaces
       */
      constexpr
      Average(const OperandType& op)
        : m_operand(op.copy())
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Average operator to copy
       */
      constexpr
      Average(const Average& other)
        : Parent(other),
          m_operand(other.m_operand->copy())
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Average operator to move
       */
      constexpr
      Average(Average&& other)
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
       * @brief Evaluates the average at a point on an interface.
       * @param[in] p Point on an interior face
       * @returns Average value @f$ \frac{1}{2}(u^+ + u^-) @f$
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
        const auto& inc = mesh.getConnectivity().getIncidence({ d, d + 1 }, face.getIndex() );
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
        return 0.5 * (this->object(getOperand().getValue(p1)) + this->object(getOperand().getValue(p2)));
      }

      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope& p) const noexcept
      {
        return this->getOperand().getOrder(p);
      }

      /**
       * @brief Creates a copy of this average operator.
       * @returns Pointer to newly allocated copy
       */
      Average* copy() const noexcept override
      {
        return new Average(*this);
      }

    private:
      std::unique_ptr<OperandType> m_operand;
  };

  /**
   * @brief Deduction guide for Average with FunctionBase.
   */
  template <class Derived>
  Average(const FunctionBase<Derived>&) -> Average<FunctionBase<Derived>>;

  /**
   * @ingroup AverageSpecializations
   * @brief Average of a ShapeFunctionBase instance across an interface.
   *
   * Computes the average of shape (trial or test) functions across element
   * interfaces in finite element formulations. This is used in DG bilinear
   * and linear form assembly.
   *
   * ## Usage Example
   * ```cpp
   * // DG consistency term
   * auto consistency = InterfaceIntegral(Average(Grad(u)), Jump(v));
   * ```
   *
   * @tparam NestedDerived Type of the shape function
   * @tparam FES Finite element space type
   * @tparam Space Trial or test function space
   */
  template <class NestedDerived, class FES, ShapeFunctionSpaceType Space>
  class Average<ShapeFunctionBase<NestedDerived, FES, Space>> final
    : public ShapeFunctionBase<Average<ShapeFunctionBase<NestedDerived, FES, Space>>, FES, Space>
  {
    public:
      using FESType = FES;

      using OperandType = ShapeFunctionBase<NestedDerived, FES, Space>;

      using Parent = ShapeFunctionBase<Average<OperandType>, FES, Space>;

      /**
       * @brief Constructs the average of a shape function.
       * @param[in] op Shape function to average
       */
      constexpr
      Average(const OperandType& op)
        : Parent(op.getFiniteElementSpace()),
          m_operand(op.copy()),
          m_ip(nullptr)
      {}

      constexpr
      Average(const Average& other)
        : Parent(other),
          m_operand(other.m_operand->copy()),
          m_ip(nullptr)
      {}

      constexpr
      Average(Average&& other)
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
       * @returns Reference to this average operator
       */
      Average& setIntegrationPoint(const IntegrationPoint& ip)
      {
        m_ip = &ip;
        return *this;
      }

      /**
       * @brief Gets the averaged basis function for local DOF.
       * @param[in] local Local DOF index
       * @returns Average of the basis function @f$ \frac{1}{2}(\phi^+(p) + \phi^-(p)) @f$
       *
       * Evaluates the operand's basis function at the same physical point
       * mapped to both adjacent elements, and returns the average.
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
        return 0.5 * (val1 + val2);
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

      Average* copy() const noexcept override
      {
        return new Average(*this);
      }

    private:
      std::unique_ptr<OperandType> m_operand;
      const IntegrationPoint* m_ip;
  };

  /**
   * @brief Deduction guide for Average with ShapeFunctionBase.
   */
  template <class NestedDerived, class FES, ShapeFunctionSpaceType Space>
  Average(const ShapeFunctionBase<NestedDerived, FES, Space>&)
    -> Average<ShapeFunctionBase<NestedDerived, FES, Space>>;
}

#endif
