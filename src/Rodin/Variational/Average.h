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

namespace Rodin::Variational
{
  /**
   * @ingroup RodinVariational
   * @brief Average operator for computing interface averages.
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
        static thread_local Math::SpatialPoint rc1;
        static thread_local Math::SpatialPoint rc2;

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
        it1->getTransformation().inverse(rc1, pc);
        it2->getTransformation().inverse(rc2, pc);
        const Geometry::Point p1(std::cref(*it1), std::cref(rc1), pc);
        const Geometry::Point p2(std::cref(*it2), std::cref(rc2), pc);
        return 0.5 * (this->object(getOperand().getValue(p1)) + this->object(getOperand().getValue(p2)));
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
   * @brief Deduction guide for Average.
   */
  template <class Derived>
  Average(const FunctionBase<Derived>&) -> Average<FunctionBase<Derived>>;

  // template <class Derived, class FES, ShapeFunctionSpaceType Space>
  // class Average<ShapeFunctionBase<Derived, FES, Space>> final
  //   : public ShapeFunctionBase<Average<ShapeFunctionBase<Derived, FES, Space>>>
  // {
  //   public:
  //     using FESType = FES;
  //     static constexpr ShapeFunctionSpaceType SpaceType = Space;

  //     using OperandType = ShapeFunctionBase<Derived, FESType, SpaceType>;

  //     using RangeType = typename FormLanguage::Traits<OperandType>::RangeType;

  //     using Parent = ShapeFunctionBase<Average<ShapeFunctionBase<Derived, FESType, SpaceType>>>;

  //     constexpr
  //     Average(const OperandType& op)
  //       : Parent(op.getFiniteElementSpace()),
  //         m_operand(op.copy())
  //     {}

  //     constexpr
  //     Average(const Average& other)
  //       : Parent(other),
  //         m_operand(other.m_operand->copy())
  //     {}

  //     constexpr
  //     Average(Average&& other)
  //       : Parent(std::move(other)),
  //         m_operand(std::move(other.m_operand))
  //     {}

  //     inline
  //     constexpr
  //     const OperandType& getOperand() const
  //     {
  //       assert(m_operand);
  //       return *m_operand;
  //     }

  //     inline
  //     constexpr
  //     const auto& getLeaf() const
  //     {
  //       return getOperand().getLeaf();
  //     }

  //     inline
  //     constexpr
  //     size_t getDOFs(const Geometry::Polytope& element) const
  //     {
  //       return getOperand().getDOFs(element);
  //     }

  //     inline
  //     const Geometry::Point& getPoint() const
  //     {
  //       return m_operand->getPoint();
  //     }

  //     inline
  //     Average& setPoint(const Geometry::Point& p)
  //     {
  //       m_operand->setPoint(p);
  //       return *this;
  //     }

  //     inline
  //     auto getBasis(size_t local, const Geometry::Point& p) const
  //     {
  //       assert(p.getPolytope().isFace());
  //       const auto& face = p.getPolytope();
  //       const size_t d = face.getDimension();
  //       const auto& mesh = face.getMesh();
  //       const auto& inc = mesh.getConnectivity().getIncidence({ d, d + 1 }, face.getIndex() );
  //       assert(inc.size() == 2);
  //       const Index idx1 = *inc.begin();
  //       const Index idx2 = *std::next(inc.begin());
  //       const auto it1 = mesh.getPolytope(d + 1, idx1);
  //       const auto it2 = mesh.getPolytope(d + 1, idx2);
  //       const auto& pc = p.getPhysicalCoordinates();
  //       const Math::SpatialVector<Real> rc1 = it1->getTransformation().inverse(pc);
  //       const Math::SpatialVector<Real> rc2 = it2->getTransformation().inverse(pc);
  //       const Geometry::Point p1(std::cref(*it1), std::cref(rc1), pc);
  //       const Geometry::Point p2(std::cref(*it2), std::cref(rc2), pc);
  //       const auto& lhs = this->object(getOperand().getBasis(local, p1));
  //       const auto& rhs = this->object(getOperand().getBasis(local, p2));
  //       return 0.5 * (lhs + rhs);
  //     }

  //     inline
  //     constexpr
  //     const auto& getFiniteElementSpace() const
  //     {
  //       return getOperand().getFiniteElementSpace();
  //     }

  //     inline Average* copy() const noexcept override
  //     {
  //       return new Average(*this);
  //     }

  //   private:
  //     std::unique_ptr<OperandType> m_operand;
  // };

  // template <class Derived, class FESType, ShapeFunctionSpaceType SpaceType>
  // Average(const ShapeFunctionBase<Derived, FESType, SpaceType>&)
  //   -> Average<ShapeFunctionBase<Derived, FESType, SpaceType>>;
}

#endif
