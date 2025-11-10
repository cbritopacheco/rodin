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
 *
 * @note This file contains placeholder code that is currently commented out.
 *       The Jump operator implementation is under development.
 */
#ifndef RODIN_VARIATIONAL_JUMP_H
#define RODIN_VARIATIONAL_JUMP_H

#include "ForwardDecls.h"
#include "ShapeFunction.h"

namespace Rodin::Variational
{
  /**
   * @ingroup RodinVariational
   * @brief Jump operator for computing interface jumps.
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
   * ## Convention
   * The sign of the jump depends on the chosen normal direction:
   * - For scalar functions: @f$ [\![u]\!] = u^+ n^+ + u^- n^- @f$
   * - For vector functions: @f$ [\![u]\!] = u^+ \otimes n^+ + u^- \otimes n^- @f$
   *
   * ## Usage Example (when implemented)
   * ```cpp
   * // DG interior penalty term
   * auto penalty = InterfaceIntegral(Jump(u), Jump(v));
   * ```
   *
   * @note This class is currently under development. The implementation is
   *       commented out pending completion.
   *
   * @see Average
   */
  // template <class Derived, class FES, ShapeFunctionSpaceType Space>
  // class Jump<ShapeFunctionBase<Derived, FES, Space>> final
  //   : public ShapeFunctionBase<Jump<ShapeFunctionBase<Derived, FES, Space>>>
  // {
  //   public:
  //     using FESType = FES;
  //     static constexpr ShapeFunctionSpaceType SpaceType = Space;

  //     using OperandType = ShapeFunctionBase<Derived, FESType, SpaceType>;

  //     using RangeType = typename FormLanguage::Traits<OperandType>::RangeType;

  //     using Parent = ShapeFunctionBase<Jump<ShapeFunctionBase<Derived, FESType, SpaceType>>>;

  //     constexpr
  //     Jump(const OperandType& op)
  //       : Parent(op.getFiniteElementSpace()),
  //         m_operand(op.copy())
  //     {}

  //     constexpr
  //     Jump(const Jump& other)
  //       : Parent(other),
  //         m_operand(other.m_operand->copy())
  //     {}

  //     constexpr
  //     Jump(Jump&& other)
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
  //       return m_operanderand->getPoint();
  //     }

  //     inline
  //     Transpose& setPoint(const Geometry::Point& p)
  //     {
  //       m_operanderand->setPoint(p);
  //       return *this;
  //     }

  //     inline
  //     auto getBasis(size_t local) const
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
  //       const auto& lhs = this->object(getOperand().getBasis(p1));
  //       const auto& rhs = this->object(getOperand().getBasis(p2));
  //       return lhs - rhs;
  //     }

  //     inline
  //     constexpr
  //     const auto& getFiniteElementSpace() const
  //     {
  //       return getOperand().getFiniteElementSpace();
  //     }

  //     inline Jump* copy() const noexcept override
  //     {
  //       return new Jump(*this);
  //     }

  //   private:
  //     std::unique_ptr<OperandType> m_operand;
  // };

  // template <class Derived, class FESType, ShapeFunctionSpaceType SpaceType>
  // Jump(const ShapeFunctionBase<Derived, FESType, SpaceType>&)
  //   -> Jump<ShapeFunctionBase<Derived, FESType, SpaceType>>;
}

#endif
