/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_VARIATIONAL_INTERFACESHAPEFUNCTION_H
#define RODIN_VARIATIONAL_INTERFACESHAPEFUNCTION_H

/**
 * @file InterfaceShapeFunction.h
 * @brief Two-sided finite element traces with consistent local numbering.
 */
#include "Rodin/Utility/OptionalReference.h"
#include "ShapeFunction.h"
#include "DOFMapping.h"
#include "Rodin/Geometry/Mesh.h"

namespace Rodin::Variational
{
  /**
   * @brief Evaluation state shared by the average and jump of a shape function.
   * @tparam Derived Average or jump expression.
   * @tparam Operand Cell-local shape expression.
   *
   * # Architecture
   *
   * A global basis restricted to an interior face has two cell-local
   * representations. For the incident cells @f$K^+,K^-@f$, the local support is
   * @f$I=I(K^+)\cup I(K^-)@f$. An index absent from one cell denotes the zero
   * function on that side. The trace is
   * @f$w^+\phi_i|_{K^+}+w^-\phi_i|_{K^-}@f$.
   *
   * # Implementation
   *
   * Entity binding constructs the union and its side-local numbering once.
   * Point binding maps the face point into each cell and prepares two
   * independent owned operand evaluations. Basis reads neither rebind an
   * operand nor retain references to temporary points. Mapped cell coordinates
   * do not inherit the face quadrature's reference association.
   */
  template <class Derived, class Operand>
  class InterfaceShapeFunctionBase;

  /**
   * @brief Two-sided evaluation of a cell-local shape expression.
   * @tparam Derived Concrete average or jump node.
   * @tparam NestedDerived Cell-local expression node.
   * @tparam FES Finite element space providing the cell DOF maps.
   * @tparam Space Trial or test grading.
   */
  template <class Derived, class NestedDerived, class FES, ShapeFunctionSpaceType Space>
  class InterfaceShapeFunctionBase<Derived, ShapeFunctionBase<NestedDerived, FES, Space>>
    : public ShapeFunctionBase<Derived, FES, Space>
  {
    public:
      /// @brief Cell-local operand type.
      using OperandType = ShapeFunctionBase<NestedDerived, FES, Space>;

      /// @brief Shape expression base.
      using Parent = ShapeFunctionBase<Derived, FES, Space>;

      /**
       * @brief Constructs a weighted two-sided trace.
       * @param operand Cell-local expression to trace.
       * @param lhs First side weight.
       * @param rhs Second side weight.
       */
      InterfaceShapeFunctionBase(const OperandType& operand, Real lhs, Real rhs)
        : Parent(operand.getFiniteElementSpace()),
          m_lhs(operand.copy()),
          m_rhs(operand.copy()),
          m_lhsWeight(lhs),
          m_rhsWeight(rhs)
      {}

      /**
       * @brief Copies the expression without retaining borrowed evaluation state.
       * @param other Expression to copy.
       */
      InterfaceShapeFunctionBase(const InterfaceShapeFunctionBase& other)
        : Parent(other),
          m_lhs(other.m_lhs->copy()),
          m_rhs(other.m_rhs->copy()),
          m_lhsWeight(other.m_lhsWeight),
          m_rhsWeight(other.m_rhsWeight)
      {}

      /**
       * @brief Moves the expression without retaining borrowed evaluation state.
       * @param other Expression to move.
       */
      InterfaceShapeFunctionBase(InterfaceShapeFunctionBase&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)),
          m_rhs(std::move(other.m_rhs)),
          m_lhsWeight(other.m_lhsWeight),
          m_rhsWeight(other.m_rhsWeight)
      {}

      /// @returns The underlying cell-local expression.
      const OperandType& getOperand() const
      {
        return *m_lhs;
      }

      /// @returns The underlying trial or test function.
      const auto& getLeaf() const
      {
        return getOperand().getLeaf();
      }

      /// @returns Borrowed union of the bound incident-cell DOFs.
      Utility::OptionalReference<const IndexArray> getDOFIndices() const
      {
        assert(m_lhsCell && m_rhsCell);
        return m_mapping.getIndices();
      }

      /**
       * @brief Counts the union of the incident-cell DOFs.
       * @param face Interior face whose support is requested.
       * @returns Number of distinct global DOFs in its two cells.
       */
      size_t getDOFs(const Geometry::Polytope& face) const
      {
        auto working = static_cast<const Derived&>(*this);
        working.setPolytope(face);
        return working.getDOFIndices()->size();
      }

      /**
       * @brief Prepares independent cell-local operands and their union numbering.
       * @param face Interior face of the finite element mesh.
       * @returns This expression after binding.
       */
      Derived& setPolytope(const Geometry::Polytope& face)
      {
        m_ip = nullptr;
        m_faceIndex.reset();
        const auto& mesh = getMesh();
        if (&face.getMesh() != &mesh || !face.isFace())
        {
          Alert::MemberFunctionException(*this, __func__)
            << "A two-sided trace requires a face of its finite element mesh."
            << Alert::Raise;
        }
        const auto dimension = face.getDimension();
        const auto& incident = mesh.getConnectivity().getIncidence(
          {dimension, dimension + 1}, face.getIndex());
        if (incident.size() != 2)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "A two-sided trace requires exactly two incident cells." << Alert::Raise;
        }
        const auto first = incident.begin();
        m_lhsCell.emplace(*mesh.getPolytope(dimension + 1, *first));
        m_rhsCell.emplace(*mesh.getPolytope(dimension + 1, *std::next(first)));
        m_lhs->setPolytope(*m_lhsCell);
        m_rhs->setPolytope(*m_rhsCell);
        const auto& fes = this->getFiniteElementSpace();
        m_mapping.setLHS(fes.getDOFs(dimension + 1, m_lhsCell->getIndex()));
        m_mapping.setRHS(fes.getDOFs(dimension + 1, m_rhsCell->getIndex()));
        m_faceIndex = face.getIndex();
        return static_cast<Derived&>(*this);
      }

      /**
       * @brief Maps one face point into both incident-cell charts.
       * @param ip Face integration point retained until the next binding.
       * @returns This expression after preparing both traces.
       */
      Derived& setIntegrationPoint(const IntegrationPoint& ip)
      {
        const auto& point = ip.getPoint();
        const auto& face = point.getPolytope();
        if (!m_faceIndex || &face.getMesh() != &getMesh() || !face.isFace() ||
          face.getIndex() != *m_faceIndex)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "The integration point must belong to the bound interior face."
            << Alert::Raise;
        }
        const auto& physical = point.getPhysicalCoordinates();
        Math::SpatialPoint lhs, rhs;
        m_lhsCell->getTransformation().inverse(lhs, physical);
        m_rhsCell->getTransformation().inverse(rhs, physical);
        m_lhsPoint.emplace(*m_lhsCell, lhs, physical);
        m_rhsPoint.emplace(*m_rhsCell, rhs, physical);
        m_lhsIP.emplace(*m_lhsPoint);
        m_rhsIP.emplace(*m_rhsPoint);
        m_lhs->setIntegrationPoint(*m_lhsIP);
        m_rhs->setIntegrationPoint(*m_rhsIP);
        m_ip = &ip;
        return static_cast<Derived&>(*this);
      }

      /// @returns The bound face integration point.
      const IntegrationPoint& getIntegrationPoint() const
      {
        assert(m_ip);
        return *m_ip;
      }

      /**
       * @brief Evaluates a union basis using its independent side-local indices.
       * @param local Position in the union numbering.
       * @returns Owning weighted two-sided basis value.
       */
      auto getBasis(size_t local) const
      {
        assert(m_ip);
        using Value = typename FormLanguage::Traits<OperandType>::RangeType;
        const auto& lhs = m_mapping.getLHS(local);
        const auto& rhs = m_mapping.getRHS(local);
        Value value = lhs ? Value(m_lhs->getBasis(*lhs)) : Value(m_rhs->getBasis(*rhs));
        value *= lhs ? m_lhsWeight : m_rhsWeight;
        if (lhs && rhs)
          value += Value(m_rhs->getBasis(*rhs)) * m_rhsWeight;
        return value;
      }

      /**
       * @brief Returns the underlying cell expression's polynomial order.
       * @param face Face on which the trace is integrated.
       * @returns Trace polynomial order, when known.
       */
      Optional<size_t> getOrder(const Geometry::Polytope& face) const noexcept
      {
        return getOperand().getOrder(face);
      }

    private:
      /// @returns The finite-element mesh whose iterators bind both cell traces.
      const auto& getMesh() const
      {
        return this->getFiniteElementSpace().getMesh();
      }

      std::unique_ptr<OperandType> m_lhs, m_rhs;
      Real m_lhsWeight, m_rhsWeight;
      DOFMapping<> m_mapping;
      Optional<Geometry::Polytope> m_lhsCell, m_rhsCell;
      Optional<Geometry::Point> m_lhsPoint, m_rhsPoint;
      Optional<IntegrationPoint> m_lhsIP, m_rhsIP;
      Optional<Index> m_faceIndex;
      const IntegrationPoint* m_ip = nullptr;
  };
}
#endif
