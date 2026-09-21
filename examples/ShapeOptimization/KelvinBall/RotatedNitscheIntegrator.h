/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef KELVIN_BALL_ROTATED_NITSCHE_INTEGRATOR_H
#define KELVIN_BALL_ROTATED_NITSCHE_INTEGRATOR_H

#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <Rodin/Geometry.h>
#include <Rodin/Location.h>
#include <Rodin/Advection/Lagrangian.h>
#include <Rodin/Alert/Info.h>
#include <Rodin/Alert/Notation.h>
#include <Rodin/Variational.h>

#include "Common.h"

namespace KelvinBall
{
  using namespace Rodin;
  using namespace Rodin::Geometry;
  using namespace Rodin::Variational;

  template <class Mesh>
  class AttributeFaceLocator
  {
    public:
      AttributeFaceLocator(const Mesh& mesh, const FlatSet<Attribute>& attributes,
        Real physicalTolerance = 1e-10, Real referenceTolerance = 1e-10)
        : m_mesh(mesh)
      {
        const size_t faceDimension = mesh.getDimension() - 1;
        for (const Attribute attribute : attributes)
        {
          typename SubMesh<Context::Local>::Builder builder;
          builder.initialize(mesh);
          for (auto face = mesh.getPolytope(faceDimension); face; ++face)
          {
            if (face->getAttribute() == attribute)
              builder.include(faceDimension, face->getIndex());
          }
          Surface surface;
          surface.mesh = std::make_unique<SubMesh<Context::Local>>(builder.finalize());
          surface.strict =
            std::make_unique<Location::AABB<SubMesh<Context::Local>>>(*surface.mesh);
          surface.relaxed =
            std::make_unique<Location::AABB<SubMesh<Context::Local>>>(*surface.mesh);
          surface.relaxed->setTolerance(physicalTolerance)
            .setReferenceTolerance(referenceTolerance);
          m_surfaces.emplace(attribute, std::move(surface));
        }
      }

      Optional<Geometry::Point> locate(
        Attribute attribute, const Math::SpatialPoint& point) const
      {
        const auto position = m_surfaces.find(attribute);
        if (position == m_surfaces.end())
          return {};
        const auto& surface = *position->second.mesh;
        const auto toParent = [&](const Geometry::Point& mapped) {
          const Index parent = surface.getPolytopeMap(surface.getDimension())
                                 .left[mapped.getPolytope().getIndex()];
          const auto polytope = m_mesh.get().getPolytope(surface.getDimension(), parent);
          if (!polytope)
            throw std::runtime_error("A located face has no parent polytope.");
          return Geometry::Point(
            *polytope, mapped.getReferenceCoordinates(), mapped.getPhysicalCoordinates());
        };
        if (const auto mapped =
              position->second.strict->locate(surface.getDimension(), point))
          return toParent(*mapped);
        if (const auto mapped =
              position->second.relaxed->locate(surface.getDimension(), point))
        {
          ++m_relaxedHits;
          return toParent(*mapped);
        }
        return {};
      }

      size_t getRelaxedHitCount() const
      {
        return m_relaxedHits;
      }

    private:
      struct Surface
      {
          std::unique_ptr<SubMesh<Context::Local>> mesh;
          std::unique_ptr<Location::AABB<SubMesh<Context::Local>>> strict;
          std::unique_ptr<Location::AABB<SubMesh<Context::Local>>> relaxed;
      };

      std::reference_wrapper<const Mesh> m_mesh;
      std::map<Attribute, Surface> m_surfaces;
      mutable size_t m_relaxedHits = 0;
  };

  /**
   * @brief Mapped-face assembly of the rotational chamber conditions.
   *
   * The integrator owns the AABB map of the master cuts. Each assembly method
   * traverses slave quadrature points, locates the rotated master trace, and
   * inserts the resulting nonlocal coupling into the supplied linear system.
   */
  class RotatedNitscheIntegrator
  {
    public:
      using Locator = AttributeFaceLocator<Mesh>;

      RotatedNitscheIntegrator(const Mesh& mesh,
        const FlatSet<Attribute>& masterCuts = MasterCuts, Real physicalTolerance = 1e-10,
        Real referenceTolerance = 1e-10);

      const Locator& getLocator() const;

      template <size_t VelocityOrder, class VelocitySpace, class PressureSpace,
        class Offsets, class LinearSystem>
      void assembleStokes(const VelocitySpace& velocity, const PressureSpace& pressure,
        const Offsets& offsets, LinearSystem& system, Real viscosity, Real penalty,
        Real pressureDiffusion, const FlatSet<Attribute>& fixedVelocityBoundaries) const;

      template <class VectorSpace, class LinearSystem>
      void assembleVector(const VectorSpace& space, LinearSystem& system, Real diffusion,
        Real penalty, const FlatSet<Attribute>& fixedBoundaries) const;

      template <class ScalarSpace, class LinearSystem>
      void assembleScalarTracePenalty(const ScalarSpace& space, LinearSystem& system,
        Real penalty, const FlatSet<Attribute>& fixedBoundaries = {}) const;

      template <class ScalarSpace, class LinearSystem, class Reference>
      void assembleScalarTracePenalty(const ScalarSpace& space, LinearSystem& system,
        Real penalty, const Reference& reference) const;

      template <class U>
      Real vectorJump(const U& field) const;

      template <class U>
      Real scalarJump(const U& field) const;

      template <class U0, class U1, class U2>
      Real familyJump(const U0& u0, const U1& u1, const U2& u2) const;

    private:
      class Internal;

      static Alert::Text<Alert::YellowT> diagnosticHeading(const std::string& text)
      {
        Alert::Text<Alert::YellowT> heading(Alert::Yellow, text);
        return heading.setBold();
      }

      static std::string diagnosticLabel(const std::string& text)
      {
        constexpr size_t width = 36;
        return text + std::string(text.size() < width ? width - text.size() : 1, ' ');
      }

      Locator m_locator;
  };

  class RotatedNitscheIntegrator::Internal
  {
    public:
      struct ScalarCellBasis
      {
          IndexVector dofs;
          std::vector<Real> values;
          std::vector<Math::SpatialVector<Real>> gradients;
      };

      struct StateBasis
      {
          Index dof;
          Math::SpatialMatrix<Real> jump;
          Math::SpatialMatrix<Real> traction;
          Math::SpatialVector<Real> pressureJump;
          Math::SpatialVector<Real> pressureFlux;
      };

      struct VectorBasis
      {
          Index dof;
          Math::SpatialVector<Real> jump;
          Math::SpatialVector<Real> flux;
      };

      static Math::SpatialMatrix<Real> zeroMatrix()
      {
        Math::SpatialMatrix<Real> result(3, 3);
        result.setZero();
        return result;
      }

      static Math::SpatialVector<Real> zeroVector()
      {
        Math::SpatialVector<Real> result(3);
        result.setZero();
        return result;
      }

      static void setColumn(Math::SpatialMatrix<Real>& matrix, size_t column,
        const Math::SpatialVector<Real>& value)
      {
        for (size_t row = 0; row < 3; ++row)
          matrix(row, column) = value(row);
      }

      static Index getIncidentCell(const Polytope& face)
      {
        const auto& mesh = face.getMesh();
        const size_t faceDimension = mesh.getDimension() - 1;
        const auto& incidence = mesh.getConnectivity().getIncidence(
          {faceDimension, mesh.getDimension()}, face.getIndex());
        if (incidence.size() != 1)
          throw std::runtime_error(
            "A chamber cut face must have one incident fluid cell.");
        return incidence[0];
      }

      template <size_t Order, class Space>
      static ScalarCellBasis evaluateScalarBasis(
        const Space& space, const Polytope& face, const Math::SpatialPoint& physical)
      {
        const Index cellIndex = getIncidentCell(face);
        const auto cellIterator =
          face.getMesh().getPolytope(face.getMesh().getDimension(), cellIndex);
        if (!cellIterator)
          throw std::runtime_error("A chamber cut face has no incident fluid cell.");
        const auto& cell = *cellIterator;
        Math::SpatialPoint reference;
        cell.getTransformation().inverse(reference, physical);
        const Geometry::Point point(cell, std::cref(reference), physical);
        const H1Element<Order, Real> element(cell.getGeometry());
        ScalarCellBasis result;
        const auto dofs = space.getDOFs(cell.getDimension(), cell.getIndex());
        result.dofs.assign(dofs.begin(), dofs.end());
        result.values.resize(element.getCount());
        result.gradients.resize(element.getCount());
        for (size_t local = 0; local < element.getCount(); ++local)
        {
          const auto& basis = element.getBasis(local);
          result.values[local] = basis(reference);
          result.gradients[local] =
            point.getJacobianInverse().transpose() * basis.getGradient()(reference);
        }
        return result;
      }

      static Real diameter(const Polytope& face)
      {
        const auto& mesh = face.getMesh();
        Real result = 0;
        const auto& vertices = face.getVertices();
        for (size_t i = 0; i < vertices.size(); ++i)
        {
          for (size_t j = i + 1; j < vertices.size(); ++j)
          {
            result = std::max(result,
              (mesh.getVertexCoordinates(vertices[i]) -
                mesh.getVertexCoordinates(vertices[j]))
                .norm());
          }
        }
        if (!(result > 0))
          throw std::runtime_error("A chamber cut face has zero diameter.");
        return result;
      }

      static Real frobenius(
        const Math::SpatialMatrix<Real>& lhs, const Math::SpatialMatrix<Real>& rhs)
      {
        return lhs.dot(rhs);
      }

      template <class LinearSystem>
      static void addWithEliminatedColumns(LinearSystem& system,
        const std::vector<Eigen::Triplet<Real>>& entries, const IndexSet& fixed)
      {
        auto& matrix = system.getOperator();
        auto& rhs = system.getVector();
        Math::SparseMatrix<Real> addition(matrix.rows(), matrix.cols());
        addition.setFromTriplets(entries.begin(), entries.end());
        std::vector<Eigen::Triplet<Real>> freeEntries;
        freeEntries.reserve(addition.nonZeros());
        for (Eigen::Index column = 0; column < addition.outerSize(); ++column)
        {
          for (Math::SparseMatrix<Real>::InnerIterator coefficient(addition, column);
               coefficient; ++coefficient)
          {
            if (fixed.contains(coefficient.row()))
              continue;
            if (fixed.contains(coefficient.col()))
              rhs(coefficient.row()) -= coefficient.value() * rhs(coefficient.col());
            else
              freeEntries.emplace_back(
                coefficient.row(), coefficient.col(), coefficient.value());
          }
        }
        Math::SparseMatrix<Real> freeAddition(matrix.rows(), matrix.cols());
        freeAddition.setFromTriplets(freeEntries.begin(), freeEntries.end());
        matrix += freeAddition;
      }

      template <class Space>
      static IndexSet getFixedDOFs(const Space& space,
        const FlatSet<Attribute>& attributes, const std::array<size_t, 3>& blocks,
        const std::vector<size_t>& offsets)
      {
        const auto& mesh = space.getMesh();
        const size_t faceDimension = mesh.getDimension() - 1;
        IndexSet fixed;
        for (auto face = mesh.getPolytope(faceDimension); face; ++face)
        {
          const auto attribute = face->getAttribute();
          if (!attribute || !attributes.contains(*attribute))
            continue;
          const auto dofs = space.getDOFs(faceDimension, face->getIndex());
          for (const size_t block : blocks)
          {
            for (const Index dof : dofs)
              fixed.insert(offsets[block] + dof);
          }
        }
        return fixed;
      }

      template <class Space>
      static IndexSet getFixedDOFs(
        const Space& space, const FlatSet<Attribute>& attributes)
      {
        const auto& mesh = space.getMesh();
        const size_t faceDimension = mesh.getDimension() - 1;
        IndexSet fixed;
        for (auto face = mesh.getPolytope(faceDimension); face; ++face)
        {
          const auto attribute = face->getAttribute();
          if (!attribute || !attributes.contains(*attribute))
            continue;
          const auto dofs = space.getDOFs(faceDimension, face->getIndex());
          fixed.insert(dofs.begin(), dofs.end());
        }
        return fixed;
      }
  };

  template <size_t VelocityOrder, class VelocitySpace, class PressureSpace, class Offsets,
    class LinearSystem>
  void RotatedNitscheIntegrator::assembleStokes(const VelocitySpace& velocity,
    const PressureSpace& pressure, const Offsets& offsets, LinearSystem& system,
    Real viscosity, Real penalty, Real pressureDiffusion,
    const FlatSet<Attribute>& fixedVelocityBoundaries) const
  {
    const auto& mesh = velocity.getMesh();
    const size_t faceDimension = mesh.getDimension() - 1;
    const std::array<size_t, 3> velocityBlocks{0, 2, 4};
    const std::array<size_t, 3> pressureBlocks{1, 3, 5};
    const std::vector<size_t> blockOffsets(offsets.begin(), offsets.end());
    const IndexSet fixed = Internal::getFixedDOFs(
      velocity, fixedVelocityBoundaries, velocityBlocks, blockOffsets);
    FaceNormal normal(mesh);
    size_t quadraturePoints = 0;
    Real normalResidual = 0;

    for (const RotationPair& pair : RotationPairs)
    {
      std::vector<Eigen::Triplet<Real>> entries;
      for (auto face = mesh.getPolytope(faceDimension); face; ++face)
      {
        if (face->getAttribute() != pair.slave)
          continue;
        const auto& qf = QF::PolytopeQuadratureFormula::get(4, face->getGeometry());
        const auto& quadrature = face->getQuadrature(qf);
        for (size_t qp = 0; qp < quadrature.getSize(); ++qp)
        {
          ++quadraturePoints;
          const auto& slavePoint = quadrature.getPoint(qp);
          const auto mapped =
            m_locator.locate(pair.master, pair.rotation * slavePoint.vector());
          if (!mapped)
            throw std::runtime_error(
              "A rotated Nitsche quadrature point was not located.");
          const auto slaveVelocity = Internal::evaluateScalarBasis<VelocityOrder>(
            velocity, *face, slavePoint.getPhysicalCoordinates());
          const auto slavePressure = Internal::evaluateScalarBasis<1>(
            pressure, *face, slavePoint.getPhysicalCoordinates());
          const auto& masterFace = mapped->getPolytope();
          const auto masterVelocity = Internal::evaluateScalarBasis<VelocityOrder>(
            velocity, masterFace, mapped->getPhysicalCoordinates());
          const auto masterPressure = Internal::evaluateScalarBasis<1>(
            pressure, masterFace, mapped->getPhysicalCoordinates());
          const Math::SpatialVector<Real> slaveNormal = normal.getValue(slavePoint);
          const Math::SpatialVector<Real> masterNormal = normal.getValue(*mapped);
          normalResidual =
            std::max(normalResidual, (masterNormal + pair.rotation * slaveNormal).norm());
          const Real weight = qf.getWeight(qp) * slavePoint.getDistortion();
          const Real stabilization = penalty * viscosity / Internal::diameter(*face);
          const Real pressureStabilization =
            penalty * pressureDiffusion / Internal::diameter(*face);

          std::vector<Internal::StateBasis> basis;
          basis.reserve(96);
          for (size_t load = 0; load < 3; ++load)
          {
            for (size_t node = 0; node < masterVelocity.values.size(); ++node)
            {
              for (size_t component = 0; component < 3; ++component)
              {
                Math::SpatialMatrix<Real> gradient(3, 3);
                gradient.setZero();
                for (size_t direction = 0; direction < 3; ++direction)
                  gradient(component, direction) =
                    masterVelocity.gradients[node](direction);
                const Math::SpatialVector<Real> traction =
                  viscosity * (gradient + gradient.transpose()) * masterNormal;
                Internal::StateBasis value{
                  static_cast<Index>(offsets[velocityBlocks[load]] +
                    masterVelocity.dofs[3 * node + component]),
                  Internal::zeroMatrix(), Internal::zeroMatrix(), Internal::zeroVector(),
                  Internal::zeroVector()};
                value.jump(component, load) = masterVelocity.values[node];
                Internal::setColumn(value.traction, load, 0.5 * traction);
                basis.emplace_back(std::move(value));
              }
            }
            for (size_t node = 0; node < masterPressure.values.size(); ++node)
            {
              Internal::StateBasis value{
                static_cast<Index>(
                  offsets[pressureBlocks[load]] + masterPressure.dofs[node]),
                Internal::zeroMatrix(), Internal::zeroMatrix(), Internal::zeroVector(),
                Internal::zeroVector()};
              Internal::setColumn(
                value.traction, load, -0.5 * masterPressure.values[node] * masterNormal);
              value.pressureJump(load) = masterPressure.values[node];
              value.pressureFlux(load) = 0.5 * pressureDiffusion *
                masterPressure.gradients[node].dot(masterNormal);
              basis.emplace_back(std::move(value));
            }
          }
          for (size_t source = 0; source < 3; ++source)
          {
            for (size_t node = 0; node < slaveVelocity.values.size(); ++node)
            {
              for (size_t component = 0; component < 3; ++component)
              {
                Math::SpatialMatrix<Real> gradient(3, 3);
                gradient.setZero();
                for (size_t direction = 0; direction < 3; ++direction)
                  gradient(component, direction) =
                    slaveVelocity.gradients[node](direction);
                const Math::SpatialVector<Real> transformedValue =
                  pair.rotation.col(component) * slaveVelocity.values[node];
                const Math::SpatialVector<Real> transformedTraction = pair.rotation *
                  (viscosity * (gradient + gradient.transpose()) * slaveNormal);
                Internal::StateBasis value{
                  static_cast<Index>(offsets[velocityBlocks[source]] +
                    slaveVelocity.dofs[3 * node + component]),
                  Internal::zeroMatrix(), Internal::zeroMatrix(), Internal::zeroVector(),
                  Internal::zeroVector()};
                for (size_t load = 0; load < 3; ++load)
                {
                  const Real coefficient = pair.rotation(load, source);
                  Internal::setColumn(value.jump, load, -coefficient * transformedValue);
                  Internal::setColumn(
                    value.traction, load, -0.5 * coefficient * transformedTraction);
                }
                basis.emplace_back(std::move(value));
              }
            }
            for (size_t node = 0; node < slavePressure.values.size(); ++node)
            {
              const Math::SpatialVector<Real> transformedTraction =
                pair.rotation * (-slavePressure.values[node] * slaveNormal);
              Internal::StateBasis value{
                static_cast<Index>(
                  offsets[pressureBlocks[source]] + slavePressure.dofs[node]),
                Internal::zeroMatrix(), Internal::zeroMatrix(), Internal::zeroVector(),
                Internal::zeroVector()};
              for (size_t load = 0; load < 3; ++load)
              {
                const Real coefficient = pair.rotation(load, source);
                Internal::setColumn(
                  value.traction, load, -0.5 * coefficient * transformedTraction);
                value.pressureJump(load) = -coefficient * slavePressure.values[node];
                value.pressureFlux(load) = -0.5 * coefficient * pressureDiffusion *
                  slavePressure.gradients[node].dot(slaveNormal);
              }
              basis.emplace_back(std::move(value));
            }
          }

          for (const auto& test : basis)
          {
            for (const auto& trial : basis)
            {
              const Real coefficient = weight *
                (-Internal::frobenius(trial.traction, test.jump) -
                  Internal::frobenius(test.traction, trial.jump) +
                  stabilization * Internal::frobenius(trial.jump, test.jump) -
                  trial.pressureFlux.dot(test.pressureJump) -
                  test.pressureFlux.dot(trial.pressureJump) +
                  pressureStabilization * trial.pressureJump.dot(test.pressureJump));
              if (std::abs(coefficient) > 1e-14)
                entries.emplace_back(test.dof, trial.dof, coefficient);
            }
          }
        }
      }
      Internal::addWithEliminatedColumns(system, entries, fixed);
    }
    Alert::Info() << diagnosticHeading("Rotated Nitsche assembly") << Alert::NewLine
                  << diagnosticLabel("Quadrature points:")
                  << Alert::Notation::Number(quadraturePoints) << Alert::NewLine
                  << diagnosticLabel("Normal residual:")
                  << Alert::Notation::Number(normalResidual) << Alert::Raise;
    if (normalResidual > 1e-8)
      throw std::runtime_error("The rotated cut-face normals are inconsistent.");
  }

  template <class VectorSpace, class LinearSystem>
  void RotatedNitscheIntegrator::assembleVector(const VectorSpace& space,
    LinearSystem& system, Real diffusion, Real penalty,
    const FlatSet<Attribute>& fixedBoundaries) const
  {
    const auto& mesh = space.getMesh();
    const size_t faceDimension = mesh.getDimension() - 1;
    const IndexSet fixed = Internal::getFixedDOFs(space, fixedBoundaries);
    FaceNormal normal(mesh);
    for (const RotationPair& pair : RotationPairs)
    {
      std::vector<Eigen::Triplet<Real>> entries;
      for (auto face = mesh.getPolytope(faceDimension); face; ++face)
      {
        if (face->getAttribute() != pair.slave)
          continue;
        const auto& qf = QF::PolytopeQuadratureFormula::get(4, face->getGeometry());
        const auto& quadrature = face->getQuadrature(qf);
        for (size_t qp = 0; qp < quadrature.getSize(); ++qp)
        {
          const auto& slavePoint = quadrature.getPoint(qp);
          const auto mapped =
            m_locator.locate(pair.master, pair.rotation * slavePoint.vector());
          if (!mapped)
            throw std::runtime_error(
              "A rotated Nitsche quadrature point was not located.");
          const auto slave = Internal::evaluateScalarBasis<1>(
            space, *face, slavePoint.getPhysicalCoordinates());
          const auto master = Internal::evaluateScalarBasis<1>(
            space, mapped->getPolytope(), mapped->getPhysicalCoordinates());
          const Math::SpatialVector<Real> slaveNormal = normal.getValue(slavePoint);
          const Math::SpatialVector<Real> masterNormal = normal.getValue(*mapped);
          const Real weight = qf.getWeight(qp) * slavePoint.getDistortion();
          const Real stabilization = penalty * diffusion / Internal::diameter(*face);
          std::vector<Internal::VectorBasis> basis;
          basis.reserve(24);
          for (size_t node = 0; node < master.values.size(); ++node)
          {
            for (size_t component = 0; component < 3; ++component)
            {
              Internal::VectorBasis value{master.dofs[3 * node + component],
                Internal::zeroVector(), Internal::zeroVector()};
              value.jump(component) = master.values[node];
              value.flux(component) =
                0.5 * diffusion * master.gradients[node].dot(masterNormal);
              basis.emplace_back(std::move(value));
            }
          }
          for (size_t node = 0; node < slave.values.size(); ++node)
          {
            for (size_t component = 0; component < 3; ++component)
            {
              Internal::VectorBasis value{slave.dofs[3 * node + component],
                -slave.values[node] * pair.rotation.col(component),
                -0.5 * diffusion * slave.gradients[node].dot(slaveNormal) *
                  pair.rotation.col(component)};
              basis.emplace_back(std::move(value));
            }
          }
          for (const auto& test : basis)
          {
            for (const auto& trial : basis)
            {
              const Real coefficient = weight *
                (-trial.flux.dot(test.jump) - test.flux.dot(trial.jump) +
                  stabilization * trial.jump.dot(test.jump));
              if (std::abs(coefficient) > 1e-14)
                entries.emplace_back(test.dof, trial.dof, coefficient);
            }
          }
        }
      }
      Internal::addWithEliminatedColumns(system, entries, fixed);
    }
  }

  template <class ScalarSpace, class LinearSystem>
  void RotatedNitscheIntegrator::assembleScalarTracePenalty(const ScalarSpace& space,
    LinearSystem& system, Real penalty, const FlatSet<Attribute>& fixedBoundaries) const
  {
    const auto& mesh = space.getMesh();
    const size_t faceDimension = mesh.getDimension() - 1;
    const IndexSet fixed = Internal::getFixedDOFs(space, fixedBoundaries);
    for (const RotationPair& pair : RotationPairs)
    {
      std::vector<Eigen::Triplet<Real>> entries;
      for (auto face = mesh.getPolytope(faceDimension); face; ++face)
      {
        if (face->getAttribute() != pair.slave)
          continue;
        const auto& qf = QF::PolytopeQuadratureFormula::get(4, face->getGeometry());
        const auto& quadrature = face->getQuadrature(qf);
        for (size_t qp = 0; qp < quadrature.getSize(); ++qp)
        {
          const auto& slavePoint = quadrature.getPoint(qp);
          const auto mapped =
            m_locator.locate(pair.master, pair.rotation * slavePoint.vector());
          if (!mapped)
            throw std::runtime_error(
              "A rotated scalar quadrature point was not located.");
          const auto slave = Internal::evaluateScalarBasis<1>(
            space, *face, slavePoint.getPhysicalCoordinates());
          const auto master = Internal::evaluateScalarBasis<1>(
            space, mapped->getPolytope(), mapped->getPhysicalCoordinates());
          std::vector<std::pair<Index, Real>> basis;
          basis.reserve(master.values.size() + slave.values.size());
          for (size_t node = 0; node < master.values.size(); ++node)
            basis.emplace_back(master.dofs[node], master.values[node]);
          for (size_t node = 0; node < slave.values.size(); ++node)
            basis.emplace_back(slave.dofs[node], -slave.values[node]);
          const Real weight = qf.getWeight(qp) * slavePoint.getDistortion() * penalty *
            Internal::diameter(*face);
          for (const auto& [testDOF, testJump] : basis)
          {
            for (const auto& [trialDOF, trialJump] : basis)
            {
              const Real coefficient = weight * testJump * trialJump;
              if (std::abs(coefficient) > 1e-14)
                entries.emplace_back(testDOF, trialDOF, coefficient);
            }
          }
        }
      }
      Internal::addWithEliminatedColumns(system, entries, fixed);
    }
  }

  template <class ScalarSpace, class LinearSystem, class Reference>
  void RotatedNitscheIntegrator::assembleScalarTracePenalty(const ScalarSpace& space,
    LinearSystem& system, Real penalty, const Reference& reference) const
  {
    assembleScalarTracePenalty(space, system, penalty);
    auto& rhs = system.getVector();
    const auto& mesh = space.getMesh();
    const size_t faceDimension = mesh.getDimension() - 1;
    for (const RotationPair& pair : RotationPairs)
    {
      for (auto face = mesh.getPolytope(faceDimension); face; ++face)
      {
        if (face->getAttribute() != pair.slave)
          continue;
        const auto& qf = QF::PolytopeQuadratureFormula::get(4, face->getGeometry());
        const auto& quadrature = face->getQuadrature(qf);
        for (size_t qp = 0; qp < quadrature.getSize(); ++qp)
        {
          const auto& slavePoint = quadrature.getPoint(qp);
          const auto mapped =
            m_locator.locate(pair.master, pair.rotation * slavePoint.vector());
          if (!mapped)
            throw std::runtime_error(
              "A rotated scalar quadrature point was not located.");
          const auto slave = Internal::evaluateScalarBasis<1>(
            space, *face, slavePoint.getPhysicalCoordinates());
          const auto master = Internal::evaluateScalarBasis<1>(
            space, mapped->getPolytope(), mapped->getPhysicalCoordinates());
          const Real referenceJump =
            reference.getValue(*mapped) - reference.getValue(slavePoint);
          const Real weight = qf.getWeight(qp) * slavePoint.getDistortion() * penalty *
            Internal::diameter(*face);
          for (size_t node = 0; node < master.values.size(); ++node)
            rhs(master.dofs[node]) += weight * master.values[node] * referenceJump;
          for (size_t node = 0; node < slave.values.size(); ++node)
            rhs(slave.dofs[node]) -= weight * slave.values[node] * referenceJump;
        }
      }
    }
  }

  template <class U0, class U1, class U2>
  Real RotatedNitscheIntegrator::familyJump(
    const U0& u0, const U1& u1, const U2& u2) const
  {
    const auto& mesh = u0.getFiniteElementSpace().getMesh();
    const size_t faceDimension = mesh.getDimension() - 1;
    Real residual = 0;
    for (const RotationPair& pair : RotationPairs)
    {
      for (auto face = mesh.getPolytope(faceDimension); face; ++face)
      {
        if (face->getAttribute() != pair.slave)
          continue;
        const auto& qf = QF::PolytopeQuadratureFormula::get(4, face->getGeometry());
        const auto& quadrature = face->getQuadrature(qf);
        for (size_t qp = 0; qp < quadrature.getSize(); ++qp)
        {
          const auto& point = quadrature.getPoint(qp);
          const auto mapped =
            m_locator.locate(pair.master, pair.rotation * point.vector());
          if (!mapped)
            throw std::runtime_error("A rotated trace diagnostic point was not located.");
          Math::SpatialMatrix<Real> slave(3, 3), master(3, 3);
          Internal::setColumn(slave, 0, u0.getValue(point));
          Internal::setColumn(slave, 1, u1.getValue(point));
          Internal::setColumn(slave, 2, u2.getValue(point));
          Internal::setColumn(master, 0, u0.getValue(*mapped));
          Internal::setColumn(master, 1, u1.getValue(*mapped));
          Internal::setColumn(master, 2, u2.getValue(*mapped));
          const auto jump = master - pair.rotation * slave * pair.rotation.transpose();
          for (size_t i = 0; i < 3; ++i)
            for (size_t j = 0; j < 3; ++j)
              residual = std::max(residual, std::abs(jump(i, j)));
        }
      }
    }
    return residual;
  }

  template <class U>
  Real RotatedNitscheIntegrator::vectorJump(const U& u) const
  {
    const auto& mesh = u.getFiniteElementSpace().getMesh();
    const size_t faceDimension = mesh.getDimension() - 1;
    Real residual = 0;
    for (const RotationPair& pair : RotationPairs)
    {
      for (auto face = mesh.getPolytope(faceDimension); face; ++face)
      {
        if (face->getAttribute() != pair.slave)
          continue;
        const auto& qf = QF::PolytopeQuadratureFormula::get(4, face->getGeometry());
        const auto& quadrature = face->getQuadrature(qf);
        for (size_t qp = 0; qp < quadrature.getSize(); ++qp)
        {
          const auto& point = quadrature.getPoint(qp);
          const auto mapped =
            m_locator.locate(pair.master, pair.rotation * point.vector());
          if (!mapped)
            throw std::runtime_error("A rotated trace diagnostic point was not located.");
          const auto jump = u.getValue(*mapped) - pair.rotation * u.getValue(point);
          for (size_t i = 0; i < 3; ++i)
            residual = std::max(residual, std::abs(jump(i)));
        }
      }
    }
    return residual;
  }

  template <class U>
  Real RotatedNitscheIntegrator::scalarJump(const U& u) const
  {
    const auto& mesh = u.getFiniteElementSpace().getMesh();
    const size_t faceDimension = mesh.getDimension() - 1;
    Real residual = 0;
    for (const RotationPair& pair : RotationPairs)
    {
      for (auto face = mesh.getPolytope(faceDimension); face; ++face)
      {
        if (face->getAttribute() != pair.slave)
          continue;
        const auto& qf = QF::PolytopeQuadratureFormula::get(4, face->getGeometry());
        const auto& quadrature = face->getQuadrature(qf);
        for (size_t qp = 0; qp < quadrature.getSize(); ++qp)
        {
          const auto& point = quadrature.getPoint(qp);
          const auto mapped =
            m_locator.locate(pair.master, pair.rotation * point.vector());
          if (!mapped)
            throw std::runtime_error(
              "A rotated scalar diagnostic point was not located.");
          residual =
            std::max(residual, std::abs(u.getValue(*mapped) - u.getValue(point)));
        }
      }
    }
    return residual;
  }

}

#endif
