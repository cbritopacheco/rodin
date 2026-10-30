/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <chrono>
#include <cmath>
#include <limits>
#include <boost/mpi/collectives.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include "FSI.h"
#include "../Fluid/Output.h"
#include "../Log.h"

namespace Rodin::Examples::Heart::Atrium
{
  using namespace Geometry;
  using namespace Variational;

  template <class SourceField, class TargetField>
  ::VecScatter FSI::makeScatter(const Fluid::MeshType& sourceMesh,
    const SourceField& source, const Fluid::MeshType& targetMesh,
    const TargetField& target)
  {
    const auto nodes = [](const Fluid::MeshType& mesh) {
      IndexSet vertices;
      if (mesh.getDimension() == 2)
      {
        for (auto vertex = mesh.getVertex(); vertex; ++vertex)
          vertices.insert(vertex->getIndex());
      }
      else
      {
        for (auto face = mesh.getBoundary(); face; ++face)
        {
          if (face->getAttribute() == 4)
            vertices.insert(face->getVertices().begin(), face->getVertices().end());
        }
      }
      return vertices;
    };
    using Record = std::pair<Index, std::array<PetscInt, 3>>;
    std::vector<Record> local;
    const auto& sourceSpace = source.getFiniteElementSpace();
    for (const Index vertex : nodes(sourceMesh))
    {
      Record record;
      record.first = sourceMesh.getShard().getPolytopeMap(0).left.at(vertex);
      for (Index c = 0; c < 3; ++c)
        record.second[c] = sourceSpace.getGlobalIndex({0, vertex}, c);
      local.push_back(record);
    }
    std::vector<std::vector<Record>> gathered;
    const auto& comm = sourceMesh.getContext().getCommunicator();
    boost::mpi::all_gather(comm, local, gathered);
    UnorderedMap<Index, std::array<PetscInt, 3>> indices;
    for (const auto& records : gathered)
    {
      for (const auto& record : records)
        indices.emplace(record);
    }
    std::vector<PetscInt> from, to;
    int missing = 0;
    const auto& targetSpace = target.getFiniteElementSpace();
    const auto targetNodes = nodes(targetMesh);
    // An interface vertex owner need not own an interface face.
    // Global trace membership must therefore select owned target vertices.
    for (auto node = targetMesh.getVertex(); node; ++node)
    {
      const Index vertex = node->getIndex();
      if (!targetMesh.getShard().isOwned(0, vertex))
        continue;
      const Index parent = targetMesh.getShard().getPolytopeMap(0).left.at(vertex);
      const auto found = indices.find(parent);
      if (found == indices.end())
      {
        if (targetNodes.count(vertex))
          missing = 1;
        continue;
      }
      for (Index c = 0; c < 3; ++c)
      {
        from.push_back(found->second[c]);
        to.push_back(targetSpace.getGlobalIndex({0, vertex}, c));
      }
    }
    const auto count = boost::mpi::all_reduce(comm, to.size() / 3, std::plus<size_t>());
    if (count != indices.size())
      missing = 1;
    if (boost::mpi::all_reduce(comm, missing, std::plus<int>()) != 0)
    {
      Alert::Exception() << "The solid and fluid interfaces do not share parent vertices"
                         << Alert::Raise;
    }
    ::IS sourceIndices = nullptr, targetIndices = nullptr;
    PetscErrorCode ierr =
      ISCreateGeneral(comm, from.size(), from.data(), PETSC_COPY_VALUES, &sourceIndices);
    assert(ierr == PETSC_SUCCESS);
    ierr = ISCreateGeneral(comm, to.size(), to.data(), PETSC_COPY_VALUES, &targetIndices);
    assert(ierr == PETSC_SUCCESS);
    ::VecScatter result = nullptr;
    ierr = VecScatterCreate(
      source.getData(), sourceIndices, target.getData(), targetIndices, &result);
    assert(ierr == PETSC_SUCCESS);
    ierr = ISDestroy(&sourceIndices);
    assert(ierr == PETSC_SUCCESS);
    ierr = ISDestroy(&targetIndices);
    assert(ierr == PETSC_SUCCESS);
    (void)ierr;
    return result;
  }

  template <class SourceField, class TargetField>
  void FSI::transfer(::VecScatter plan, const SourceField& source, TargetField& target,
    Real scale, ::ScatterMode direction)
  {
    source.flush();
    static_cast<const TargetField&>(target).flush();
    target.flush();
    PetscErrorCode ierr = VecSet(target.getData(), 0.0);
    assert(ierr == PETSC_SUCCESS);
    ierr =
      VecScatterBegin(plan, source.getData(), target.getData(), INSERT_VALUES, direction);
    assert(ierr == PETSC_SUCCESS);
    ierr =
      VecScatterEnd(plan, source.getData(), target.getData(), INSERT_VALUES, direction);
    assert(ierr == PETSC_SUCCESS);
    target *= scale;
    (void)ierr;
  }

  FSI::FSI(Solid& solid, Fluid& fluid, const Parameters& parameters)
    : m_solid(solid),
      m_fluid(fluid),
      m_parameters(parameters),
      m_motionTrial(fluid.m_vh),
      m_motionTest(fluid.m_vh),
      m_motion(m_motionTrial, m_motionTest),
      m_motionSolver(m_motion),
      m_tractionLoad(fluid.m_v),
      m_tractionWeights(fluid.m_v),
      m_previousMeshDisplacement(fluid.m_vh),
      m_importedTraction(solid.m_space),
      m_trialTraction(solid.m_space),
      m_importedVelocity(solid.m_space),
      m_interfaceVelocity(solid.m_space),
      m_previousRobinData(solid.m_space),
      m_previousDisplacementIteration(solid.m_space),
      m_residualDifference(solid.m_space),
      m_tractionResidual(solid.m_space),
      m_robinCoefficient(parameters.robinCoefficient)
  {
    if (parameters.maximumIterations < 1 ||
      !std::isfinite(parameters.relativeTolerance) || parameters.relativeTolerance < 0 ||
      !std::isfinite(m_robinCoefficient) || m_robinCoefficient < 0)
      Alert::Exception() << "Invalid FSI parameters" << Alert::Raise;
    if (fluid.m_cfg.useVentricularModel || fluid.m_cfg.useVMS ||
      fluid.m_cfg.solveKinetics)
    {
      Alert::Exception() << "Atrium FSI currently uses Oseen SUPG/PSPG flow without "
                            "0D ventricular coupling, orthogonal subscales, or species"
                         << Alert::Raise;
    }
    fluid.m_coupled = true;
    fluid.m_cfg.inletPressureMean = 0;
    fluid.m_cfg.inletPressureAmplitude = 0;
    fluid.m_cfg.prescribedOutletPressure = 0;
    fluid.m_cfg.prescribedOutletResistance = 0;
    fluid.m_cfg.inletImpedance = 0;
    fluid.m_cfg.inletTangentialDamping = 0;
    if (m_robinCoefficient == 0)
    {
      // Matrix-based isotropic scale; fiber reinforcement is not isotropic.
      constexpr Real KilopascalsToPascals = 1e3;
      const auto& passive = solid.m_parameters.passive;
      const Real shear = KilopascalsToPascals * passive.matrixScale;
      const Real bulk = KilopascalsToPascals * 2 * passive.volumePenalty;
      const Real young = 9 * bulk * shear / (3 * bulk + shear);
      m_robinCoefficient = std::sqrt(solid.m_parameters.density * young);
      if (!(m_robinCoefficient > 0))
      {
        Alert::Exception() << "A positive FSI Robin coefficient is required"
                           << Alert::Raise;
      }
    }
    fluid.m_robinCoefficient = m_robinCoefficient;
    m_diagnostics.robinCoefficient = m_robinCoefficient;
    // Rebuild the retained form after setting its constant boundary coefficients.
    fluid.setupFlow();
    // Pa s/m -> kPa s/mm. Both sides use the same reference-area impedance.
    solid.m_robinCoefficient = 1e-6 * m_robinCoefficient;
    const Real interfaceVelocityDerivative = 1 / solid.m_timeStep;
    const auto velocity =
      interfaceVelocityDerivative * (solid.m_displacement - solid.m_previousDisplacement);
    solid.m_mechanics += solid.m_robinCoefficient * interfaceVelocityDerivative *
      BoundaryIntegral(solid.m_increment, solid.m_test).over(Attribute{4});
    solid.m_mechanics += BoundaryIntegral(
      solid.m_robinCoefficient * velocity - solid.m_robinData, solid.m_test)
                           .over(Attribute{4});

    const VectorFunction zero(
      3, [](const Point&) { return Math::SpatialVector<Real>{0, 0, 0}; });
    m_motion = Integral(Jacobian(m_motionTrial), Jacobian(m_motionTest)) +
      DirichletBC(m_motionTrial, fluid.m_wallDisplacement).on(fluid.m_cfg.wall) +
      DirichletBC(m_motionTrial, zero)
        .on(Attribute{1}, Attribute{2}, Attribute{3}, Attribute{5}, Attribute{6});
    m_motionSolver.setPrefix(std::string("atrium_ale_")).setType(KSPCG);
    m_motionSolver.setTolerances(1e-12, 1e-30, 1e5, 10000);
    m_referenceVertices.resize(fluid.m_mesh.getShard().getVertexCount());
    for (auto vertex = fluid.m_mesh.getVertex(); vertex; ++vertex)
      m_referenceVertices[vertex->getIndex()] = vertex->getCoordinates();
    m_referenceCellDeterminants.resize(fluid.m_mesh.getShard().getPolytopeCount(3));
    m_previousCellMeasures.resize(m_referenceCellDeterminants.size());
    fluid.m_massRatios.assign(m_referenceCellDeterminants.size(), 1);
    for (auto cell = fluid.m_mesh.getCell(); cell; ++cell)
    {
      const Point point(*cell, Math::SpatialPoint{.25, .25, .25});
      m_referenceCellDeterminants[cell->getIndex()] = point.getJacobian().determinant();
      m_previousCellMeasures[cell->getIndex()] = cell->getMeasure();
    }
    m_referenceFaceMeasures.resize(fluid.m_mesh.getShard().getPolytopeCount(2));
    fluid.m_inverseInterfaceStretches.assign(m_referenceFaceMeasures.size(), 1);
    for (auto face = fluid.m_mesh.getFace(); face; ++face)
      m_referenceFaceMeasures[face->getIndex()] = face->getMeasure();
    m_previousMeshDisplacement = zero;
    m_importedTraction = zero;
    m_trialTraction = zero;
    m_interfaceVelocity = zero;
    m_tractionResidual = zero;
    m_importedVelocity = zero;
    m_previousRobinData = zero;
    m_previousDisplacementIteration = zero;
    m_residualDifference = zero;
    m_solidToFluid = makeScatter(
      solid.m_mesh, solid.m_displacement, fluid.m_mesh, fluid.m_wallDisplacement);
  }

  FSI::~FSI()
  {
    if (m_solidToFluid)
    {
      const PetscErrorCode ierr = VecScatterDestroy(&m_solidToFluid);
      assert(ierr == PETSC_SUCCESS);
      (void)ierr;
    }
  }

  void FSI::initialize()
  {
    moveFluidMesh();
    m_previousMeshDisplacement.setData(m_fluid.m_meshDisplacement.getData());
    m_fluid.m_meshVelocity = Math::SpatialVector<Real>{0, 0, 0};
    m_fluid.m_wallVelocity = Math::SpatialVector<Real>{0, 0, 0};
    m_fluid.m_output->write(0);
  }

  void FSI::restoreReference()
  {
    for (Index i = 0; i < m_referenceVertices.size(); ++i)
      m_fluid.m_mesh.setVertexCoordinates(i, m_referenceVertices[i]);
  }

  void FSI::moveFluidMesh()
  {
    Stage stage("ALE mesh motion", m_fluid.isRoot());
    transfer(m_solidToFluid, m_solid.m_displacement, m_fluid.m_wallDisplacement, 1e-3);
    restoreReference();
    m_motion.assemble();
    m_motion.solve(m_motionSolver);
    ::KSPConvergedReason reason;
    const PetscErrorCode ierr =
      KSPGetConvergedReason(m_motionSolver.getHandle(), &reason);
    assert(ierr == PETSC_SUCCESS);
    (void)ierr;
    if (reason <= 0)
      Alert::Exception() << "Harmonic ALE extension did not converge" << Alert::Raise;
    m_fluid.m_meshDisplacement.setData(m_motionTrial.getSolution().getData());
    m_fluid.m_meshVelocity =
      (m_fluid.m_meshDisplacement - m_previousMeshDisplacement) / m_fluid.m_run.timeStep;
    m_fluid.m_wallVelocity.setData(m_fluid.m_meshVelocity.getData());
    std::vector<Math::SpatialPoint> moved(m_referenceVertices.size());
    const Polytope::Traits vertexTraits(Polytope::Type::Point);
    for (auto vertex = m_fluid.m_mesh.getVertex(); vertex; ++vertex)
    {
      const Index i = vertex->getIndex();
      const Point point(*vertex, vertexTraits.getVertex(0), m_referenceVertices[i]);
      moved[i] = m_referenceVertices[i] + m_fluid.m_meshDisplacement(point);
    }
    for (Index i = 0; i < moved.size(); ++i)
      m_fluid.m_mesh.setVertexCoordinates(i, moved[i]);
    for (auto face = m_fluid.m_mesh.getBoundary(); face; ++face)
    {
      if (face->getAttribute() == m_fluid.m_cfg.wall)
      {
        m_fluid.m_inverseInterfaceStretches[face->getIndex()] =
          m_referenceFaceMeasures[face->getIndex()] / face->getMeasure();
      }
    }
    Real minimum = std::numeric_limits<Real>::infinity();
    for (auto cell = m_fluid.m_mesh.getCell(); cell; ++cell)
    {
      const Index i = cell->getIndex();
      const Point point(*cell, Math::SpatialPoint{.25, .25, .25});
      const Real determinant = point.getJacobian().determinant();
      if (!std::isfinite(determinant) || m_referenceCellDeterminants[i] == 0)
        minimum = -std::numeric_limits<Real>::infinity();
      else
        minimum = std::min(minimum, determinant / m_referenceCellDeterminants[i]);
      m_fluid.m_massRatios[i] = m_previousCellMeasures[i] / (std::abs(determinant) / 6);
    }
    m_diagnostics.minimumMeshJacobian =
      boost::mpi::all_reduce(m_fluid.m_mesh.getContext().getCommunicator(), minimum,
        boost::mpi::minimum<Real>());
    if (!(m_diagnostics.minimumMeshJacobian > 0))
      Alert::Exception() << "ALE mesh has a nonpositive Jacobian" << Alert::Raise;
  }

  void FSI::updateTraction()
  {
    Stage stage("Interface traction exchange", m_fluid.isRoot());
    transfer(m_solidToFluid, m_fluid.m_u.getSolution(), m_importedVelocity, 1e3,
      SCATTER_REVERSE);
    transfer(
      m_solidToFluid, m_fluid.m_wallVelocity, m_interfaceVelocity, 1e3, SCATTER_REVERSE);
    m_fluid.recoverGradient();
    const auto normal = BoundaryNormal(m_fluid.m_mesh);
    const auto& cy = m_fluid.m_cfg.viscosity;
    const VectorFunction traction(3, [&](const Point& point) {
      Math::SpatialMatrix<Real> gradient(3, 3);
      const auto g0 = m_fluid.m_gradRec0(point);
      const auto g1 = m_fluid.m_gradRec1(point);
      const auto g2 = m_fluid.m_gradRec2(point);
      for (Index j = 0; j < 3; ++j)
      {
        gradient(0, j) = g0(j);
        gradient(1, j) = g1(j);
        gradient(2, j) = g2(j);
      }
      const Math::SpatialMatrix<Real> strain = .5 * (gradient + gradient.transpose());
      const Real shear = std::sqrt(
        cy.gammaRegularization * cy.gammaRegularization + 2 * strain.squaredNorm());
      const Real viscosity = cy.muInf +
        (cy.mu0 - cy.muInf) *
          std::pow(1 + std::pow(cy.lambda * shear, cy.yasuda), (cy.n - 1) / cy.yasuda);
      const Math::SpatialVector<Real> n = normal.getValue(point);
      Math::SpatialVector<Real> value =
        m_fluid.m_p.getSolution()(point) * n - 2 * viscosity * strain * n;
      for (Index j = 0; j < 3; ++j)
      {
        if (!std::isfinite(value(j)))
          Alert::Exception() << "Non-finite recovered interface traction" << Alert::Raise;
      }
      return value;
    });
    const VectorFunction ones(
      3, [](const Point&) { return Math::SpatialVector<Real>{1, 1, 1}; });
    const RealFunction inverseStretch([&](const Point& point) {
      return m_fluid.m_inverseInterfaceStretches.at(point.getPolytope().getIndex());
    });
    m_tractionWeights =
      BoundaryIntegral(inverseStretch * ones, m_fluid.m_v).over(m_fluid.m_cfg.wall);
    m_tractionWeights.assemble();
    // Lumped reference trace mass averages corner normals across ranks and
    // preserves the assembled current-area resultant on the solid reference wall.
    auto& target = m_fluid.m_wallTraction;
    auto load = BoundaryIntegral(traction, m_fluid.m_v);
    load.setOrder(4);
    m_tractionLoad = load.over(m_fluid.m_cfg.wall);
    m_tractionLoad.assemble();
    target.flush();
    const PetscScalar *weights = nullptr, *values = nullptr;
    PetscScalar* result = nullptr;
    PetscInt size = 0;
    PetscErrorCode ierr = VecGetLocalSize(target.getData(), &size);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecGetArrayRead(m_tractionWeights.getVector(), &weights);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecGetArrayRead(m_tractionLoad.getVector(), &values);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecGetArray(target.getData(), &result);
    assert(ierr == PETSC_SUCCESS);
    for (PetscInt i = 0; i < size; ++i)
      result[i] = PetscRealPart(weights[i]) > 0 ? values[i] / weights[i] : 0;
    ierr = VecRestoreArray(target.getData(), &result);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecRestoreArrayRead(m_tractionLoad.getVector(), &values);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecRestoreArrayRead(m_tractionWeights.getVector(), &weights);
    assert(ierr == PETSC_SUCCESS);
    target *= Real(1);
    (void)ierr;
    transfer(
      m_solidToFluid, m_fluid.m_wallTraction, m_trialTraction, 1e-3, SCATTER_REVERSE);
    m_solid.m_robinData =
      m_trialTraction + m_solid.m_robinCoefficient * m_importedVelocity;
  }

  void FSI::step()
  {
    Stage stage("FSI", m_fluid.isRoot());
    m_diagnostics = Diagnostics{};
    m_diagnostics.robinCoefficient = m_robinCoefficient;
    m_previousMeshDisplacement.setData(m_fluid.m_meshDisplacement.getData());
    for (auto cell = m_fluid.m_mesh.getCell(); cell; ++cell)
      m_previousCellMeasures[cell->getIndex()] = cell->getMeasure();
    m_solid.beginStep();
    m_fluid.beginStep();
    m_previousDisplacementIteration.setData(m_solid.m_displacement.getData());
    for (size_t iteration = 1; iteration <= m_parameters.maximumIterations; ++iteration)
    {
      Stage coupling("Coupling iteration " + std::to_string(iteration), m_fluid.isRoot());
      if (m_fluid.isRoot())
      {
        Log("FSI/Coupling") << "iteration=" << iteration << '/'
                            << m_parameters.maximumIterations << Alert::Raise;
      }
      m_previousRobinData.setData(m_solid.m_robinData.getData());
      m_solid.solveStep();
      m_residualDifference = m_solid.m_displacement - m_previousDisplacementIteration;
      const Real change = m_residualDifference.norm(NORM_2);
      const Real state = m_solid.m_displacement.norm(NORM_2);
      m_diagnostics.displacementChange = state > 0 ? change / state : 1e-3 * change;
      m_previousDisplacementIteration.setData(m_solid.m_displacement.getData());
      moveFluidMesh();
      if (!m_fluid.solveFlow())
        Alert::Exception() << "FSI fluid velocity exceeded its limit" << Alert::Raise;
      updateTraction();
      Stage interface("Interface convergence and traction diagnostics", m_fluid.isRoot());
      PetscReal difference = 0, tractionDifference = 0;
      m_tractionResidual = m_importedTraction - m_trialTraction;
      m_previousRobinData.axpy(-1, m_solid.m_robinData);
      PetscErrorCode ierr =
        VecNorm(m_previousRobinData.getData(), NORM_INFINITY, &difference);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecNorm(m_tractionResidual.getData(), NORM_INFINITY, &tractionDifference);
      assert(ierr == PETSC_SUCCESS);
      (void)ierr;
      m_diagnostics.iterations = iteration;
      m_diagnostics.tractionDataUpdate = 1e3 * std::max(difference, tractionDifference);
      m_importedTraction.setData(m_trialTraction.getData());
      m_diagnostics.converged = iteration > 1 &&
        m_diagnostics.displacementChange < m_parameters.relativeTolerance;
      if (m_fluid.isRoot())
      {
        Log("FSI/Coupling") << "iteration=" << iteration
                            << " displacement_change=" << m_diagnostics.displacementChange
                            << " tolerance=" << m_parameters.relativeTolerance
                            << " traction_data_update="
                            << m_diagnostics.tractionDataUpdate << " Pa" << Alert::Raise;
      }
      if (m_diagnostics.converged)
        break;
    }
    m_diagnostics.budgetExhausted = !m_diagnostics.converged;
    {
      Stage stage("Accepted interface checks", m_fluid.isRoot());
    Real velocityMismatch = 0, displacementMismatch = 0, newmarkMismatch = 0;
    transfer(m_solidToFluid, m_fluid.m_meshDisplacement, m_residualDifference, 1e3,
      SCATTER_REVERSE);
    const Polytope::Traits vertexTraits(Polytope::Type::Point);
    for (auto face = m_solid.m_mesh.getBoundary(); face; ++face)
    {
      if (face->getAttribute() != 4)
        continue;
      for (const Index vertex : face->getVertices())
      {
        const Point point(*m_solid.m_mesh.getVertex(vertex), vertexTraits.getVertex(0));
        velocityMismatch = std::max(velocityMismatch,
          1e-3 * (m_interfaceVelocity(point) - m_importedVelocity(point)).norm());
        newmarkMismatch = std::max(newmarkMismatch,
          1e-3 * (m_solid.m_velocity(point) - m_importedVelocity(point)).norm());
        displacementMismatch = std::max(displacementMismatch,
          1e-3 * (m_solid.m_displacement(point) - m_residualDifference(point)).norm());
      }
    }
    for (auto face = m_fluid.m_mesh.getBoundary(); face; ++face)
    {
      if (face->getAttribute() != m_fluid.m_cfg.wall)
        continue;
      for (const Index vertex : face->getVertices())
      {
        const auto v = m_fluid.m_mesh.getVertex(vertex);
        const Point point(*v, vertexTraits.getVertex(0));
        displacementMismatch = std::max(displacementMismatch,
          (v->getCoordinates() - m_referenceVertices[vertex] -
            m_fluid.m_wallDisplacement(point))
            .norm());
      }
    }
    const auto& comm = m_fluid.m_mesh.getContext().getCommunicator();
    m_diagnostics.velocityMismatch =
      boost::mpi::all_reduce(comm, velocityMismatch, boost::mpi::maximum<Real>());
    m_diagnostics.newmarkVelocityMismatch =
      boost::mpi::all_reduce(comm, newmarkMismatch, boost::mpi::maximum<Real>());
    m_diagnostics.displacementMismatch =
      boost::mpi::all_reduce(comm, displacementMismatch, boost::mpi::maximum<Real>());
    const auto normal = BoundaryNormal(m_fluid.m_mesh);
    m_fluid.m_flux =
      BoundaryIntegral(Dot(m_fluid.m_u.getSolution(), normal), m_fluid.m_qFlux)
        .over(m_fluid.m_cfg.wall);
    m_fluid.m_flux.assemble();
    m_diagnostics.wallFlux = m_fluid.m_flux(m_fluid.m_one);
    }
    m_fluid.computeWallShear();
    m_solid.commitStep();
    m_fluid.commitStep();
    if (m_fluid.isRoot())
    {
      Log("FSI/State", false,
        m_diagnostics.converged ? Log::Status::Success : Log::Status::Retry)
        << (m_diagnostics.converged ? "converged"
                                    : "budget exhausted; last sweep accepted")
        << ": coupling_iterations=" << m_diagnostics.iterations
        << " velocity_mismatch=" << m_diagnostics.velocityMismatch << " m/s"
        << " displacement_mismatch=" << m_diagnostics.displacementMismatch << " m"
        << " min_mesh_J=" << m_diagnostics.minimumMeshJacobian << Alert::Raise;
    }
  }
}
