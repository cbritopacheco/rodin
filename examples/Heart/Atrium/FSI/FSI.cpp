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
      m_previousMeshDisplacement(fluid.m_vh),
      m_importedTraction(solid.m_space),
      m_trialTraction(solid.m_space),
      m_importedVelocity(solid.m_space),
      m_interfaceVelocity(solid.m_space),
      m_previousRobinData(solid.m_space),
      m_previousResidual(solid.m_space),
      m_residualDifference(solid.m_space),
      m_tractionResidual(solid.m_space),
      m_previousTractionResidual(solid.m_space),
      m_tractionResidualDifference(solid.m_space),
      m_robinCoefficient(parameters.robinCoefficient)
  {
    if (parameters.maximumIterations < 1 ||
      !std::isfinite(parameters.relativeTolerance) || parameters.relativeTolerance <= 0 ||
      !std::isfinite(parameters.absoluteTolerance) || parameters.absoluteTolerance <= 0 ||
      !std::isfinite(m_robinCoefficient) || m_robinCoefficient < 0 ||
      !std::isfinite(parameters.initialRelaxation) ||
      !std::isfinite(parameters.minimumRelaxation) ||
      !std::isfinite(parameters.maximumRelaxation) || parameters.minimumRelaxation <= 0 ||
      parameters.minimumRelaxation > parameters.initialRelaxation ||
      parameters.initialRelaxation > parameters.maximumRelaxation ||
      parameters.maximumRelaxation > 1)
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
    m_previousTractionResidual = zero;
    m_tractionResidualDifference = zero;
    m_importedVelocity = zero;
    m_previousRobinData = zero;
    m_previousResidual = zero;
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
    // The Robin stress update is the fluid variational reaction, including
    // the interface convective correction. All three traces are P1.
    m_trialTraction = m_importedTraction +
      m_solid.m_robinCoefficient * (m_importedVelocity - m_interfaceVelocity);
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
    bool converged = false;
    Real relaxation = m_parameters.initialRelaxation;
    for (size_t iteration = 1; iteration <= m_parameters.maximumIterations; ++iteration)
    {
      Stage coupling("Coupling iteration " + std::to_string(iteration), m_fluid.isRoot());
      if (m_fluid.isRoot())
      {
        Log("FSI/Coupling") << "iteration=" << iteration << '/'
                            << m_parameters.maximumIterations
                            << " relaxation=" << relaxation << Alert::Raise;
      }
      m_previousRobinData.setData(m_solid.m_robinData.getData());
      m_solid.solveStep();
      moveFluidMesh();
      if (!m_fluid.solveFlow())
        Alert::Exception() << "FSI fluid velocity exceeded its limit" << Alert::Raise;
      updateTraction();
      Stage interface("Interface convergence and relaxation", m_fluid.isRoot());
      PetscReal difference = 0, norm = 0, tractionDifference = 0, tractionNorm = 0;
      m_tractionResidual = m_importedTraction - m_trialTraction;
      m_previousRobinData.axpy(-1, m_solid.m_robinData);
      PetscErrorCode ierr =
        VecNorm(m_previousRobinData.getData(), NORM_INFINITY, &difference);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecNorm(m_solid.m_robinData.getData(), NORM_INFINITY, &norm);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecNorm(m_tractionResidual.getData(), NORM_INFINITY, &tractionDifference);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecNorm(m_trialTraction.getData(), NORM_INFINITY, &tractionNorm);
      assert(ierr == PETSC_SUCCESS);
      (void)ierr;
      m_diagnostics.iterations = iteration;
      m_diagnostics.robinResidual = 1e3 * std::max(difference, tractionDifference);
      m_diagnostics.robinTolerance = m_parameters.absoluteTolerance +
        m_parameters.relativeTolerance * 1e3 * std::max(norm, tractionNorm);
      if (m_fluid.isRoot())
      {
        Log("FSI/Coupling") << "iteration=" << iteration
                            << " Robin residual=" << m_diagnostics.robinResidual
                            << " Pa tolerance=" << m_diagnostics.robinTolerance << " Pa"
                            << Alert::Raise;
      }
      if (m_diagnostics.robinResidual <= m_diagnostics.robinTolerance)
      {
        m_importedTraction.setData(m_trialTraction.getData());
        transfer(m_solidToFluid, m_importedTraction, m_fluid.m_wallTraction, 1e3);
        converged = true;
        break;
      }
      if (iteration > 1)
      {
        m_residualDifference.setData(m_previousRobinData.getData());
        m_residualDifference.axpy(-1, m_previousResidual);
        m_tractionResidualDifference = m_tractionResidual - m_previousTractionResidual;
        PetscScalar numerator = 0, denominator = 0, tractionNumerator = 0,
                    tractionDenominator = 0;
        ierr = VecDot(
          m_previousResidual.getData(), m_residualDifference.getData(), &numerator);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecDot(
          m_residualDifference.getData(), m_residualDifference.getData(), &denominator);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecDot(m_previousTractionResidual.getData(),
          m_tractionResidualDifference.getData(), &tractionNumerator);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecDot(m_tractionResidualDifference.getData(),
          m_tractionResidualDifference.getData(), &tractionDenominator);
        assert(ierr == PETSC_SUCCESS);
        numerator += tractionNumerator;
        denominator += tractionDenominator;
        if (PetscRealPart(denominator) > 0)
        {
          const Real proposal =
            -relaxation * PetscRealPart(numerator) / PetscRealPart(denominator);
          if (std::isfinite(proposal))
          {
            relaxation = std::clamp(
              proposal, m_parameters.minimumRelaxation, m_parameters.maximumRelaxation);
          }
        }
      }
      m_previousResidual.setData(m_previousRobinData.getData());
      m_previousTractionResidual.setData(m_tractionResidual.getData());
      m_solid.m_robinData.axpy(1 - relaxation, m_previousRobinData);
      m_trialTraction.axpy(1 - relaxation, m_tractionResidual);
      m_importedTraction.setData(m_trialTraction.getData());
      transfer(m_solidToFluid, m_importedTraction, m_fluid.m_wallTraction, 1e3);
      m_diagnostics.relaxation = relaxation;
    }
    if (!converged)
      Alert::Exception() << "Atrium interface iteration did not converge" << Alert::Raise;
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
      Log("FSI/State", false, Log::Status::Success)
        << "converged: coupling_iterations=" << m_diagnostics.iterations
        << " velocity_mismatch=" << m_diagnostics.velocityMismatch << " m/s"
        << " displacement_mismatch=" << m_diagnostics.displacementMismatch << " m"
        << " min_mesh_J=" << m_diagnostics.minimumMeshJacobian << Alert::Raise;
    }
  }
}
