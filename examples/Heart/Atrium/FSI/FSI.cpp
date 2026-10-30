/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <chrono>
#include <cmath>
#include <limits>
#include <Eigen/Dense>
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
      m_previousDisplacementIteration(solid.m_space),
      m_residualDifference(solid.m_space),
      m_tractionResidual(solid.m_space),
      m_robinCoefficient(parameters.robinCoefficient)
  {
    if (parameters.maximumIterations < 1 ||
      !std::isfinite(parameters.relativeTolerance) || parameters.relativeTolerance < 0 ||
      !std::isfinite(m_robinCoefficient) || m_robinCoefficient < 0)
      Alert::Exception() << "Invalid FSI parameters" << Alert::Raise;
    if (fluid.m_cfg.useVentricularModel || fluid.m_cfg.solveKinetics)
    {
      Alert::Exception() << "Atrium FSI currently excludes 0D ventricular coupling "
                            "and species transport"
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
    const ::Vec layout = solid.m_robinData.sync().getData();
    for (::Vec* v : {&m_inputData, &m_inputTraction, &m_outputData, &m_outputTraction,
           &m_residualData, &m_residualTraction})
    {
      const PetscErrorCode ierr = VecDuplicate(layout, v);
      assert(ierr == PETSC_SUCCESS);
      (void)ierr;
    }
  }

  void FSI::clearHistory(bool all)
  {
    // Pairs from the current and the previous `quasiNewtonReuse` steps remain.
    for (size_t i = m_secants.size(); i-- > 0;)
    {
      if (!all && m_secantSteps[i] + m_parameters.quasiNewtonReuse >= m_stepIndex)
        continue;
      for (::Vec& v : m_secants[i])
      {
        const PetscErrorCode ierr = VecDestroy(&v);
        assert(ierr == PETSC_SUCCESS);
        (void)ierr;
      }
      m_secants.erase(m_secants.begin() + i);
      m_secantSteps.erase(m_secantSteps.begin() + i);
    }
    m_haveSecantBase = false;
  }

  void FSI::accelerate()
  {
    // IQN-ILS (Degroote et al. 2009) on x = (g, T): with sweep input x^k,
    // output H(x^k) and residual r^k = H(x^k) - x^k, the columns
    // V = [r^k - r^{k-1}, ...], W = [H(x^k) - H(x^{k-1}), ...] give
    // c = argmin |V c + r^k| and x^{k+1} = H(x^k) + W c. Both blocks are in kPa.
    const ::Vec gOut = m_solid.m_robinData.sync().getData();
    const ::Vec tOut = m_trialTraction.sync().getData();
    PetscErrorCode ierr;
    std::array<::Vec, 4> column{};
    if (m_haveSecantBase)
    {
      for (::Vec& v : column)
      {
        ierr = VecDuplicate(gOut, &v);
        assert(ierr == PETSC_SUCCESS);
      }
    }
    // r^k into the input vectors, after forming the new differences.
    ierr = VecAYPX(m_inputData, -1, gOut);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecAYPX(m_inputTraction, -1, tOut);
    assert(ierr == PETSC_SUCCESS);
    if (m_haveSecantBase)
    {
      ierr = VecWAXPY(column[0], -1, m_residualData, m_inputData);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecWAXPY(column[1], -1, m_residualTraction, m_inputTraction);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecWAXPY(column[2], -1, m_outputData, gOut);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecWAXPY(column[3], -1, m_outputTraction, tOut);
      assert(ierr == PETSC_SUCCESS);
      m_secants.insert(m_secants.begin(), column);
      m_secantSteps.insert(m_secantSteps.begin(), m_stepIndex);
      while (m_secants.size() > m_parameters.quasiNewtonColumns)
      {
        for (::Vec& v : m_secants.back())
        {
          ierr = VecDestroy(&v);
          assert(ierr == PETSC_SUCCESS);
        }
        m_secants.pop_back();
        m_secantSteps.pop_back();
      }
    }
    ierr = VecCopy(m_inputData, m_residualData);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecCopy(m_inputTraction, m_residualTraction);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecCopy(gOut, m_outputData);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecCopy(tOut, m_outputTraction);
    assert(ierr == PETSC_SUCCESS);
    m_haveSecantBase = true;
    const size_t n = m_secants.size();
    if (n == 0)
      return;
    Eigen::MatrixXd gram(n, n);
    Eigen::VectorXd rhs(n);
    for (size_t i = 0; i < n; ++i)
    {
      PetscScalar a = 0, b = 0;
      ierr = VecDot(m_secants[i][0], m_residualData, &a);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecDot(m_secants[i][1], m_residualTraction, &b);
      assert(ierr == PETSC_SUCCESS);
      rhs(i) = -PetscRealPart(a + b);
      for (size_t j = 0; j <= i; ++j)
      {
        ierr = VecDot(m_secants[i][0], m_secants[j][0], &a);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecDot(m_secants[i][1], m_secants[j][1], &b);
        assert(ierr == PETSC_SUCCESS);
        gram(i, j) = gram(j, i) = PetscRealPart(a + b);
      }
    }
    // Rank-revealing solve of the normal equations; nearly dependent
    // secants (relative pivot below the threshold) are discarded.
    constexpr Real DependencyThreshold = 1e-12;
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(gram);
    qr.setThreshold(DependencyThreshold);
    const Eigen::VectorXd coefficients = qr.solve(rhs);
    for (size_t i = 0; i < n; ++i)
    {
      if (!std::isfinite(coefficients(i)))
        continue;
      ierr = VecAXPY(gOut, coefficients(i), m_secants[i][2]);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecAXPY(tOut, coefficients(i), m_secants[i][3]);
      assert(ierr == PETSC_SUCCESS);
    }
    (void)ierr;
    // Refresh ghost values of the two modified fields.
    m_solid.m_robinData *= Real(1);
    m_trialTraction *= Real(1);
  }

  FSI::~FSI()
  {
    clearHistory(true);
    for (::Vec* v : {&m_inputData, &m_inputTraction, &m_outputData, &m_outputTraction,
           &m_residualData, &m_residualTraction})
    {
      if (*v)
      {
        const PetscErrorCode ierr = VecDestroy(v);
        assert(ierr == PETSC_SUCCESS);
        (void)ierr;
      }
    }
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
    // the interface convective correction. All three traces are P1. At a
    // fixed point u_f = v_Gamma exactly, unlike recovered-gradient tractions.
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
    m_previousDisplacementIteration.setData(m_solid.m_displacement.getData());
    ++m_stepIndex;
    clearHistory();
    for (size_t iteration = 1; iteration <= m_parameters.maximumIterations; ++iteration)
    {
      Stage coupling("Coupling iteration " + std::to_string(iteration), m_fluid.isRoot());
      if (m_fluid.isRoot())
      {
        Log("FSI/Coupling") << "iteration=" << iteration << '/'
                            << m_parameters.maximumIterations << Alert::Raise;
      }
      m_previousRobinData.setData(m_solid.m_robinData.getData());
      if (m_parameters.quasiNewton)
      {
        PetscErrorCode ierr = VecCopy(m_solid.m_robinData.sync().getData(), m_inputData);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecCopy(m_importedTraction.sync().getData(), m_inputTraction);
        assert(ierr == PETSC_SUCCESS);
        (void)ierr;
      }
      m_solid.solveStep();
      m_residualDifference = m_solid.m_displacement - m_previousDisplacementIteration;
      const Real change = m_residualDifference.norm(NORM_2);
      // Normalized by the motion of this step, not by the accumulated
      // displacement: |d| includes all previous motion and prestretch, which
      // makes a state-relative test increasingly lax as the cycle proceeds.
      m_tractionResidual = m_solid.m_displacement - m_solid.m_previousDisplacement;
      const Real increment = m_tractionResidual.norm(NORM_2);
      m_diagnostics.displacementChange = increment > 0 ? change / increment
        : change > 0 ? std::numeric_limits<Real>::infinity()
                     : 0;
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
      m_diagnostics.converged = iteration > 1 &&
        m_diagnostics.displacementChange < m_parameters.relativeTolerance;
      // A converged sweep commits its own (unaccelerated) fluid reaction.
      if (m_parameters.quasiNewton && !m_diagnostics.converged &&
        iteration < m_parameters.maximumIterations)
        accelerate();
      m_importedTraction.setData(m_trialTraction.getData());
      transfer(m_solidToFluid, m_importedTraction, m_fluid.m_wallTraction, 1e3);
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
    m_fluid.recoverGradient();
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
