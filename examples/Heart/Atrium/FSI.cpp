/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <cmath>
#include <limits>
#include <boost/mpi/collectives.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include "FSI.h"
#include "CavityVolume.h"

namespace Rodin::Examples::Heart::Atrium
{
  using namespace Geometry;
  using namespace Variational;

  Fluid::MeshType FSI::makeInterface(const Fluid& fluid)
  {
    Fluid::MeshType::Builder builder;
    builder.initialize(fluid.m_mesh);
    for (auto face = fluid.m_mesh.getBoundary(); face; ++face)
    {
      if (face->getAttribute() == fluid.m_cfg.wall &&
        fluid.m_mesh.getShard().isOwned(2, face->getIndex()))
        builder.include(2, face->getIndex());
    }
    auto mesh = builder.finalize();
    mesh.getConnectivity().compute(2, 0);
    return mesh;
  }

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
      m_interface(makeInterface(fluid)),
      m_interfaceSpace(std::integral_constant<size_t, 1>{}, m_interface, 3),
      m_motionTrial(fluid.m_vh),
      m_tractionTrial(m_interfaceSpace),
      m_motionTest(fluid.m_vh),
      m_tractionTest(m_interfaceSpace),
      m_motion(m_motionTrial, m_motionTest),
      m_traction(m_tractionTrial, m_tractionTest),
      m_motionSolver(m_motion),
      m_tractionSolver(m_traction),
      m_previousMeshDisplacement(fluid.m_vh),
      m_importedTraction(solid.m_space),
      m_importedVelocity(solid.m_space),
      m_previousRobinData(solid.m_space),
      m_previousResidual(solid.m_space),
      m_residualDifference(solid.m_space),
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
    // Boundary impedances are value coefficients in the retained flow form.
    fluid.setupFlow();
    if (m_robinCoefficient == 0)
    {
      m_robinCoefficient =
        fluid.m_cfg.rho * std::cbrt(CavityVolume::measure(fluid.m_mesh)) / fluid.m_cfg.dt;
    }
    // Pa s/m -> kPa s/mm. Robin data are expressed per reference solid area.
    solid.m_robinCoefficient = 1e-6 * m_robinCoefficient;
    const auto velocity = solid.m_predictedVelocity +
      solid.m_velocityCoefficient *
        (solid.m_displacement - solid.m_predictedDisplacement);
    solid.m_mechanics += solid.m_robinCoefficient * solid.m_velocityCoefficient *
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
    m_tractionSolver.setPrefix(std::string("atrium_traction_")).setType(KSPCG);
    m_tractionSolver.setTolerances(1e-12, 1e-30, 1e5, 10000);

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
    m_previousMeshDisplacement = zero;
    m_importedTraction = zero;
    m_importedVelocity = zero;
    m_previousRobinData = zero;
    m_previousResidual = zero;
    m_residualDifference = zero;
    m_solidToFluid = makeScatter(
      solid.m_mesh, solid.m_displacement, fluid.m_mesh, fluid.m_wallDisplacement);
    m_tractionToSolid = makeScatter(
      m_interface, m_tractionTrial.getSolution(), solid.m_mesh, m_importedTraction);
  }

  FSI::~FSI()
  {
    if (m_solidToFluid)
    {
      const PetscErrorCode ierr = VecScatterDestroy(&m_solidToFluid);
      assert(ierr == PETSC_SUCCESS);
      (void)ierr;
    }
    if (m_tractionToSolid)
    {
      const PetscErrorCode ierr = VecScatterDestroy(&m_tractionToSolid);
      assert(ierr == PETSC_SUCCESS);
      (void)ierr;
    }
  }

  void FSI::initialize()
  {
    moveFluidMesh();
    m_previousMeshDisplacement.setData(m_fluid.m_meshDisplacement.getData());
    m_fluid.m_meshVelocity = Math::SpatialVector<Real>{0, 0, 0};
  }

  void FSI::restoreReference()
  {
    for (Index i = 0; i < m_referenceVertices.size(); ++i)
      m_fluid.m_mesh.setVertexCoordinates(i, m_referenceVertices[i]);
  }

  void FSI::moveFluidMesh()
  {
    transfer(m_solidToFluid, m_solid.m_displacement, m_fluid.m_wallDisplacement, 1e-3);
    transfer(m_solidToFluid, m_solid.m_velocity, m_fluid.m_wallVelocity, 1e-3);
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
      (m_fluid.m_meshDisplacement - m_previousMeshDisplacement) / m_fluid.m_cfg.dt;
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
    const auto normal = BoundaryNormal(m_fluid.m_mesh);
    VectorFunction traction(3, [this, normal](const Point& reference) {
      const Index faceIndex =
        m_interface.getPolytopeMap(2).left.at(reference.getPolytope().getIndex());
      const auto face = m_fluid.m_mesh.getFace(faceIndex);
      const Point current(*face, reference.getReferenceCoordinates());
      const auto n = normal.getValue(current);
      Math::SpatialMatrix<Real> gradient(3, 3);
      const std::array rows{m_fluid.m_gradRec0(current), m_fluid.m_gradRec1(current),
        m_fluid.m_gradRec2(current)};
      for (Index i = 0; i < 3; ++i)
      {
        for (Index j = 0; j < 3; ++j)
          gradient(i, j) = rows[i](j);
      }
      const auto& cells =
        m_fluid.m_mesh.getConnectivity().getIncidence(2, 3).at(faceIndex);
      assert(cells.size() == 1);
      const auto cell = m_fluid.m_mesh.getCell(cells.front());
      const Point center(*cell, Math::SpatialPoint{.25, .25, .25});
      const Real mu = m_fluid.viscosityAt(center);
      const Real stretch = face->getMeasure() / reference.getPolytope().getMeasure();
      return Math::SpatialVector<Real>(stretch *
        (m_fluid.m_p.getSolution()(current) * n -
          mu * (gradient + gradient.transpose()) * n));
    });
    // The interface itself is a two-dimensional mesh: integrate its cells.
    // Opaque coefficient: traction is P1 on each affine interface triangle.
    constexpr size_t TractionOrder = 2;
    auto load = Integral(traction, m_tractionTest);
    load.setOrder(TractionOrder);
    m_traction = Integral(m_tractionTrial, m_tractionTest) - load;
    m_traction.assemble();
    m_traction.solve(m_tractionSolver);
    ::KSPConvergedReason reason;
    const PetscErrorCode ierr =
      KSPGetConvergedReason(m_tractionSolver.getHandle(), &reason);
    assert(ierr == PETSC_SUCCESS);
    (void)ierr;
    if (reason <= 0)
    {
      Alert::Exception() << "Interface traction projection did not converge"
                         << Alert::Raise;
    }
    transfer(m_tractionToSolid, m_tractionTrial.getSolution(), m_importedTraction, 1e-3);
    transfer(m_solidToFluid, m_fluid.m_u.getSolution(), m_importedVelocity, 1e3,
      SCATTER_REVERSE);
    m_solid.m_robinData =
      m_importedTraction + m_solid.m_robinCoefficient * m_importedVelocity;
  }

  void FSI::step()
  {
    m_diagnostics = Diagnostics{};
    m_previousMeshDisplacement.setData(m_fluid.m_meshDisplacement.getData());
    for (auto cell = m_fluid.m_mesh.getCell(); cell; ++cell)
      m_previousCellMeasures[cell->getIndex()] = cell->getMeasure();
    m_solid.beginStep();
    m_fluid.beginStep();
    bool converged = false;
    Real relaxation = m_parameters.initialRelaxation;
    for (size_t iteration = 1; iteration <= m_parameters.maximumIterations; ++iteration)
    {
      m_previousRobinData.setData(m_solid.m_robinData.getData());
      m_solid.solveStep();
      moveFluidMesh();
      if (!m_fluid.solveFlow())
        Alert::Exception() << "FSI fluid velocity exceeded its limit" << Alert::Raise;
      m_fluid.computeWallShear();
      updateTraction();
      PetscReal difference = 0, norm = 0;
      m_previousRobinData.axpy(-1, m_solid.m_robinData);
      PetscErrorCode ierr =
        VecNorm(m_previousRobinData.getData(), NORM_INFINITY, &difference);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecNorm(m_solid.m_robinData.getData(), NORM_INFINITY, &norm);
      assert(ierr == PETSC_SUCCESS);
      (void)ierr;
      m_diagnostics.iterations = iteration;
      m_diagnostics.robinResidual = 1e3 * difference;
      m_diagnostics.robinTolerance =
        m_parameters.absoluteTolerance + m_parameters.relativeTolerance * 1e3 * norm;
      if (m_fluid.isRoot())
      {
        Alert::Info() << "FSI iteration " << iteration
                      << ": Robin residual=" << m_diagnostics.robinResidual << " Pa"
                      << Alert::Raise;
      }
      if (m_diagnostics.robinResidual <= m_diagnostics.robinTolerance)
      {
        converged = true;
        break;
      }
      if (iteration > 1)
      {
        m_residualDifference.setData(m_previousRobinData.getData());
        m_residualDifference.axpy(-1, m_previousResidual);
        PetscScalar numerator = 0, denominator = 0;
        ierr = VecDot(
          m_previousResidual.getData(), m_residualDifference.getData(), &numerator);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecDot(
          m_residualDifference.getData(), m_residualDifference.getData(), &denominator);
        assert(ierr == PETSC_SUCCESS);
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
      m_solid.m_robinData.axpy(1 - relaxation, m_previousRobinData);
      m_diagnostics.relaxation = relaxation;
    }
    if (!converged)
      Alert::Exception() << "Atrium interface iteration did not converge" << Alert::Raise;
    Real velocityMismatch = 0, displacementMismatch = 0;
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
    m_diagnostics.displacementMismatch =
      boost::mpi::all_reduce(comm, displacementMismatch, boost::mpi::maximum<Real>());
    const auto normal = BoundaryNormal(m_fluid.m_mesh);
    m_fluid.m_flux =
      BoundaryIntegral(Dot(m_fluid.m_u.getSolution(), normal), m_fluid.m_qFlux)
        .over(m_fluid.m_cfg.wall);
    m_fluid.m_flux.assemble();
    m_diagnostics.wallFlux = m_fluid.m_flux(m_fluid.m_one);
    m_solid.commitStep();
    m_fluid.commitStep();
  }
}
