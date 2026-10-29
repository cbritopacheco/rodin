/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Output.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include "../Log.h"
#include <Rodin/Alert.h>
#include <boost/mpi/collectives.hpp>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Examples::Heart::Atrium
{
  Solid::Output::Output(const SolidParameters& cfg, const RunParameters& run,
    const ProblemData& data, Optional<DiagnosticState>& diagnostics,
    const VectorField& displacement, const VectorField& velocity,
    const VectorField& acceleration, const CellScalarField& extension,
    const CellScalarField& gammaField, const CellScalarField& betaField, Real timeStep,
    Real density, Real beta, Real gamma)
    : m_cfg(cfg),
      m_run(run),
      m_data(data),
      m_diagnostics(diagnostics),
      m_displacement(displacement),
      m_velocity(velocity),
      m_acceleration(acceleration),
      m_extension(extension),
      m_gamma(gammaField),
      m_beta(betaField),
      m_activation(data.scalarSpace),
      m_electrical(data.scalarSpace),
      m_contractileStress(data.diagnosticSpace),
      m_contractileStiffness(data.diagnosticSpace),
      m_output(
        data.mesh.getContext().getCommunicator(), (run.outputDirectory + "/Solid"), 0)
  {
    const auto parent =
      std::filesystem::path((run.outputDirectory + "/Solid")).parent_path();
    if (!parent.empty())
      std::filesystem::create_directories(parent);
    m_activation.setName("Activation");
    m_electrical.setName("ElectricalActivation");
    m_contractileStress.setName("ContractileStress");
    m_contractileStiffness.setName("ContractileStiffness");
    const bool isRoot = data.mesh.getContext().getCommunicator().rank() == 0;
    if (isRoot)
    {
      std::ofstream metadata((m_run.outputDirectory + "/Solid") + ".geometry.json");
    metadata << std::setprecision(17) << "{\n  \"base_label\": " << m_cfg.baseLabel
             << ",\n  \"base_origin_mm\": [" << m_data.baseOrigin(0) << ','
             << m_data.baseOrigin(1) << ',' << m_data.baseOrigin(2)
             << "],\n  \"wave_direction\": [" << m_data.waveDirection(0) << ','
             << m_data.waveDirection(1) << ',' << m_data.waveDirection(2)
             << "],\n  \"appendage_apex_mm\": [" << m_cfg.appendageApex(0) << ','
             << m_cfg.appendageApex(1) << ',' << m_cfg.appendageApex(2)
             << "],\n  \"apex_radius_mm\": " << m_cfg.apexRadius
             << ",\n  \"apex_patch_faces\": " << m_data.apexPatchFaces
             << ",\n  \"fiber_fallback_cells\": " << m_data.fallbackCells
             << ",\n  \"fiber_coordinate_residual\": " << m_data.coordinateResidual
             << ",\n  \"transmural_coordinate_residual\": " << m_data.transmuralResidual
             << ",\n  \"max_endocardial_fiber_normal_component\": " << m_data.endocardialTangencyError
             << ",\n  \"max_epicardial_fiber_normal_component\": " << m_data.epicardialTangencyError
             << ",\n  \"minimum_wave_coordinate_mm\": " << m_data.minimumProjection
             << ",\n  \"maximum_wave_coordinate_mm\": " << m_data.maximumProjection
             << ",\n  \"speed_mm_per_s\": " << m_cfg.speed
             << ",\n  \"width_mm\": " << m_cfg.width
             << ",\n  \"activation_peak_per_s\": " << m_cfg.activation
             << ",\n  \"relaxation_drive_per_s\": " << m_cfg.relaxation
             << ",\n  \"matrix_scale_kPa\": " << m_cfg.passive.matrixScale
             << ",\n  \"matrix_exponent\": " << m_cfg.passive.matrixExponent
             << ",\n  \"fiber_scale_kPa\": " << m_cfg.passive.fiberScale
             << ",\n  \"fiber_exponent\": " << m_cfg.passive.fiberExponent
             << ",\n  \"volume_penalty_kPa\": " << m_cfg.passive.volumePenalty
             << ",\n  \"tissue_viscosity_kPa_s\": " << m_cfg.viscosity
             << ",\n  \"active_stiffness_kPa\": " << m_cfg.activeStiffness
             << ",\n  \"active_damping_kPa_s\": " << m_cfg.activeDamping
             << ",\n  \"cross_bridge_stiffness_kPa\": " << m_cfg.crossBridgeStiffness
             << ",\n  \"contractility_kPa\": " << m_cfg.contractility
             << ",\n  \"destruction_rate\": " << m_cfg.destructionRate
             << ",\n  \"reference_density_kg_per_m3\": " << m_cfg.density
             << ",\n  \"density_kPa_s2_per_mm2\": " << density
             << ",\n  \"time_step_s\": " << timeStep << ",\n  \"newmark_beta\": " << beta
             << ",\n  \"newmark_gamma\": " << gamma << "\n}\n";
    metadata.close();
      Log("Solid/Fibers") << "base_label=" << m_cfg.baseLabel
              << ": origin=" << m_data.baseOrigin.transpose()
              << "; inward normal=" << m_data.waveDirection.transpose()
              << "; appendage projection="
              << m_data.waveDirection.dot(m_cfg.appendageApex - m_data.baseOrigin)
              << " mm; harmonic residual=" << m_data.coordinateResidual
              << "; fiber fallback cells=" << m_data.fallbackCells
              << "; transmural residual=" << m_data.transmuralResidual
              << "; tangency (endo, epi)=" << m_data.endocardialTangencyError
                          << ", " << m_data.epicardialTangencyError << Alert::Raise;
    }
    auto grid = m_output.grid();
    grid.setMesh(m_data.mesh);
    grid.add(m_displacement);
    grid.add(m_velocity);
    grid.add(m_acceleration);
    grid.add(m_activation);
    grid.add(m_electrical);
    grid.add(m_extension, IO::XDMF::Center::Cell);
    grid.add(m_gamma, IO::XDMF::Center::Cell);
    grid.add(m_beta, IO::XDMF::Center::Cell);
    grid.add(m_contractileStress, IO::XDMF::Center::Cell);
    grid.add(m_contractileStiffness, IO::XDMF::Center::Cell);
    grid.add(
      m_data.fiberField, IO::XDMF::Center::Cell, IO::XDMF::AttributePolicy::Static);
    grid.add(
      m_data.normalField, IO::XDMF::Center::Cell, IO::XDMF::AttributePolicy::Static);
    grid.add(m_data.longitudinalField, IO::XDMF::Center::Cell,
      IO::XDMF::AttributePolicy::Static);
    grid.add(
      m_data.fallbackField, IO::XDMF::Center::Cell, IO::XDMF::AttributePolicy::Static);
    grid.add(
      m_data.fiberCoordinate, IO::XDMF::Center::Node, IO::XDMF::AttributePolicy::Static);
    grid.add(
      m_data.transmuralCoordinate, IO::XDMF::Center::Node, IO::XDMF::AttributePolicy::Static);
    grid.add(
      m_data.waveCoordinate, IO::XDMF::Center::Node, IO::XDMF::AttributePolicy::Static);
    const size_t cells = m_data.mesh.getCellCount(),
                 vertices = m_data.mesh.getVertexCount();
    if (isRoot)
    {
      Log("Solid/Mesh") << "wall: " << cells << " wedges, " << vertices
              << " vertices; fixed rim faces=" << m_data.rimFaces
                        << "; minimum reference scaled Jacobian="
                        << m_data.minimumReferenceJacobian << "; dt=" << timeStep
                        << " s; run end=" << timeStep * run.steps << " s"
                        << "; wave end=" << m_data.endTime << " s" << Alert::Raise;
    }
  }

  void Solid::Output::write(Real time, size_t iterations, Real residual)
  {
    m_activation = [&](const Point& point) { return m_data.activationAt(point, time); };
    m_electrical = [&](const Point& point) { return m_data.electricalAt(point, time); };
    m_contractileStress = m_gamma * m_beta;
    m_contractileStiffness = m_gamma * m_gamma;
    if (m_cfg.fibersOnly)
    {
      m_output.write(time).flush();
      return;
    }
    Real maximum = 0, maximumVelocity = 0, maximumAcceleration = 0;
    for (auto vertex = m_data.mesh.getVertex(); vertex; ++vertex)
    {
      const auto& dofs =
        m_displacement.getFiniteElementSpace().getDOFs(0, vertex->getIndex());
      Vec u, v, a;
      for (size_t j = 0; j < 3; ++j)
      {
        u(j) = m_displacement[dofs(j)];
        v(j) = m_velocity[dofs(j)];
        a(j) = m_acceleration[dofs(j)];
      }
      maximum = std::max(maximum, u.norm());
      maximumVelocity = std::max(maximumVelocity, v.norm());
      maximumAcceleration = std::max(maximumAcceleration, a.norm());
    }
    const auto& comm = m_data.mesh.getContext().getCommunicator();
    maximum = boost::mpi::all_reduce(comm, maximum, boost::mpi::maximum<Real>());
    maximumVelocity =
      boost::mpi::all_reduce(comm, maximumVelocity, boost::mpi::maximum<Real>());
    maximumAcceleration =
      boost::mpi::all_reduce(comm, maximumAcceleration, boost::mpi::maximum<Real>());
    const auto minimum = m_data.minimumDeformedJacobian(m_displacement);
    if (!minimum || *minimum <= 0)
      Alert::Exception() << "The contraction produced a folded wedge" << Alert::Raise;
    const Real minimumCurrentJacobian = *minimum;
    m_time = time;
    const size_t step = static_cast<size_t>(std::llround(time / m_run.timeStep));
    m_pendingOutput =
      step != 0 && (m_run.outputEvery == 0 || step % m_run.outputEvery != 0);
    if (!m_pendingOutput)
      m_output.write(time).flush();
    Index index;
    const Real activationMax = m_activation.max(index);
    const Real extensionMin = m_extension.min(index),
               extensionMax = m_extension.max(index);
    const Real stressMax = m_contractileStress.max(index);
    m_diagnostics = DiagnosticState{maximum, activationMax, iterations, residual,
      minimumCurrentJacobian, maximumVelocity, maximumAcceleration, extensionMin,
      extensionMax, stressMax};
    if (comm.rank() == 0)
    {
      Log("Solid/State") << "max_displacement=" << maximum << " mm" << Alert::Raise;
    }
  }

  void Solid::Output::close()
  {
    if (m_pendingOutput)
      m_output.write(m_time).flush();
    m_output.close();
  }
}
