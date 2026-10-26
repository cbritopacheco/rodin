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
#include <iostream>
#include <Rodin/Alert.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Examples::Heart::AtriumContraction
{
  Output::Output(const Options& cfg, const ProblemData& data,
    const VectorField& displacement, const VectorField& velocity,
    const VectorField& acceleration, Real timeStep, Real density, Real beta, Real gamma)
    : m_cfg(cfg),
      m_data(data),
      m_displacement(displacement),
      m_velocity(velocity),
      m_acceleration(acceleration),
      m_activation(data.scalarSpace),
      m_extension(data.scalarSpace),
      m_output(cfg.output)
  {
    const auto parent = std::filesystem::path(cfg.output).parent_path();
    if (!parent.empty())
      std::filesystem::create_directories(parent);
    m_activation.setName("Activation");
    m_extension.setName("ActiveExtension");
    std::ofstream metadata(m_cfg.output + ".geometry.json");
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
             << ",\n  \"minimum_wave_coordinate_mm\": " << m_data.minimumProjection
             << ",\n  \"maximum_wave_coordinate_mm\": " << m_data.maximumProjection
             << ",\n  \"speed_mm_per_s\": " << m_cfg.speed
             << ",\n  \"width_mm\": " << m_cfg.width
             << ",\n  \"rayleigh_alpha_per_s\": " << m_cfg.rayleighAlpha
             << ",\n  \"rayleigh_beta_s\": " << m_cfg.rayleighBeta
             << ",\n  \"reference_density_kg_per_m3\": " << m_cfg.density
             << ",\n  \"density_kPa_s2_per_mm2\": " << density
             << ",\n  \"time_step_s\": " << timeStep << ",\n  \"newmark_beta\": " << beta
             << ",\n  \"newmark_gamma\": " << gamma << "\n}\n";
    metadata.close();
    std::cout << "Base label " << m_cfg.baseLabel
              << ": origin=" << m_data.baseOrigin.transpose()
              << "; inward normal=" << m_data.waveDirection.transpose()
              << "; appendage projection="
              << m_data.waveDirection.dot(m_cfg.appendageApex - m_data.baseOrigin)
              << " mm; harmonic residual=" << m_data.coordinateResidual
              << "; fiber fallback cells=" << m_data.fallbackCells << std::endl;
    auto grid = m_output.grid();
    grid.setMesh(m_data.mesh);
    grid.add(m_displacement);
    grid.add(m_velocity);
    grid.add(m_acceleration);
    grid.add(m_activation);
    grid.add(m_extension);
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
      m_data.waveCoordinate, IO::XDMF::Center::Node, IO::XDMF::AttributePolicy::Static);
    if (!cfg.fibersOnly)
    {
      m_log.open(cfg.output + ".csv");
      m_log << "time_s,max_displacement_mm,activation_max,newton_iterations,residual,min_"
               "scaled_jacobian_mm,max_velocity_mm_per_s,max_acceleration_mm_per_s2\n";
      m_log.flush();
    }
    std::cout << std::setprecision(9) << "Wall: " << m_data.mesh.getCellCount()
              << " wedges, " << m_data.mesh.getVertexCount()
              << " vertices; fixed rim faces=" << m_data.rimFaces
              << "; minimum reference scaled Jacobian=" << m_data.minimumReferenceJacobian
              << "; dt=" << timeStep << "; run end=" << timeStep * cfg.steps
              << "; wave end=" << m_data.endTime << std::endl;
  }

  void Output::write(Real time, size_t iterations, Real residual)
  {
    m_activation = [&](const Point& point) { return m_data.activationAt(point, time); };
    m_extension = [&](const Point& point) {
      return -m_cfg.contraction * m_data.activationAt(point, time);
    };
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
        u(j) = m_displacement.getData()(dofs(j));
        v(j) = m_velocity.getData()(dofs(j));
        a(j) = m_acceleration.getData()(dofs(j));
      }
      maximum = std::max(maximum, u.norm());
      maximumVelocity = std::max(maximumVelocity, v.norm());
      maximumAcceleration = std::max(maximumAcceleration, a.norm());
    }
    Real minimumCurrentJacobian = std::numeric_limits<Real>::infinity();
    for (auto cell = m_data.mesh.getCell(); cell; ++cell)
    {
      std::array<Vec, 6> q;
      for (size_t j = 0; j < 6; ++j)
      {
        const Index vertex = cell->getVertices()(j);
        const auto x = m_data.mesh.getVertexCoordinates(vertex);
        const auto& dofs = m_displacement.getFiniteElementSpace().getDOFs(0, vertex);
        for (size_t k = 0; k < 3; ++k)
          q[j](k) = x(k) + m_displacement.getData()(dofs(k));
      }
      const Real determinant = minimumJacobian(q);
      if (!(determinant > 0) || !std::isfinite(determinant))
      {
        Alert::Exception() << "The contraction produced a folded wedge" << Alert::Raise;
      }
      const Real area = (q[1] - q[0]).cross(q[2] - q[0]).norm();
      minimumCurrentJacobian = std::min(minimumCurrentJacobian, determinant / area);
    }
    m_output.write(time).flush();
    m_log << std::setprecision(17) << time << ',' << maximum << ','
          << m_activation.getData().maxCoeff() << ',' << iterations << ',' << residual
          << ',' << minimumCurrentJacobian << ',' << maximumVelocity << ','
          << maximumAcceleration << '\n';
    m_log.flush();
    std::cout << "  max displacement=" << maximum << " mm" << std::endl;
  }

  void Output::close()
  {
    m_output.close();
  }
}
