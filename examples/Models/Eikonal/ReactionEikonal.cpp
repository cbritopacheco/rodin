/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Courtemanche RE+ driven by anisotropic FIM, with tensor FMM comparison.
 *
 * Implements the monodomain structure of Nagel et al. (2023), equation (12),
 * on a rotated, stretched tissue sheet. See ReactionEikonal.md for units, parameter
 * provenance, velocity calibration, time discretization and comparison limits.
 */
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <Eigen/SparseCholesky>
#include <Rodin/Assembly.h>
#include <Rodin/Eikonal/FIM.h>
#include <Rodin/Eikonal/FMM.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/Variational.h>
#include "Courtemanche1998.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Examples
{
  /**
   * @brief First-order split, mass-lumped P1 reaction-diffusion simulation.
   *
   * The normalized equation is @f$\partial_t V=\nabla\cdot(D\nabla V)
   * -I_{ion}/C_m+I_s/C_m@f$, with natural zero-flux boundary conditions.
   * Ionic states are explicit nodal grid functions. Rodin assembles the mass
   * and tensor diffusion forms; the fixed backward-Euler diffusion matrix is
   * factored once. A diffusion step preserves constants and the lumped integral.
   */
  class ReactionDiffusion
  {
    public:
      using Space = P1<Real, LocalMesh>;
      using Field = GridFunction<Space, Math::Vector<Real>>;

      /**
       * @brief Constructs the fixed diffusion solve and all ionic fields.
       * @param mesh Borrowed native mesh, fixed throughout this object's lifetime.
       * @param diffusion Normalized tensor in square millimetres per millisecond.
       * @param dt Time step in milliseconds.
       * @param initial Uniform initial ionic state.
       */
      ReactionDiffusion(const LocalMesh& mesh, const Math::SpatialMatrix<Real>& diffusion,
        Real dt, const Courtemanche1998::State& initial)
        : m_mesh(mesh),
          m_space(mesh),
          m_dt(dt),
          m_activation(m_space),
          m_maximumSlope(m_space),
          m_peakVoltage(m_space)
      {
        for (Real value : initial)
        {
          m_states.emplace_back(m_space);
          m_states.back() = RealFunction(value);
        }
        m_activation = RealFunction(std::numeric_limits<Real>::infinity());
        m_maximumSlope = RealFunction(0);
        m_peakVoltage = RealFunction(initial[0]);
        TrialFunction u(m_space);
        TestFunction v(m_space);
        BilinearForm mass(u, v), stiffness(u, v);
        mass = Integral(u, v);
        stiffness = Integral(MatrixFunction(diffusion) * Grad(u), Grad(v));
        mass.assemble();
        stiffness.assemble();
        m_mass = mass.getOperator() * Math::Vector<Real>::Ones(m_space.getSize());
        if (!(m_mass.array() > 0).all())
          throw std::runtime_error("Nonpositive lumped mass");
        Math::SparseMatrix<Real> matrix = dt * stiffness.getOperator();
        for (Index i = 0; i < static_cast<size_t>(m_mass.size()); ++i)
          matrix.coeffRef(i, i) += m_mass[i];
        m_diffusion.compute(matrix);
        if (m_diffusion.info() != Eigen::Success)
          throw std::runtime_error("Diffusion factorization failed");
      }

      /**
       * @brief Runs a single paced activation and records steepest-upstroke LATs.
       * @param endTime Requested end time in milliseconds, rounded up to a step.
       * @param stimulated Stimulated global vertex indices.
       * @param prefix Optional XDMF time-series basename; empty suppresses output.
       * @param schedule Optional nodal foot-current onset times; null runs monodomain.
       */
      void run(Real endTime, const IndexVector& stimulated,
        const std::string& prefix = "", const Field* schedule = nullptr)
      {
        std::vector<bool> stimulus(m_mesh.getVertexCount(), false);
        for (Index vertex : stimulated)
          stimulus[vertex] = true;
        std::vector<bool> footFinished(m_mesh.getVertexCount(), false);
        std::unique_ptr<IO::XDMF> output;
        if (!prefix.empty())
        {
          output = std::make_unique<IO::XDMF>(prefix);
          output->grid()
            .setMesh(m_mesh)
            .add("voltage_mV", m_states[0])
            .add("calcium_mM", m_states[Courtemanche1998::Ca_i]);
          output->write(0);
        }
        const size_t steps = static_cast<size_t>(std::ceil(endTime / m_dt));
        const size_t outputStride = std::max(size_t(1), size_t(std::llround(1 / m_dt)));
        for (size_t step = 0; step < steps; ++step)
        {
          const Real time = step * m_dt;
          const auto previous = m_states[0].getData();
          // Average the pulse over this step, including its final partial overlap.
          const Real pulse = 20 * std::max(Real(0), std::min(m_dt, 2 - time)) / m_dt;
          for (Index vertex = 0; vertex < m_mesh.getVertexCount(); ++vertex)
          {
            const Index dof = m_space.getDOFs(0, vertex)[0];
            Courtemanche1998::State state{};
            for (size_t component = 0; component < state.size(); ++component)
              state[component] = m_states[component][dof];
            Real current = stimulus[vertex] ? pulse : 0;
            if (schedule)
            {
              // Positive depolarizing convention; openCARP uses negative I_diff.
              constexpr Real amplitude = Real(2.464218456421841); // mV
              constexpr Real tau = Real(0.43917740479584627); // ms
              footFinished[vertex] = footFinished[vertex] || state[0] >= -30;
              current = 0;
              const Real onset = (*schedule)[dof];
              if (!footFinished[vertex] && time + m_dt > onset)
              {
                const Real start = std::max(time, onset);
                current = amplitude / m_dt *
                  (std::exp((time + m_dt - onset) / tau) -
                    std::exp((start - onset) / tau));
              }
            }
            state = Courtemanche1998::step(state, m_dt, current);
            for (size_t component = 0; component < state.size(); ++component)
            {
              if (!std::isfinite(state[component]) ||
                (component >= Courtemanche1998::Na_i && !(state[component] > 0)) ||
                (component > 0 && component < Courtemanche1998::Na_i &&
                  !(state[component] >= 0 && state[component] <= 1)))
                throw std::runtime_error("Invalid ionic state; refine the time step");
              m_states[component][dof] = state[component];
            }
          }
          const Math::Vector<Real> rhs = m_mass.cwiseProduct(m_states[0].getData());
          m_states[0].getData() = m_diffusion.solve(rhs);
          if (m_diffusion.info() != Eigen::Success || !m_states[0].getData().allFinite())
            throw std::runtime_error("Diffusion solve failed");
          for (Index dof = 0; dof < static_cast<size_t>(m_mass.size()); ++dof)
          {
            const Real slope = (m_states[0][dof] - previous[dof]) / m_dt;
            m_peakVoltage[dof] = std::max(m_peakVoltage[dof], m_states[0][dof]);
            if (slope > m_maximumSlope[dof])
            {
              m_maximumSlope[dof] = slope;
              m_activation[dof] = time + m_dt / 2;
            }
          }
          if (output && ((step + 1) % outputStride == 0 || step + 1 == steps))
            output->write(time + m_dt);
        }
        for (Index dof = 0; dof < static_cast<size_t>(m_mass.size()); ++dof)
        {
          if (m_peakVoltage[dof] <= -40)
            m_activation[dof] = std::numeric_limits<Real>::infinity();
        }
        if (output)
          output->close();
      }

      /// @returns Nodal steepest-upstroke times; unreached nodes retain infinity.
      const Field& getActivation() const
      {
        return m_activation;
      }

      /// @returns Simulation space, with its vertex-to-DOF map.
      const Space& getSpace() const
      {
        return m_space;
      }

    private:
      const LocalMesh& m_mesh;
      Space m_space;
      Real m_dt;
      std::vector<Field> m_states;
      Field m_activation, m_maximumSlope, m_peakVoltage;
      Math::Vector<Real> m_mass;
      Eigen::SimplicialLDLT<Math::SparseMatrix<Real>> m_diffusion;
  };
}

int main(int argc, char** argv)
{
  using namespace Rodin::Examples;
  try
  {
    if (argc > 7)
    {
      throw std::runtime_error("Usage: Eikonal_ReactionEikonal [intervals dt_ms end_ms "
                               "pace_beats velocity_ratio workers]");
    }
    const size_t intervals = argc > 1 ? std::stoul(argv[1]) : 40;
    const Real dt = argc > 2 ? std::stod(argv[2]) : Real(0.02);
    const Real endTime = argc > 3 ? std::stod(argv[3]) : Real(80);
    const size_t paceBeats = argc > 4 ? std::stoul(argv[4]) : 10;
    const Real ratio = argc > 5 ? std::stod(argv[5]) : Real(2);
    const size_t workers = argc > 6 ? std::stoul(argv[6]) : 1;
    if (intervals < 20 || !(dt > 0 && dt <= Real(0.02)) ||
      !(endTime >= 40 && endTime <= 1000) || !(ratio >= 1 && ratio <= 4) ||
      !std::isfinite(dt) || !std::isfinite(endTime) || !std::isfinite(ratio) ||
      workers == 0)
    {
      throw std::runtime_error("Require intervals >= 20, 0 < dt <= 0.02 ms, 40 <= end <= "
                               "1000 ms, 1 <= ratio <= 4, workers >= 1");
    }
    constexpr Real length = 10; // Unstretched sheet and strand length in mm.
    // Illustrative tissue values: sigma_t=0.1 S/m, beta=140000 1/m,
    // membrane capacitance=0.01 F/m^2. Convert m^2/s to mm^2/ms by 1000.
    constexpr Real diffusionT = 1000 * Real(0.1) / (140000 * Real(0.01));
    auto initial = Courtemanche1998::initialState();
    std::ofstream cellTrace("ReactionEikonalCell.csv");
    cellTrace << std::setprecision(17) << "time_ms,voltage_mV,calcium_mM\n";
    const size_t cycleSteps = static_cast<size_t>(std::ceil(1000 / dt));
    const Real paceDt = Real(1000) / cycleSteps;
    for (size_t beat = 0; beat < paceBeats; ++beat)
    {
      for (size_t step = 0; step < cycleSteps; ++step)
      {
        const Real time = step * paceDt;
        const Real pulse = 20 * std::max(Real(0), std::min(paceDt, 2 - time)) / paceDt;
        initial = Courtemanche1998::step(initial, paceDt, pulse);
        for (size_t component = 0; component < initial.size(); ++component)
        {
          if (!std::isfinite(initial[component]) ||
            (component >= Courtemanche1998::Na_i && !(initial[component] > 0)))
            throw std::runtime_error("Single-cell conditioning failed");
        }
        if (beat + 1 == paceBeats && step % std::max(size_t(1), cycleSteps / 1000) == 0)
        {
          cellTrace << time + paceDt << ',' << initial[0] << ','
                    << initial[Courtemanche1998::Ca_i] << '\n';
        }
      }
    }
    std::cout << std::setprecision(10) << "conditioning_beats=" << paceBeats
              << " initial_voltage_mV=" << initial[0] << '\n';

    // Calibrate both principal velocities from independent planar strand runs.
    std::array<Real, 2> velocity{};
    for (size_t axis = 0; axis < 2; ++axis)
    {
      const Real stretch = axis == 0 ? ratio : Real(1);
      auto strand = LocalMesh::UniformGrid(Polytope::Type::Segment, {intervals + 1});
      strand.scale(length * stretch / intervals);
      strand.getConnectivity().compute(1, 0);
      strand.getConnectivity().compute(0, 1);
      IndexVector seeds;
      Index first = 0, second = 0;
      for (auto vertex = strand.getVertex(); vertex; ++vertex)
      {
        const Real coordinate = vertex->getCoordinates()(0) / stretch;
        if (coordinate <= Real(0.5) + Real(1e-12))
          seeds.push_back(vertex->getIndex());
        if (std::abs(coordinate - 3) <
          std::abs(strand.getVertexCoordinates(first)(0) / stretch - 3))
          first = vertex->getIndex();
        if (std::abs(coordinate - 7) <
          std::abs(strand.getVertexCoordinates(second)(0) / stretch - 7))
          second = vertex->getIndex();
      }
      Math::SpatialMatrix<Real> diffusion(1, 1);
      diffusion(0, 0) = diffusionT * stretch * stretch;
      ReactionDiffusion calibration(strand, diffusion, dt, initial);
      calibration.run(endTime, seeds);
      const auto& activation = calibration.getActivation();
      const auto& space = calibration.getSpace();
      const Real difference =
        activation[space.getDOFs(0, second)[0]] - activation[space.getDOFs(0, first)[0]];
      velocity[axis] =
        (strand.getVertexCoordinates(second)(0) - strand.getVertexCoordinates(first)(0)) /
        difference;
      if (!(difference > 0) || !std::isfinite(velocity[axis]))
        throw std::runtime_error("Strand did not activate both velocity probes");
    }
    std::cout << "longitudinal_velocity_mm_per_ms=" << velocity[0]
              << " transverse_velocity_mm_per_ms=" << velocity[1] << '\n';

    // Rotate the fiber axes, stretching the physical mesh consistently with D.
    auto mesh =
      LocalMesh::UniformGrid(Polytope::Type::Triangle, {intervals + 1, intervals + 1});
    mesh.scale(length / intervals);
    Math::SpatialMatrix<Real> rotation(2, 2);
    rotation(0, 0) = Real(0.8);
    rotation(0, 1) = Real(-0.6);
    rotation(1, 0) = Real(0.6);
    rotation(1, 1) = Real(0.8);
    Math::SpatialMatrix<Real> stretch = Math::SpatialMatrix<Real>::Identity(2, 2);
    stretch(0, 0) = ratio;
    const Math::SpatialMatrix<Real> transform = rotation * stretch;
    Math::SpatialMatrix<Real> squaredVelocity(2, 2);
    squaredVelocity.setZero();
    squaredVelocity(0, 0) = velocity[0] * velocity[0];
    squaredVelocity(1, 1) = velocity[1] * velocity[1];
    const Math::SpatialMatrix<Real> tensor =
      rotation * squaredVelocity * rotation.transpose();
    const Math::SpatialMatrix<Real> diffusion =
      diffusionT * transform * transform.transpose();
    IndexVector seeds;
    std::vector<bool> source(mesh.getVertexCount(), false);
    for (auto vertex = mesh.getVertex(); vertex; ++vertex)
    {
      Math::SpatialPoint reference = vertex->getCoordinates();
      // Equilateral reference triangles provide a strict causal angle margin.
      reference(0) += Real(0.5) * reference(1);
      reference(1) *= std::sqrt(Real(3)) / 2;
      if (reference.norm() <= Real(0.5) + Real(1e-12))
      {
        seeds.push_back(vertex->getIndex());
        source[vertex->getIndex()] = true;
      }
      mesh.setVertexCoordinates(vertex->getIndex(), transform * reference);
    }
    mesh.getConnectivity().compute(2, 0);
    mesh.getConnectivity().compute(0, 2);
    mesh.getConnectivity().compute(0, 0);
    ReactionDiffusion tissue(mesh, diffusion, dt, initial);
    const auto& space = tissue.getSpace();
    GridFunction iterative(space), marching(space), difference(space);
    MatrixFunction coefficient(tensor);
    Eikonal::FIM solver(iterative, coefficient);
    solver.setThreadCount(workers).seed(seeds).solve();
    Eikonal::FMM(marching, coefficient).seed(seeds).solve();
    tissue.run(endTime, seeds, "ReactionEikonalVoltage", &iterative);
    std::ofstream comparison("ReactionEikonalActivation.csv");
    comparison << std::setprecision(17)
               << "x_mm,y_mm,stimulated,replus_lat_ms,fim_foot_onset_"
                  "ms,fmm_foot_onset_ms,lat_minus_foot_onset_ms\n";
    Real maximumDifference = 0, sum = 0, maximumMarchingDifference = 0;
    size_t count = 0;
    for (auto vertex = mesh.getVertex(); vertex; ++vertex)
    {
      const Index index = vertex->getIndex();
      const Index dof = space.getDOFs(0, index)[0];
      const Real activation = tissue.getActivation()[dof];
      difference[dof] = activation - iterative[dof];
      if (!std::isfinite(activation))
        throw std::runtime_error("Unactivated RE+ tissue; increase end_ms");
      if (!source[index])
      {
        maximumDifference = std::max(maximumDifference, std::abs(difference[dof]));
        sum += difference[dof] * difference[dof];
        ++count;
      }
      maximumMarchingDifference =
        std::max(maximumMarchingDifference, std::abs(iterative[dof] - marching[dof]));
      comparison << vertex->getCoordinates()(0) << ',' << vertex->getCoordinates()(1)
                 << ',' << source[index] << ',' << activation << ',' << iterative[dof]
                 << ',' << marching[dof] << ',' << difference[dof] << '\n';
    }
    IO::XDMF output("ReactionEikonalActivation");
    output.grid()
      .setMesh(mesh)
      .add("replus_lat_ms", tissue.getActivation())
      .add("fim_foot_onset_ms", iterative)
      .add("fmm_foot_onset_ms", marching)
      .add("lat_minus_foot_onset_ms", difference);
    output.write().close();
    std::cout << "fim_rounds=" << solver.getIterations()
              << " max_fim_fmm_ms=" << maximumMarchingDifference
              << " max_abs_lat_minus_onset_ms=" << maximumDifference
              << " rms_lat_minus_onset_ms=" << std::sqrt(sum / count) << '\n';
    if (maximumMarchingDifference > 1e-8)
      throw std::runtime_error("FIM and FMM disagree on the matched causal mesh");
    return cellTrace && comparison ? EXIT_SUCCESS : EXIT_FAILURE;
  }
  catch (const std::exception& error)
  {
    std::cerr << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
