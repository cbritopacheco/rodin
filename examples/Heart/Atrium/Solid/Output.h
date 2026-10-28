/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ATRIUM_SOLID_OUTPUT_H
#define RODIN_ATRIUM_SOLID_OUTPUT_H

#include <fstream>
#include <Rodin/IO/XDMF.h>
#include "ProblemData.h"

namespace Rodin::Examples::Heart::AtriumContraction
{
  /**
   * @brief XDMF export and diagnostics for the active solid state.
   *
   * ## Architecture
   *
   * The writer observes the mechanical fields and owns the scalar activation
   * fields used for visualization. Each write validates the deformed wedges,
   * exports the state, and updates cached scalar diagnostics. Standalone CSV
   * output is optional; coordinated runs use the common CSV. Metadata describes the
   * prescribed data and integration parameters. Observed inputs must outlive
   * this object; close() finalizes the XDMF document.
   */
  class Output
  {
    public:
      /**
       * @brief Binds the mechanical state and configures output.
       * @param cfg Output and model configuration.
       * @param data Prescribed mesh and fields.
       * @param diagnostics Destination for cached scalar diagnostics.
       * @param displacement Current displacement.
       * @param velocity Current velocity.
       * @param acceleration Current acceleration.
       * @param extension Solved active extension.
       * @param gammaField Square root of contractile stiffness.
       * @param betaField Contractile stress divided by square root of stiffness.
       * @param timeStep Time step in seconds.
       * @param density Density in kPa s^2/mm^2.
       * @param beta Newmark displacement parameter.
       * @param gamma Newmark velocity parameter.
       */
      Output(const Options& cfg, const ProblemData& data,
        Optional<DiagnosticState>& diagnostics, const VectorField& displacement,
        const VectorField& velocity, const VectorField& acceleration,
        const CellScalarField& extension, const CellScalarField& gammaField,
        const CellScalarField& betaField, Real timeStep, Real density, Real beta,
        Real gamma);

      /**
       * @brief Writes a mechanical state and its diagnostics.
       * @param time Physical time in seconds.
       * @param iterations Newton iteration count.
       * @param residual Final assembled residual norm.
       */
      void write(Real time, size_t iterations, Real residual);

      /// @brief Finalizes the XDMF document.
      void close();

    private:
      const Options& m_cfg;
      const ProblemData& m_data;
      Optional<DiagnosticState>& m_diagnostics;
      const VectorField& m_displacement;
      const VectorField& m_velocity;
      const VectorField& m_acceleration;
      const CellScalarField& m_extension;
      const CellScalarField& m_gamma;
      const CellScalarField& m_beta;
      ScalarField m_activation, m_electrical;
      CellScalarField m_contractileStress, m_contractileStiffness;
      std::ofstream m_log;
      IO::XDMF m_output;
  };
}

#endif
