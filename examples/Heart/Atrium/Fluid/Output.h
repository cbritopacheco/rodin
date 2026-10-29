/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ATRIUM_FLUID_OUTPUT_H
#define RODIN_ATRIUM_FLUID_OUTPUT_H

#include <Rodin/IO/XDMF.h>
#include "Fluid.h"

namespace Rodin::Examples::Heart::Atrium
{
  /**
   * @brief Output for the atrium fluid and transport state.
   *
   * ## Architecture
   *
   * The output observes the owning atrium's mesh and fields, writes its
   * distributed XDMF series, and collects scalar diagnostics for the coordinator.
   * The atrium and its MPI context outlive this object.
   */
  class Fluid::Output
  {
    public:
      /**
       * @brief Constructs output bound to an atrium and its communicator.
       * @param atrium Owning atrium.
       * @param context Communicator of the atrium mesh.
       */
      Output(Fluid& atrium, const Context::MPI& context);

      /// @brief Registers the mesh and observed fields with XDMF.
      void bind();

      /**
       * @brief Collects scalar diagnostics for the coordinator CSV.
       * @param cycle Zero-based cycle index.
       * @note Collective on the atrium communicator.
       */
      void collectDiagnostics(int cycle);

      /**
       * @brief Writes and flushes a transient state.
       * @param time Physical time in seconds.
       * @note Collective on the atrium communicator.
       */
      void write(Real time);

      /// @brief Finalizes output on the atrium communicator.
      void close();

    private:
      Fluid& m_atrium;
      IO::XDMF m_xdmf;
  };
}

#endif
