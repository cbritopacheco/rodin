#ifndef RODIN_MPI_ASSEMBLY_DEFAULT_H
#define RODIN_MPI_ASSEMBLY_DEFAULT_H

/**
 * @file
 * @brief Default assembly policy specializations for MPI contexts.
 *
 * This file wires @ref Rodin::Context::MPI to the distributed assembly
 * backend @ref Rodin::Assembly::MPI through @ref Rodin::Assembly::Default.
 */

#include "Rodin/MPI/Context.h"
#include "Rodin/Assembly/Default.h"

#include "MPI.h"

namespace Rodin::Assembly
{
  /**
   * @brief Selects MPI assembly when only the trial/test context is MPI.
   */
  template <>
  class Default<Context::MPI>
  {
    public:
      /**
       * @brief Assembly-backend selector for MPI trial/test contexts.
       *
       * @tparam LinearAlgebraType Backend linear-algebra container type.
       * @tparam Object Assembly operand type.
       */
      template <class LinearAlgebraType, class Object>
      using Type = MPI<LinearAlgebraType, Object>;
  };

  /**
   * @brief Selects MPI assembly when both operand and assembly contexts are MPI.
   */
  template <>
  class Default<Context::MPI, Context::MPI>
  {
    public:
      /**
       * @brief Assembly-backend selector for fully MPI-distributed assembly.
       *
       * @tparam LinearAlgebraType Backend linear-algebra container type.
       * @tparam Object Assembly operand type.
       */
      template <class LinearAlgebraType, class Object>
      using Type = MPI<LinearAlgebraType, Object>;
  };
}

#endif
