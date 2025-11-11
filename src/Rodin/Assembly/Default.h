/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Default.h
 * @brief Default assembly strategy selection based on compilation configuration.
 *
 * This file provides compile-time selection of assembly strategies based on
 * the execution context and available parallelization support. When OpenMP
 * is available and enabled (RODIN_USE_OPENMP), parallel assembly is used by
 * default; otherwise, sequential assembly is selected.
 *
 * The Default class acts as a policy selector that maps execution contexts
 * to appropriate assembly implementations without requiring explicit user
 * specification.
 */
#ifndef RODIN_ASSEMBLY_DEFAULT_H
#define RODIN_ASSEMBLY_DEFAULT_H

#include "Rodin/Configure.h"

#include "Rodin/Context/Local.h"

#include "ForwardDecls.h"

#ifdef RODIN_USE_OPENMP

#include "OpenMP.h"

namespace Rodin::Assembly
{
  /**
   * @ingroup RodinAssembly
   * @brief Default assembly strategy for local context with OpenMP enabled.
   *
   * When OpenMP support is compiled in (RODIN_USE_OPENMP is defined), this
   * specialization selects the OpenMP parallel assembly implementation as
   * the default for local execution contexts.
   *
   * This provides automatic parallelization of assembly operations without
   * requiring users to explicitly choose the parallel implementation.
   */
  template <>
  class Default<Context::Local>
  {
    public:
      /**
       * @brief Default assembly type for local context (OpenMP-enabled).
       *
       * Type alias that selects OpenMP parallel assembly as the default
       * implementation for assembling variational forms.
       *
       * @tparam LinearAlgebraType Target linear algebra type (matrix, vector, etc.)
       * @tparam Object Variational form type to assemble
       */
      template <class LinearAlgebraType, class Object>
      using Type = OpenMP<LinearAlgebraType, Object>;
  };

  /**
   * @ingroup RodinAssembly
   * @brief Default assembly strategy for mixed local contexts with OpenMP enabled.
   *
   * Specialization for problems with separate trial and test spaces in local
   * execution contexts. Selects OpenMP parallel assembly by default.
   */
  template <>
  class Default<Context::Local, Context::Local>
  {
    public:
      /**
       * @brief Default assembly type for mixed local contexts (OpenMP-enabled).
       *
       * Type alias for assembling bilinear forms with potentially different
       * trial and test spaces, both in local execution contexts.
       *
       * @tparam LinearAlgebraType Target linear algebra type
       * @tparam Object Variational form type to assemble
       */
      template <class LinearAlgebraType, class Object>
      using Type = OpenMP<LinearAlgebraType, Object>;
 };
}

#else

#include "Sequential.h"

namespace Rodin::Assembly
{
  /**
   * @ingroup RodinAssembly
   * @brief Default assembly strategy for local context without OpenMP.
   *
   * When OpenMP support is not available (RODIN_USE_OPENMP is not defined),
   * this specialization selects the Sequential single-threaded assembly
   * implementation as the default for local execution contexts.
   *
   * This ensures assembly operations work correctly even without parallel
   * support, providing deterministic and reproducible results.
   */
  template <>
  class Default<Context::Local>
  {
    public:
      /**
       * @brief Default assembly type for local context (sequential).
       *
       * Type alias that selects Sequential single-threaded assembly when
       * OpenMP is not available.
       *
       * @tparam LinearAlgebraType Target linear algebra type (matrix, vector, etc.)
       * @tparam Object Variational form type to assemble
       */
      template <class LinearAlgebraType, class Object>
      using Type = Sequential<LinearAlgebraType, Object>;
  };

  /**
   * @ingroup RodinAssembly
   * @brief Default assembly strategy for mixed local contexts without OpenMP.
   *
   * Specialization for problems with separate trial and test spaces in local
   * execution contexts. Selects Sequential assembly by default when OpenMP
   * is not available.
   */
  template <>
  class Default<Context::Local, Context::Local>
  {
    public:
      /**
       * @brief Default assembly type for mixed local contexts (sequential).
       *
       * Type alias for assembling bilinear forms with potentially different
       * trial and test spaces, using single-threaded Sequential assembly.
       *
       * @tparam LinearAlgebraType Target linear algebra type
       * @tparam Object Variational form type to assemble
       */
      template <class LinearAlgebraType, class Object>
      using Type = Sequential<LinearAlgebraType, Object>;
 };
}

#endif

#endif
