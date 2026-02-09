/**
 * @file ForwardDecls.h
 * @brief Forward declarations for the Assembly module.
 *
 * This file contains forward declarations of all classes in the Assembly module
 * to enable their use before complete definitions are available.
 */
#ifndef RODIN_ASSEMBLY_FORWARDDECLS_H
#define RODIN_ASSEMBLY_FORWARDDECLS_H

namespace Rodin::Assembly
{
  /**
   * @brief Base class template for assembly operations.
   *
   * Provides the interface for assembling variational forms into discrete
   * linear algebra objects (matrices and vectors). Different specializations
   * handle bilinear forms, linear forms, boundary conditions, and problems.
   *
   * @tparam LinearAlgebraType Target linear algebra type (e.g., sparse matrix, vector)
   * @tparam Operand Variational form type (e.g., BilinearForm, LinearForm)
   */
  template <class LinearAlgebraType, class Operand>
  class AssemblyBase;

  /**
   * @brief Sequential (single-threaded) assembly implementation.
   *
   * Implements assembly operations using sequential iteration over mesh elements.
   * Provides deterministic, reproducible results suitable for debugging and
   * verification purposes.
   *
   * @tparam LinearAlgebraType Target linear algebra type
   * @tparam Operand Variational form type to assemble
   */
  template <class LinearAlgebraType, class Operand>
  class Sequential;

  /**
   * @brief Sequential mesh iteration strategy.
   *
   * Provides sequential iteration over mesh polytopes for assembly operations.
   * Used internally by Sequential assembly implementations.
   *
   * @tparam Mesh Mesh type to iterate over
   */
  template <class Mesh>
  class SequentialIteration;

  /**
   * @brief Input data for bilinear form assembly.
   *
   * Encapsulates trial and test finite element spaces along with local and
   * global bilinear form integrators required for assembly.
   *
   * @tparam TrialFES Trial finite element space type
   * @tparam TestFES Test finite element space type
   */
  template <class TrialFES, class TestFES>
  class BilinearFormAssemblyInput;

  /**
   * @brief Input data for linear form assembly.
   *
   * Encapsulates the finite element space and linear form integrators
   * required for assembly.
   *
   * @tparam FES Finite element space type
   */
  template <class FES>
  class LinearFormAssemblyInput;

  /**
   * @brief Input data for tuple of bilinear forms assembly.
   *
   * Handles assembly of multiple bilinear forms organized as a tuple,
   * typically used for mixed or block formulations.
   *
   * @tparam Ts Input types for individual bilinear forms
   */
  template <class ... Ts>
  class BilinearFormTupleAssemblyInput;

  /**
   * @brief Input data for Dirichlet boundary condition assembly.
   *
   * Encapsulates the trial function, boundary value function, and essential
   * boundary attributes for Dirichlet boundary condition assembly.
   *
   * @tparam Scalar Scalar type for boundary values
   * @tparam Solution Solution variable type
   * @tparam FES Finite element space type
   * @tparam ValueDerived Type of boundary value function
   */
  template <class Scalar, class Solution, class FES, class ValueDerived>
  class DirichletBCAssemblyInput;

  /**
   * @brief Default assembly strategy selector.
   *
   * Selects the default assembly implementation based on execution context.
   * When OpenMP is available and enabled, defaults to parallel assembly;
   * otherwise, defaults to sequential assembly.
   *
   * @tparam Ts Context type parameters (e.g., Context::Local)
   */
  template <class ... Ts>
  class Default;
}

#ifdef RODIN_USE_OPENMP
namespace Rodin::Assembly
{
  /**
   * @brief OpenMP-based parallel assembly implementation.
   *
   * Implements multi-threaded assembly operations using OpenMP parallelization.
   * Distributes element computations across threads for improved performance
   * on shared-memory systems.
   *
   * @tparam LinearAlgebraType Target linear algebra type
   * @tparam Operand Variational form type to assemble
   */
  template <class LinearAlgebraType, class Operand>
  class OpenMP;

  /**
   * @brief OpenMP-based parallel mesh iteration strategy.
   *
   * Provides parallel iteration over mesh polytopes using OpenMP threading.
   * Used internally by OpenMP assembly implementations.
   *
   * @tparam Mesh Mesh type to iterate over
   */
  template <class Mesh>
  class OpenMPIteration;
}
#endif

#endif
