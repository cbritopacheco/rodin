/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ForwardDecls.h
 * @brief Forward declarations for the Assembly module.
 *
 * This file contains forward declarations of all classes in the Assembly module
 * to enable their use before complete definitions are available.
 */
#ifndef RODIN_ASSEMBLY_FORWARDDECLS_H
#define RODIN_ASSEMBLY_FORWARDDECLS_H

#include "Rodin/Variational/ForwardDecls.h"

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
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref AssemblyBase "AssemblyBase<OperatorType, BilinearForm<...>>" | Base interface for assembling one bilinear form into one operator. |
   * | @ref AssemblyBase "AssemblyBase<OperatorType, Tuple<BilinearForm<...>...>>" | Base interface for assembling block/tuple bilinear forms. |
   * | @ref AssemblyBase "AssemblyBase<VectorType, LinearForm<...>>" | Base interface for assembling one linear form into one vector. |
   * | @ref AssemblyBase "AssemblyBase<VectorType, Tuple<LinearForm<...>...>>" | Base interface for assembling block/tuple linear forms. |
   * | @ref AssemblyBase "AssemblyBase<IndexMap<Scalar>, DirichletBC<TrialFunction, FunctionBase>>" | Base interface for scalar Dirichlet elimination maps. |
   * | @ref AssemblyBase "AssemblyBase<IndexMap<pair<IndexArray, Vector>>, DirichletBC<TrialFunction, ShapeFunctionBase>>" | Base interface for identification-style Dirichlet constraints. |
   * | @ref AssemblyBase "AssemblyBase<LinearSystem, Problem<LinearSystem, TrialFunction, TestFunction>>" | Base interface for single-field problem assembly. |
   * | @ref AssemblyBase "AssemblyBase<LinearSystem, Problem<LinearSystem, U1, U2, U3, Us...>>" | Base interface for mixed multi-field problem assembly. |
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
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref Sequential "Sequential<VectorType, LinearForm<FES, VectorType>>" | Sequential assembly of one linear form. |
   * | @ref Sequential "Sequential<VectorType, Tuple<LinearForm<...>...>>" | Sequential assembly of block/tuple linear forms. |
   * | @ref Sequential "Sequential<OperatorType, BilinearForm<...>>" | Sequential assembly of one bilinear form. |
   * | @ref Sequential "Sequential<OperatorType, Tuple<BilinearForm<...>...>>" | Sequential assembly of block/tuple bilinear forms. |
   * | @ref Sequential "Sequential<Math::Vector<Real>, LinearForm<FES, Math::Vector<Real>>>" | Sequential real-vector linear form assembly helper. |
   * | @ref Sequential "Sequential<LinearSystem, Problem<LinearSystem, TrialFunction, TestFunction>>" | Sequential single-field problem assembly. |
   * | @ref Sequential "Sequential<LinearSystem, Problem<LinearSystem, U1, U2, U3, Us...>>" | Sequential mixed multi-field problem assembly. |
   * | @ref Sequential "Sequential<IndexMap<Scalar>, DirichletBC<TrialFunction, FunctionBase>>" | Sequential scalar Dirichlet elimination-map assembly. |
   * | @ref Sequential "Sequential<IndexMap<pair<IndexArray, Vector>>, DirichletBC<TrialFunction, ShapeFunctionBase>>" | Sequential identification-style Dirichlet constraint assembly. |
   * | @ref Sequential "Sequential<Vec, LinearForm<FES, Vec>>" | PETSc sequential linear-form assembly. |
   * | @ref Sequential "Sequential<Mat, BilinearForm<..., Mat>>" | PETSc sequential bilinear-form assembly. |
   * | @ref Sequential "Sequential<PETSc::Math::LinearSystem, Problem<...>>" | PETSc sequential problem assembly. |
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
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref SequentialIteration "SequentialIteration<Mesh<Context::Local>>" | Sequential iteration over local mesh polytopes. |
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
   * @brief Input data for identification-style Dirichlet BC assembly
   *        (`u = A(v)`).
   *
   * Encapsulates the slave trial function @f$ u @f$, the right-hand-side
   * shape-function expression @f$ A(v) @f$, and the essential boundary
   * attributes.
   */
  template <class Scalar, class Sol1, class FES1, class Derived2, class FES2,
    Variational::ShapeFunctionSpaceType Sp>
  class DirichletBCShapeFunctionAssemblyInput;

  /**
   * @brief Default assembly strategy selector.
   *
   * Selects the default assembly implementation based on execution context.
   * When OpenMP is available and enabled, defaults to parallel assembly;
   * otherwise, defaults to sequential assembly.
   *
   * @tparam Ts Context type parameters (e.g., Context::Local)
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref Default "Default<Context::Local>" | Default local assembly strategy. |
   * | @ref Default "Default<Context::Local, Context::Local>" | Default local/local bilinear assembly strategy. |
   * | @ref Default "Default<Context::MPI>" | Default distributed assembly strategy. |
   * | @ref Default "Default<Context::MPI, Context::MPI>" | Default distributed/distributed bilinear assembly strategy. |
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
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref OpenMP "OpenMP<VectorType, LinearForm<FES, VectorType>>" | OpenMP assembly of one linear form. |
   * | @ref OpenMP "OpenMP<VectorType, Tuple<LinearForm<...>...>>" | OpenMP assembly of block/tuple linear forms. |
   * | @ref OpenMP "OpenMP<OperatorType, BilinearForm<...>>" | OpenMP assembly of one bilinear form. |
   * | @ref OpenMP "OpenMP<OperatorType, Tuple<BilinearForm<...>...>>" | OpenMP assembly of block/tuple bilinear forms. |
   * | @ref OpenMP "OpenMP<LinearSystem, Problem<LinearSystem, TrialFunction, TestFunction>>" | OpenMP single-field problem assembly. |
   * | @ref OpenMP "OpenMP<LinearSystem, Problem<LinearSystem, U1, U2, U3, Us...>>" | OpenMP mixed multi-field problem assembly. |
   * | @ref OpenMP "OpenMP<IndexMap<Scalar>, DirichletBC<TrialFunction, FunctionBase>>" | OpenMP scalar Dirichlet elimination-map assembly. |
   * | @ref OpenMP "OpenMP<IndexMap<pair<IndexArray, Vector>>, DirichletBC<TrialFunction, ShapeFunctionBase>>" | OpenMP identification-style Dirichlet constraint assembly. |
   * | @ref OpenMP "OpenMP<Vec, LinearForm<FES, Vec>>" | PETSc/OpenMP linear-form assembly. |
   * | @ref OpenMP "OpenMP<Mat, BilinearForm<..., Mat>>" | PETSc/OpenMP bilinear-form assembly. |
   * | @ref OpenMP "OpenMP<PETSc::Math::LinearSystem, Problem<...>>" | PETSc/OpenMP problem assembly. |
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
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref OpenMPIteration "OpenMPIteration<Mesh<Context::Local>>" | OpenMP iteration over local mesh polytopes. |
   */
  template <class Mesh>
  class OpenMPIteration;
}
#endif

#endif
