/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Input.h
 * @brief Input data structures for assembly operations.
 *
 * This file defines the input data structures that encapsulate all information
 * required for assembling variational forms into discrete linear algebra objects.
 * These structures provide a consistent interface for passing finite element spaces,
 * integrators, and boundary conditions to assembly implementations.
 */
#ifndef RODIN_ASSEMBLY_INPUT_H
#define RODIN_ASSEMBLY_INPUT_H

#include <functional>

#include "Rodin/Math.h"
#include "Rodin/Tuple.h"

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/FormLanguage/List.h"
#include "Rodin/FormLanguage/Base.h"
#include "Rodin/Variational/ForwardDecls.h"

#include "ForwardDecls.h"
#include "Rodin/Variational/LinearFormIntegrator.h"

namespace Rodin::Assembly
{
  /**
   * @ingroup RodinAssembly
   * @brief Input data for bilinear form assembly.
   *
   * This class encapsulates all data required to assemble a bilinear form
   * @f$ a(u,v) @f$ into a matrix operator. It stores references to the trial
   * and test finite element spaces, as well as lists of local and global
   * bilinear form integrators.
   *
   * ## Local vs Global Integrators
   * - **Local integrators**: Operate on individual elements @f$ K @f$, computing
   *   contributions @f$ a_K(\phi_j, \psi_i) @f$ for basis functions with support
   *   on @f$ K @f$.
   * - **Global integrators**: Couple multiple elements, useful for non-local
   *   operators or interface conditions.
   *
   * @tparam TrialFES Trial finite element space type
   * @tparam TestFES Test finite element space type
   *
   * ## Usage Example
   * ```cpp
   * P1 Uh(mesh), Vh(mesh);
   * auto lbfis = // list of local integrators
   * auto gbfis = // list of global integrators  
   * BilinearFormAssemblyInput input(Uh, Vh, lbfis, gbfis);
   * ```
   */
  template <class TrialFES, class TestFES>
  class BilinearFormAssemblyInput
  {
    public:
      /// @brief Trial finite element space type
      using TrialFESType = TrialFES;

      /// @brief Test finite element space type
      using TestFESType = TestFES;

      /// @brief Scalar type for trial space
      using TrialFESScalarType  = typename FormLanguage::Traits<TrialFESType>::ScalarType;

      /// @brief Scalar type for test space
      using TestFESScalarType   = typename FormLanguage::Traits<TestFESType>::ScalarType;

      /// @brief Resulting scalar type for assembled operator
      using ScalarType = decltype(
          std::declval<TrialFESScalarType>() * std::declval<TestFESScalarType>());

      /// @brief Type for local bilinear form integrators
      using LocalBilinearFormIntegratorBaseType =
        Variational::LocalBilinearFormIntegratorBase<ScalarType>;

      /// @brief List type for local integrators
      using LocalBilinearFormIntegratorBaseListType =
        FormLanguage::List<LocalBilinearFormIntegratorBaseType>;

      /// @brief Type for global bilinear form integrators
      using GlobalBilinearFormIntegratorBaseType =
        Variational::GlobalBilinearFormIntegratorBase<ScalarType>;

      /// @brief List type for global integrators
      using GlobalBilinearFormIntegratorBaseListType =
        FormLanguage::List<GlobalBilinearFormIntegratorBaseType>;

      /**
       * @brief Constructs bilinear form assembly input.
       *
       * @param trialFES Trial finite element space
       * @param testFES Test finite element space
       * @param lbfis List of local bilinear form integrators
       * @param gbfis List of global bilinear form integrators
       */
      BilinearFormAssemblyInput(
          const TrialFES& trialFES, const TestFES& testFES,
          LocalBilinearFormIntegratorBaseListType& lbfis,
          GlobalBilinearFormIntegratorBaseListType& gbfis)
        : m_trialFES(trialFES), m_testFES(testFES), m_lbfis(lbfis), m_gbfis(gbfis)
      {}

      /**
       * @brief Gets the trial finite element space.
       * @return Reference to trial space
       */
      const TrialFES& getTrialFES() const
      {
        return m_trialFES.get();
      }

      /**
       * @brief Gets the test finite element space.
       * @return Reference to test space
       */
      const TestFES& getTestFES() const
      {
        return m_testFES.get();
      }

      /**
       * @brief Gets the list of local bilinear form integrators.
       * @return Reference to local integrator list
       */
      LocalBilinearFormIntegratorBaseListType& getLocalBFIs() const
      {
        return m_lbfis.get();
      }

      /**
       * @brief Gets the list of global bilinear form integrators.
       * @return Reference to global integrator list
       */
      GlobalBilinearFormIntegratorBaseListType& getGlobalBFIs() const
      {
        return m_gbfis.get();
      }

    private:
      std::reference_wrapper<const TrialFES>  m_trialFES;   ///< Trial space reference
      std::reference_wrapper<const TestFES>   m_testFES;    ///< Test space reference
      std::reference_wrapper<LocalBilinearFormIntegratorBaseListType>   m_lbfis;   ///< Local integrators
      std::reference_wrapper<GlobalBilinearFormIntegratorBaseListType>  m_gbfis;   ///< Global integrators
  };

  /// @brief Template argument deduction guide for BilinearFormAssemblyInput
  template <class TrialFES, class TestFES>
  BilinearFormAssemblyInput(
      const TrialFES&, const TestFES&,
      FormLanguage::List<
        Variational::LocalBilinearFormIntegratorBase<
          decltype(
            std::declval<typename FormLanguage::Traits<TrialFES>::ScalarType>() *
            std::declval<typename FormLanguage::Traits<TestFES>::ScalarType>())>>&,
      FormLanguage::List<
        Variational::GlobalBilinearFormIntegratorBase<
          decltype(
            std::declval<typename FormLanguage::Traits<TrialFES>::ScalarType>() *
            std::declval<typename FormLanguage::Traits<TestFES>::ScalarType>())>>&)
    -> BilinearFormAssemblyInput<TrialFES, TestFES>;

  /**
   * @ingroup RodinAssembly
   * @brief Input data for linear form assembly.
   *
   * This class encapsulates all data required to assemble a linear form
   * @f$ l(v) @f$ into a vector. It stores a reference to the finite element
   * space and a list of linear form integrators.
   *
   * Linear form integrators compute element contributions @f$ l_K(\psi_i) @f$
   * for test basis functions @f$ \psi_i @f$ with support on element @f$ K @f$.
   *
   * @tparam FES Finite element space type
   *
   * ## Usage Example
   * ```cpp
   * P1 Vh(mesh);
   * auto lfis = // list of linear form integrators
   * LinearFormAssemblyInput input(Vh, lfis);
   * ```
   */
  template <class FES>
  class LinearFormAssemblyInput
  {
    public:
      /// @brief Finite element space type
      using FESType = FES;

      /// @brief Scalar type for the space
      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      /// @brief List type for linear form integrators
      using LinearFormIntegratorBaseList =
        FormLanguage::List<Variational::LinearFormIntegratorBase<ScalarType>>;

      /**
       * @brief Constructs linear form assembly input.
       *
       * @param fes Finite element space
       * @param lfis List of linear form integrators
       */
      LinearFormAssemblyInput(const FES& fes, LinearFormIntegratorBaseList& lfis)
        : m_fes(fes), m_lfis(lfis)
      {}

      /**
       * @brief Gets the finite element space.
       * @return Reference to finite element space
       */
      const FES& getFES() const
      {
        return m_fes.get();
      }

      /**
       * @brief Gets the list of linear form integrators.
       * @return Reference to integrator list
       */
      LinearFormIntegratorBaseList& getLFIs() const
      {
        return m_lfis.get();
      }

    private:
      std::reference_wrapper<const FES> m_fes;                                ///< Finite element space reference
      std::reference_wrapper<LinearFormIntegratorBaseList> m_lfis;            ///< Linear form integrators
  };

  /**
   * @ingroup RodinAssembly
   * @brief Input data for tuple of bilinear forms assembly.
   *
   * This class handles assembly of multiple bilinear forms organized as a tuple,
   * typically used for mixed finite element formulations or block systems.
   * It manages the layout and offsets for assembling block matrices.
   *
   * For a mixed problem with spaces @f$ V_h = V_{h,1} \times V_{h,2} @f$ and
   * @f$ W_h = W_{h,1} \times W_{h,2} @f$, this assembles the block matrix:
   * @f[
   *   A = \begin{bmatrix}
   *     A_{11} & A_{12} \\
   *     A_{21} & A_{22}
   *   \end{bmatrix}
   * @f]
   *
   * @tparam Ts Input types for individual bilinear form blocks
   */
  template <class ... Ts>
  class BilinearFormTupleAssemblyInput
  {
    public:
      /// @brief Number of blocks in the tuple
      static constexpr size_t Size = sizeof...(Ts);

      /// @brief Array type for row and column offsets
      using Offsets = std::array<Pair<size_t, size_t>, sizeof...(Ts)>;

      /**
       * @brief Constructs tuple assembly input.
       *
       * @param rows Total number of rows in block matrix
       * @param cols Total number of columns in block matrix
       * @param offsets Row and column offsets for each block
       * @param ins Tuple of individual bilinear form inputs
       */
      BilinearFormTupleAssemblyInput(
          size_t rows, size_t cols,
          const Offsets& offsets,
          const Tuple<Ts...>& ins)
        : m_rows(rows), m_cols(cols), m_offsets(offsets), m_ins(ins)
      {}

      /**
       * @brief Gets the total number of rows.
       * @return Total row count
       */
      size_t getRows() const
      {
        return m_rows;
      }

      /**
       * @brief Gets the total number of columns.
       * @return Total column count
       */
      size_t getColumns() const
      {
        return m_cols;
      }

      /**
       * @brief Gets the block offsets.
       * @return Array of (column_offset, row_offset) pairs for each block
       */
      const Offsets& getOffsets() const
      {
        return m_offsets;
      }

      /**
       * @brief Gets the tuple of input data.
       * @return Tuple of bilinear form assembly inputs
       */
      const Tuple<Ts...>& getTuple() const
      {
        return m_ins;
      }

    private:
      size_t m_rows, m_cols;     ///< Total matrix dimensions
      Offsets m_offsets;         ///< Block offsets
      const Tuple<Ts...> m_ins;  ///< Individual block inputs
  };

  /**
   * @ingroup RodinAssembly
   * @brief Input data for tuple of linear forms assembly.
   *
   * This class handles assembly of multiple linear forms organized as a tuple,
   * typically used for mixed finite element formulations. It manages the layout
   * and offsets for assembling block vectors.
   *
   * For a mixed problem with test space @f$ W_h = W_{h,1} \times W_{h,2} @f$,
   * this assembles the block vector:
   * @f[
   *   b = \begin{bmatrix} b_1 \\ b_2 \end{bmatrix}
   * @f]
   *
   * @tparam Ts Input types for individual linear form blocks
   */
  template <class ... Ts>
  class LinearFormTupleAssemblyInput
  {
    public:
      /// @brief Number of blocks in the tuple
      static constexpr size_t Size = sizeof...(Ts);

      /// @brief Array type for block offsets
      using Offsets = std::array<size_t, sizeof...(Ts)>;

      /**
       * @brief Constructs tuple assembly input.
       *
       * @param size Total size of block vector
       * @param offsets Starting indices for each block
       * @param ins Tuple of individual linear form inputs
       */
      LinearFormTupleAssemblyInput(
          size_t size, const Offsets& offsets, const Tuple<Ts...>& ins)
        : m_size(size), m_offsets(offsets), m_ins(ins)
      {}

      /**
       * @brief Gets the total vector size.
       * @return Total size
       */
      size_t getSize() const
      {
        return m_size;
      }

      /**
       * @brief Gets the block offsets.
       * @return Array of starting indices for each block
       */
      const Offsets& getOffsets() const
      {
        return m_offsets;
      }

      /**
       * @brief Gets the tuple of input data.
       * @return Tuple of linear form assembly inputs
       */
      const Tuple<Ts...>& getTuple() const
      {
        return m_ins;
      }

    private:
      size_t m_size;             ///< Total vector size
      Offsets m_offsets;         ///< Block offsets
      const Tuple<Ts...> m_ins;  ///< Individual block inputs
  };

  /**
   * @ingroup RodinAssembly
   * @brief Input data for Dirichlet boundary condition assembly.
   *
   * This class encapsulates all data required to assemble essential (Dirichlet)
   * boundary conditions. It stores the trial function, the prescribed boundary
   * value function, and the set of boundary attributes where conditions apply.
   *
   * For a Dirichlet boundary condition @f$ u = g @f$ on @f$ \Gamma_D @f$,
   * this assembles the map of constrained degrees of freedom to their values.
   *
   * @tparam Scalar Scalar type for assembled boundary values
   * @tparam Solution Solution variable type
   * @tparam FES Finite element space type
   * @tparam Value Type of boundary value function
   */
  template <class Scalar, class Solution, class FES, class Value>
  class DirichletBCAssemblyInput
  {
    public:
      /// @brief Boundary value function type
      using ValueType = Value;

      /// @brief Trial function type
      using OperandType = Variational::TrialFunction<Solution, FES>;

      /// @brief Dirichlet boundary condition type
      using DirichletBCType = Variational::DirichletBC<OperandType, ValueType>;

      /**
       * @brief Constructs Dirichlet BC assembly input.
       *
       * @param u Trial function to constrain
       * @param value Boundary value function @f$ g @f$
       * @param essBdr Set of boundary attributes for essential boundary @f$ \Gamma_D @f$
       */
      DirichletBCAssemblyInput(
          const OperandType& u, const ValueType& value, const FlatSet<Geometry::Attribute>& essBdr)
        : m_u(u), m_value(value), m_essBdr(essBdr)
      {}

      /**
       * @brief Gets the trial function operand.
       * @return Reference to trial function
       */
      const OperandType& getOperand() const
      {
        return m_u.get();
      }

      /**
       * @brief Gets the boundary value function.
       * @return Reference to value function
       */
      const ValueType& getValue() const
      {
        return m_value.get();
      }

      /**
       * @brief Gets the essential boundary attributes.
       * @return Set of boundary attributes where condition applies
       */
      const FlatSet<Geometry::Attribute>& getEssentialBoundary() const
      {
        return m_essBdr.get();
      }

    private:
      std::reference_wrapper<const OperandType> m_u;          ///< Trial function
      std::reference_wrapper<const ValueType> m_value;        ///< Boundary value
      std::reference_wrapper<const FlatSet<Geometry::Attribute>> m_essBdr;  ///< Boundary attributes
  };

  /**
   * @ingroup RodinAssembly
   * @brief Input data for complete problem assembly.
   *
   * This class encapsulates all data required to assemble a complete variational
   * problem into a linear system. It stores references to the problem body
   * (containing all forms and boundary conditions), trial function, and test
   * function.
   *
   * The problem body contains:
   * - Local and global bilinear form integrators
   * - Linear form integrators
   * - Dirichlet boundary conditions
   * - Periodic boundary conditions
   * - Pre-assembled bilinear and linear forms
   *
   * @tparam ProblemBody Type of problem body containing forms and BCs
   * @tparam TrialFunction Trial function type
   * @tparam TestFunction Test function type
   */
  template <class ProblemBody, class TrialFunction, class TestFunction>
  class ProblemAssemblyInput
  {
    public:
      /// @brief Problem body type
      using ProblemBodyType = ProblemBody;

      /**
       * @brief Constructs problem assembly input.
       *
       * @param body Problem body containing all variational forms and boundary conditions
       * @param trialFunction Trial function for the problem
       * @param testFunction Test function for the problem
       */
      ProblemAssemblyInput(
          ProblemBody& body, const TrialFunction& trialFunction, const TestFunction& testFunction)
        : m_body(body),
          m_trialFunction(trialFunction),
          m_testFunction(testFunction)
      {}

      /**
       * @brief Gets the problem body.
       * @return Reference to problem body
       */
      ProblemBody& getProblemBody() const
      {
        return m_body.get();
      }

      /**
       * @brief Gets the trial function.
       * @return Reference to trial function
       */
      const TrialFunction& getTrialFunction() const
      {
        return m_trialFunction.get();
      }

      /**
       * @brief Gets the test function.
       * @return Reference to test function
       */
      const TestFunction& getTestFunction() const
      {
        return m_testFunction.get();
      }

    private:
      std::reference_wrapper<ProblemBody> m_body;                           ///< Problem body reference
      std::reference_wrapper<const TrialFunction> m_trialFunction;          ///< Trial function reference
      std::reference_wrapper<const TestFunction> m_testFunction;            ///< Test function reference
  };
}

#endif
