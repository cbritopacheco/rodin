/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ASSEMBLY_ASSEMBLYBASE_H
#define RODIN_ASSEMBLY_ASSEMBLYBASE_H

#include "Rodin/Math.h"
#include "Rodin/Tuple.h"

#include "ForwardDecls.h"
#include "Input.h"
#include "Rodin/Variational/ProblemBody.h"

namespace Rodin::Assembly
{
  /**
   * @brief Base class for bilinear form assembly operations.
   *
   * This template specialization handles the assembly of bilinear forms into
   * matrix operators. It provides the interface for converting variational
   * forms @f$ a(u,v) @f$ into discrete matrix representations @f$ A @f$ where
   * @f$ A_{ij} = a(\phi_j, \psi_i) @f$ for trial functions @f$ \phi_j @f$ and
   * test functions @f$ \psi_i @f$.
   *
   * @tparam OperatorType Matrix type for the assembled operator (e.g., sparse matrix)
   * @tparam Solution Solution variable type 
   * @tparam TrialFES Trial finite element space type
   * @tparam TestFES Test finite element space type
   */
  template <class OperatorType, class Solution, class TrialFES, class TestFES>
  class AssemblyBase<OperatorType, Variational::BilinearForm<Solution, TrialFES, TestFES, OperatorType>>
    : public FormLanguage::Base
  {
    public:
      /// @brief Input data type for bilinear form assembly
      using InputType = BilinearFormAssemblyInput<TrialFES, TestFES>;

      /// @brief Default constructor
      AssemblyBase() = default;

      /// @brief Copy constructor
      AssemblyBase(const AssemblyBase&) = default;

      /// @brief Move constructor  
      AssemblyBase(AssemblyBase&&) = default;

      /// @brief Virtual destructor
      virtual ~AssemblyBase() = default;

      /**
       * @brief Executes the assembly operation.
       *
       * @param out Output matrix to store the assembled bilinear form
       * @param data Assembly input data containing finite element spaces and integration data
       *
       * Performs the core assembly operation, computing the matrix entries
       * @f$ A_{ij} = a(\phi_j, \psi_i) @f$ where @f$ a(\cdot,\cdot) @f$ is
       * the bilinear form, @f$ \phi_j @f$ are trial basis functions, and
       * @f$ \psi_i @f$ are test basis functions.
       */
      virtual void execute(OperatorType& out, const InputType& data) const = 0;

      /**
       * @brief Creates a polymorphic copy of this assembly object.
       * 
       * @return AssemblyBase* Pointer to a new copy of this object
       */
      virtual AssemblyBase* copy() const noexcept = 0;
  };

  template <class OperatorType, class ... Solution, class ... TrialFES, class ... TestFES, class ... BlockType>
  class AssemblyBase<OperatorType, Tuple<Variational::BilinearForm<Solution, TrialFES, TestFES, BlockType>...>>
  {
    public:
      static_assert(sizeof...(TrialFES) == sizeof...(TestFES));
      using InputType =
        BilinearFormTupleAssemblyInput<BilinearFormAssemblyInput<TrialFES, TestFES>...>;

      AssemblyBase() = default;

      AssemblyBase(const AssemblyBase&) = default;

      AssemblyBase(AssemblyBase&&) = default;

      virtual ~AssemblyBase() = default;

      virtual void execute(OperatorType& out, const InputType& data) const = 0;

      virtual AssemblyBase* copy() const noexcept = 0;
  };

  /**
   * @brief Base class for linear form assembly operations.
   *
   * This template specialization handles the assembly of linear forms into
   * vector objects. It provides the interface for converting variational
   * forms @f$ l(v) @f$ into discrete vector representations @f$ b @f$ where
   * @f$ b_i = l(\psi_i) @f$ for test functions @f$ \psi_i @f$.
   *
   * @tparam VectorType Vector type for the assembled right-hand side
   * @tparam FES Finite element space type
   */
  template <class VectorType, class FES>
  class AssemblyBase<VectorType, Variational::LinearForm<FES, VectorType>>
    : public FormLanguage::Base
  {
    public:
      /// @brief Input data type for linear form assembly
      using InputType = LinearFormAssemblyInput<FES>;

      /// @brief Default constructor
      AssemblyBase() = default;

      /// @brief Copy constructor
      AssemblyBase(const AssemblyBase&) = default;

      /// @brief Move constructor
      AssemblyBase(AssemblyBase&&) = default;

      /// @brief Virtual destructor
      virtual ~AssemblyBase() = default;

      /**
       * @brief Executes the linear form assembly operation.
       *
       * @param out Output vector to store the assembled linear form
       * @param data Assembly input data containing finite element space and integration data
       *
       * Performs the core assembly operation, computing the vector entries
       * @f$ b_i = l(\psi_i) @f$ where @f$ l(\cdot) @f$ is the linear form
       * and @f$ \psi_i @f$ are test basis functions.
       */
      virtual void execute(VectorType& out, const InputType& data) const = 0;

      /**
       * @brief Creates a polymorphic copy of this assembly object.
       * 
       * @return AssemblyBase* Pointer to a new copy of this object
       */
      virtual AssemblyBase* copy() const noexcept = 0;
  };

  template <class VectorType, class ... FES, class ... BlockType>
  class AssemblyBase<VectorType, Tuple<Variational::LinearForm<FES, BlockType>...>>
  {
    public:
      static_assert(sizeof...(FES) == sizeof...(FES));
      using InputType =
        LinearFormTupleAssemblyInput<LinearFormAssemblyInput<FES>...>;

      AssemblyBase() = default;

      AssemblyBase(const AssemblyBase&) = default;

      AssemblyBase(AssemblyBase&&) = default;

      virtual ~AssemblyBase() = default;

      virtual void execute(VectorType& out, const InputType& data) const = 0;

      virtual AssemblyBase* copy() const noexcept = 0;
  };

  template <class Scalar, class Solution, class FES, class ValueDerived>
  class AssemblyBase<
    IndexMap<Scalar>,
    Variational::DirichletBC<
      Variational::TrialFunction<Solution, FES>, Variational::FunctionBase<ValueDerived>>>
    : public FormLanguage::Base
  {
    public:
      using ScalarType = Scalar;

      using ValueType = Variational::FunctionBase<ValueDerived>;

      using InputType = DirichletBCAssemblyInput<Scalar, Solution, FES, ValueType>;

      AssemblyBase() = default;

      AssemblyBase(const AssemblyBase&) = default;

      AssemblyBase(AssemblyBase&&) = default;

      virtual ~AssemblyBase() = default;

      virtual void execute(IndexMap<ScalarType>& out, const InputType& data) const = 0;

      virtual AssemblyBase* copy() const noexcept = 0;
  };

  template <class LinearSystem, class TrialFunction, class TestFunction>
  class AssemblyBase<LinearSystem, Variational::Problem<LinearSystem, TrialFunction, TestFunction>>
    : public FormLanguage::Base
  {
    public:
      using OperatorType =
        typename FormLanguage::Traits<LinearSystem>::OperatorType;

      using VectorType =
        typename FormLanguage::Traits<LinearSystem>::VectorType;

      using ScalarType =
        typename FormLanguage::Traits<LinearSystem>::ScalarType;

      using ProblemBodyType =
        Variational::ProblemBody<OperatorType, VectorType, ScalarType>;

      using InputType =
        ProblemAssemblyInput<ProblemBodyType, TrialFunction, TestFunction>;

      using Parent =
        FormLanguage::Base;

      AssemblyBase() = default;

      AssemblyBase(const AssemblyBase& other)
        : Parent(other)
      {}

      AssemblyBase(AssemblyBase&& other)
        : Parent(std::move(other))
      {}

      virtual ~AssemblyBase() = default;

      virtual void execute(LinearSystem& out, const InputType& input) const = 0;

      virtual AssemblyBase* copy() const noexcept = 0;
  };
}

#endif
