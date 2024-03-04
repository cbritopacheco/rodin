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

namespace Rodin::Assembly
{
  template <class OperatorType, class TrialFES, class TestFES>
  class AssemblyBase<OperatorType, Variational::BilinearForm<TrialFES, TestFES, OperatorType>>
    : public FormLanguage::Base
  {
    public:
      using InputType = BilinearFormAssemblyInput<TrialFES, TestFES>;

      AssemblyBase() = default;

      AssemblyBase(const AssemblyBase&) = default;

      AssemblyBase(AssemblyBase&&) = default;

      virtual ~AssemblyBase() = default;

      virtual OperatorType execute(const InputType& data) const = 0;

      virtual AssemblyBase* copy() const noexcept = 0;
  };

  template <class OperatorType, class ... TrialFES, class ... TestFES>
  class AssemblyBase<OperatorType, Tuple<Variational::BilinearForm<TrialFES, TestFES, OperatorType>...>>
  {
    public:
      static_assert(sizeof...(TrialFES) == sizeof...(TestFES));
      using InputType = BilinearFormTupleAssemblyInput<BilinearFormAssemblyInput<TrialFES, TestFES>...>;

      AssemblyBase() = default;

      AssemblyBase(const AssemblyBase&) = default;

      AssemblyBase(AssemblyBase&&) = default;

      virtual ~AssemblyBase() = default;

      virtual OperatorType execute(const InputType& data) const = 0;

      virtual AssemblyBase* copy() const noexcept = 0;
  };

  template <class VectorType, class FES>
  class AssemblyBase<VectorType, Variational::LinearForm<FES, VectorType>>
    : public FormLanguage::Base
  {
    public:
      using InputType = LinearFormAssemblyInput<FES>;

      AssemblyBase() = default;

      AssemblyBase(const AssemblyBase&) = default;

      AssemblyBase(AssemblyBase&&) = default;

      virtual ~AssemblyBase() = default;

      virtual VectorType execute(const InputType& data) const = 0;

      virtual AssemblyBase* copy() const noexcept = 0;
  };
}

#endif
