/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ASSEMBLY_INPUT_H
#define RODIN_ASSEMBLY_INPUT_H

#include "Rodin/Math.h"
#include "Rodin/Tuple.h"

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/FormLanguage/List.h"
#include "Rodin/Variational/ForwardDecls.h"

#include "ForwardDecls.h"

namespace Rodin::Assembly
{
  template <class TrialFES, class TestFES>
  class BilinearFormAssemblyInput
  {
    public:
      BilinearFormAssemblyInput(
          const TrialFES& trialFES, const TestFES& testFES,
          FormLanguage::List<Variational::LocalBilinearFormIntegratorBase>& lbfis,
          FormLanguage::List<Variational::GlobalBilinearFormIntegratorBase>& gbfis)
        : m_trialFES(trialFES), m_testFES(testFES), m_lbfis(lbfis), m_gbfis(gbfis)
      {}

      const TrialFES& getTrialFES() const
      {
        return m_trialFES.get();
      }

      const TestFES& getTestFES() const
      {
        return m_testFES.get();
      }

      FormLanguage::List<Variational::LocalBilinearFormIntegratorBase>& getLocalBFIs() const
      {
        return m_lbfis.get();
      }

      FormLanguage::List<Variational::GlobalBilinearFormIntegratorBase>& getGlobalBFIs() const
      {
        return m_gbfis.get();
      }

    private:
      std::reference_wrapper<const TrialFES> m_trialFES;
      std::reference_wrapper<const TestFES> m_testFES;
      std::reference_wrapper<FormLanguage::List<Variational::LocalBilinearFormIntegratorBase>> m_lbfis;
      std::reference_wrapper<FormLanguage::List<Variational::GlobalBilinearFormIntegratorBase>> m_gbfis;
  };

  template <class TrialFES, class TestFES>
  BilinearFormAssemblyInput(
      const TrialFES&, const TestFES&,
      FormLanguage::List<Variational::LocalBilinearFormIntegratorBase>&,
      FormLanguage::List<Variational::GlobalBilinearFormIntegratorBase>&)
    -> BilinearFormAssemblyInput<TrialFES, TestFES>;

  template <class FES>
  class LinearFormAssemblyInput
  {
    public:
      LinearFormAssemblyInput(
          const FES& fes, FormLanguage::List<Variational::LinearFormIntegratorBase>& lfis)
        : m_fes(fes), m_lfis(lfis)
      {}

      const FES& getFES() const
      {
        return m_fes.get();
      }

      FormLanguage::List<Variational::LinearFormIntegratorBase>& getLFIs() const
      {
        return m_lfis.get();
      }

    private:
      std::reference_wrapper<const FES> m_fes;
      std::reference_wrapper<FormLanguage::List<Variational::LinearFormIntegratorBase>> m_lfis;
  };

  template <class ... Ts>
  class BilinearFormTupleAssemblyInput
  {
    public:
      static constexpr size_t Size = sizeof...(Ts);

      using Sizes = std::array<Pair<size_t, size_t>, Size>;

      using Offsets = Sizes;

      BilinearFormTupleAssemblyInput(
          const Sizes& sizes, const Offsets& offsets, const Tuple<Ts...>& ins)
        : m_sizes(sizes), m_offsets(offsets), m_ins(ins)
      {}

      const Sizes& getSizes() const
      {
        return m_sizes;
      }

      const Offsets& getOffsets() const
      {
        return m_offsets;
      }

      const Tuple<Ts...>& getTuple() const
      {
        return m_ins;
      }

    private:
      Sizes m_sizes;
      Offsets m_offsets;
      const Tuple<Ts...> m_ins;
  };

  template <class ... Ts>
  class LinearFormTupleAssemblyInput
  {
    public:
      static constexpr size_t Size = sizeof...(Ts);

      using Sizes = std::array<size_t, Size>;

      using Offsets = Sizes;

      LinearFormTupleAssemblyInput(
          const Sizes& sizes, const Offsets& offsets, const Tuple<Ts...>& ins)
        : m_sizes(sizes), m_offsets(offsets), m_ins(ins)
      {}

      const Sizes& getSizes() const
      {
        return m_sizes;
      }

      const Offsets& getOffsets() const
      {
        return m_offsets;
      }

      const Tuple<Ts...>& getTuple() const
      {
        return m_ins;
      }

    private:
      Sizes m_sizes;
      Offsets m_offsets;
      const Tuple<Ts...> m_ins;
  };
}

#endif
