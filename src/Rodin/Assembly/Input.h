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
#include "Rodin/FormLanguage/Base.h"
#include "Rodin/Variational/ForwardDecls.h"

#include "ForwardDecls.h"

namespace Rodin::Assembly
{
  template <class TrialFES, class TestFES>
  class BilinearFormAssemblyInput
  {
    public:
      using TrialFESType = TrialFES;

      using TestFESType = TestFES;

      using TrialFESScalarType  = typename FormLanguage::Traits<TrialFESType>::ScalarType;

      using TestFESScalarType   = typename FormLanguage::Traits<TestFESType>::ScalarType;

      using ScalarType = decltype(
          std::declval<TrialFESScalarType>() * std::declval<TestFESScalarType>());

      using LocalBilinearFormIntegratorBaseType =
        Variational::LocalBilinearFormIntegratorBase<ScalarType>;

      using LocalBilinearFormIntegratorBaseListType =
        FormLanguage::List<LocalBilinearFormIntegratorBaseType>;

      using GlobalBilinearFormIntegratorBaseType =
        Variational::GlobalBilinearFormIntegratorBase<ScalarType>;

      using GlobalBilinearFormIntegratorBaseListType =
        FormLanguage::List<GlobalBilinearFormIntegratorBaseType>;

      BilinearFormAssemblyInput(
          const TrialFES& trialFES, const TestFES& testFES,
          LocalBilinearFormIntegratorBaseListType& lbfis,
          GlobalBilinearFormIntegratorBaseListType& gbfis)
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

      LocalBilinearFormIntegratorBaseListType& getLocalBFIs() const
      {
        return m_lbfis.get();
      }

      GlobalBilinearFormIntegratorBaseListType& getGlobalBFIs() const
      {
        return m_gbfis.get();
      }

    private:
      std::reference_wrapper<const TrialFES>  m_trialFES;
      std::reference_wrapper<const TestFES>   m_testFES;
      std::reference_wrapper<LocalBilinearFormIntegratorBaseListType>   m_lbfis;
      std::reference_wrapper<GlobalBilinearFormIntegratorBaseListType>  m_gbfis;
  };

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

  template <class FES>
  class LinearFormAssemblyInput
  {
    public:
      using FESType = FES;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using LinearFormIntegratorBaseList =
        FormLanguage::List<Variational::LinearFormIntegratorBase<ScalarType>>;

      LinearFormAssemblyInput(const FES& fes, LinearFormIntegratorBaseList& lfis)
        : m_fes(fes), m_lfis(lfis)
      {}

      const FES& getFES() const
      {
        return m_fes.get();
      }

      LinearFormIntegratorBaseList& getLFIs() const
      {
        return m_lfis.get();
      }

    private:
      std::reference_wrapper<const FES> m_fes;
      std::reference_wrapper<LinearFormIntegratorBaseList> m_lfis;
  };

  template <class ... Ts>
  class BilinearFormTupleAssemblyInput
  {
    public:
      static constexpr size_t Size = sizeof...(Ts);

      using Offsets = std::array<Pair<size_t, size_t>, sizeof...(Ts)>;

      BilinearFormTupleAssemblyInput(
          size_t rows, size_t cols,
          const Offsets& offsets,
          const Tuple<Ts...>& ins)
        : m_rows(rows), m_cols(cols), m_offsets(offsets), m_ins(ins)
      {}

      inline
      size_t getRows() const
      {
        return m_rows;
      }

      inline
      size_t getColumns() const
      {
        return m_cols;
      }

      inline
      const Offsets& getOffsets() const
      {
        return m_offsets;
      }

      inline
      const Tuple<Ts...>& getTuple() const
      {
        return m_ins;
      }

    private:
      size_t m_rows, m_cols;
      Offsets m_offsets;
      const Tuple<Ts...> m_ins;
  };

  template <class ... Ts>
  class LinearFormTupleAssemblyInput
  {
    public:
      static constexpr size_t Size = sizeof...(Ts);

      using Offsets = std::array<size_t, sizeof...(Ts)>;

      LinearFormTupleAssemblyInput(
          size_t size, const Offsets& offsets, const Tuple<Ts...>& ins)
        : m_size(size), m_offsets(offsets), m_ins(ins)
      {}

      inline
      size_t getSize() const
      {
        return m_size;
      }

      inline
      const Offsets& getOffsets() const
      {
        return m_offsets;
      }

      inline
      const Tuple<Ts...>& getTuple() const
      {
        return m_ins;
      }

    private:
      size_t m_size;
      Offsets m_offsets;
      const Tuple<Ts...> m_ins;
  };
}

#endif
