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
  template <class OperatorType, class TrialFES, class TestFES>
  class Input<Variational::BilinearForm<TrialFES, TestFES, OperatorType>>
  {
    public:
      Input(Variational::BilinearForm<TrialFES, TestFES, OperatorType>& bf)
        : m_mesh(bf.getTrialFunction().getFiniteElementSpace().getMesh()),
          m_trialFES(bf.getTrialFunction().getFiniteElementSpace()),
          m_testFES(bf.getTestFunction().getFiniteElementSpace()),
          m_lbfis(bf.getLocalIntegrators()), m_gbfis(bf.getGlobalIntegrators())
      {}

      Input(const Geometry::MeshBase& mesh, const TrialFES& trialFES, const TestFES& testFES,
          FormLanguage::List<Variational::LocalBilinearFormIntegratorBase>& lbfis,
          FormLanguage::List<Variational::GlobalBilinearFormIntegratorBase>& gbfis)
        : m_mesh(mesh), m_trialFES(trialFES), m_testFES(testFES), m_lbfis(lbfis), m_gbfis(gbfis)
      {}

      const Geometry::MeshBase& getMesh() const
      {
        return m_mesh.get();
      }

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
      std::reference_wrapper<const Geometry::MeshBase> m_mesh;
      std::reference_wrapper<const TrialFES> m_trialFES;
      std::reference_wrapper<const TestFES> m_testFES;
      std::reference_wrapper<FormLanguage::List<Variational::LocalBilinearFormIntegratorBase>> m_lbfis;
      std::reference_wrapper<FormLanguage::List<Variational::GlobalBilinearFormIntegratorBase>> m_gbfis;
  };

  template <class OperatorType, class TrialFES, class TestFES>
  Input(Variational::BilinearForm<TrialFES, TestFES, OperatorType>&)
    -> Input<Variational::BilinearForm<TrialFES, TestFES, OperatorType>>;

  template <class VectorType, class FES>
  class Input<Variational::LinearForm<FES, VectorType>>
  {
    public:
      Input(Variational::LinearForm<FES, VectorType>& lf)
        : m_mesh(lf.getTestFunction().getFiniteElementSpace().getMesh()),
          m_fes(lf.getFiniteElementSpace()), m_lfis(lf.getIntegrators())
      {}

      Input(const Geometry::MeshBase& mesh, const FES& fes,
          FormLanguage::List<Variational::LinearFormIntegratorBase>& lfis)
        : m_mesh(mesh), m_fes(fes), m_lfis(lfis)
      {}

      const Geometry::MeshBase& getMesh() const
      {
        return m_mesh.get();
      }

      const FES& getFES() const
      {
        return m_fes.get();
      }

      FormLanguage::List<Variational::LinearFormIntegratorBase>& getLFIs() const
      {
        return m_lfis.get();
      }

    private:
      std::reference_wrapper<const Geometry::MeshBase> m_mesh;
      std::reference_wrapper<const FES> m_fes;
      std::reference_wrapper<FormLanguage::List<Variational::LinearFormIntegratorBase>> m_lfis;
  };

  template <class VectorType, class FES>
  Input(Variational::LinearForm<FES, VectorType>&)
    -> Input<Variational::LinearForm<FES, VectorType>>;

  template <class OperatorType, class ... TrialFES, class ... TestFES>
  class Input<Tuple<Variational::BilinearForm<TrialFES, TestFES, OperatorType>...>>
    : public Tuple<Input<Variational::BilinearForm<TrialFES, TestFES, OperatorType>>...>
  {
    public:
      using Parent = Tuple<Input<Variational::BilinearForm<TrialFES, TestFES, OperatorType>>...>;
      using Parent::Parent;
  };
}

#endif
