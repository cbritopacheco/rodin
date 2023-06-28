/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_BILINEARFORM_HPP
#define RODIN_VARIATIONAL_BILINEARFORM_HPP

#include <cassert>

#include "Rodin/Alert.h"

#include "Assembly/AssemblyBase.h"

#include "BilinearForm.h"
#include "BilinearFormIntegrator.h"
#include "FiniteElementSpace.h"


namespace Rodin::Variational
{
  template <class TrialFES, class TestFES>
  void
  BilinearForm<TrialFES, TestFES, Context::Serial, Math::SparseMatrix>::assemble()
  {
     assert(&getTrialFunction().getFiniteElementSpace().getMesh() ==
           &getTestFunction().getFiniteElementSpace().getMesh());
     const auto& trialFES = getTrialFunction().getFiniteElementSpace();
     const auto& testFES = getTestFunction().getFiniteElementSpace();
     const auto& mesh = getTrialFunction().getFiniteElementSpace().getMesh();
     m_operator = getAssembly().execute({mesh, trialFES, testFES, getIntegrators()});
  }

  template <class TrialFES, class TestFES>
  void
  BilinearForm<TrialFES, TestFES, Context::Serial, std::vector<Eigen::Triplet<Scalar>>>::assemble()
  {
     assert(&getTrialFunction().getFiniteElementSpace().getMesh() ==
           &getTestFunction().getFiniteElementSpace().getMesh());
     const auto& trialFES = getTrialFunction().getFiniteElementSpace();
     const auto& testFES = getTestFunction().getFiniteElementSpace();
     const auto& mesh = getTrialFunction().getFiniteElementSpace().getMesh();
     m_operator = getAssembly().execute({mesh, trialFES, testFES, getIntegrators()});
  }
}

#endif
