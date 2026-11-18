/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_VARIATIONAL_BILINEARFORM_H
#define RODIN_KOKKOS_VARIATIONAL_BILINEARFORM_H

/**
 * @file
 * @brief Backend algebra deduction for bilinear forms.
 */

#include "Rodin/Variational/BilinearForm.h"
#include "Rodin/Variational/TrialFunction.h"
#include "Rodin/Variational/TestFunction.h"
#include "Rodin/Kokkos/Math/SparseMatrix.h"
#include "DOFs.h"

namespace Rodin::Variational
{
  /// @param u Trial function. @param v Test function on a Kokkos space.
  template <class Solution, ::Rodin::Kokkos::Variational::Space TrialFES,
    ::Rodin::Kokkos::Variational::Space TestFES>
  BilinearForm(const TrialFunction<Solution, TrialFES>& u, const TestFunction<TestFES>& v)
    -> BilinearForm<Solution, TrialFES, TestFES,
      ::Rodin::Kokkos::Math::SparseMatrix<
        typename FormLanguage::Mult<typename FormLanguage::Traits<TrialFES>::ScalarType,
          typename FormLanguage::Traits<TestFES>::ScalarType>::Type,
        typename TrialFES::MeshType::ContextType::BackendType>>;

}

#endif
