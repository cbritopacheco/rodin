/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef RODIN_EXTERNAL_MMG_CAST_HPP
#define RODIN_EXTERNAL_MMG_CAST_HPP

#include "Cast.h"

#include "ScalarSolution2D.h"

namespace Rodin
{
  // ---- Cast<MMG::ScalarSolution2D, Rodin::GridFunction> -------------------
  template <class FEC, bool HasMesh>
  Cast<External::MMG::ScalarSolution2D<HasMesh>, Variational::GridFunction<FEC>>
  ::Cast(Variational::FiniteElementSpace<FEC>& fes)
    : m_fes(fes)
  {}

  template <class FEC, bool HasMesh>
  Variational::GridFunction<FEC>
  Cast<External::MMG::ScalarSolution2D<HasMesh>, Variational::GridFunction<FEC>>
  ::cast(const External::MMG::ScalarSolution2D<HasMesh>& sol)
  const
  {
    MMG5_pSol mmgSol = sol.getHandle();
    assert(mmgSol->type == MMG5_Scalar);

    Variational::GridFunction<FEC> res(m_fes);
    std::unique_ptr<double[]> data(new double[mmgSol->np]);
    // MMG5_pSol is 1 indexed. We must start at m + 1 and finish at m + np + 1.
    std::copy(mmgSol->m + 1, mmgSol->m + mmgSol->np + 1, data.get());
    res.setData(std::move(data), mmgSol->np);

    return res;
  }

  // ---- Cast<Rodin::GridFunction, MMG::ScalarSolution2D> -------------------
  // template <class FEC>
  // External::MMG::ScalarSolution2D
  // Cast<Variational::GridFunction<FEC>, External::MMG::ScalarSolution2D>
  // ::cast(const Variational::GridFunction<FEC>& gf)
  // const
  // {
  // }
}

#endif

