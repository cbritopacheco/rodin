#ifndef RODIN_EXTERNAL_MMG_CAST_HPP
#define RODIN_EXTERNAL_MMG_CAST_HPP

#include "Cast.h"

#include "ScalarSolution2D.h"

namespace Rodin::Cast
{
  // ---- <MMG::ScalarSolution2D, Rodin::GridFunction> -----------------------
  template <class FEC>
  Cast<External::MMG::ScalarSolution2D, Variational::GridFunction<FEC>>
  ::Cast(Variational::FiniteElementSpace<FEC>& fes)
    : m_fes(fes)
  {}

  template <class FEC>
  Variational::GridFunction<FEC>
  Cast<External::MMG::ScalarSolution2D, Variational::GridFunction<FEC>>
  ::cast(const External::MMG::ScalarSolution2D& sol)
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
}

#endif

