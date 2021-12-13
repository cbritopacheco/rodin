/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_L2_H
#define RODIN_VARIATIONAL_L2_H

#include "FiniteElementSpace.h"

namespace Rodin::Variational
{
   /**
    * @brief Arbitrary order @f$ L^2(\Omega)^d @f$ conforming (discontinuous)
    * finite element space.
    *
    * Given some discretization @f$ \mathcal{T}_h @f$ (e.g. a triangulation)
    * of @f$ \Omega @f$, instances of this class will represent the space
    * @f[
    *    V_h := \left\{ v \in L^2(\Omega)^d \mid v_{|\tau} \in \mathcal{P},
    *    \ \forall \tau \in \mathcal{T}_h \right\}
    * @f]
    * where @f$ \mathcal{P} @f$ denotes a vector space of functions from
    * @f$ \tau @f$ to @f$ \mathbb{R}^d @f$.
    *
    */
   class L2 : public FiniteElementSpace<L2>
   {

   };
}

#endif
