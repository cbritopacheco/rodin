/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Class @ref NeumannBC
 */
#ifndef RODIN_VARIATIONAL_NEUMANNBC_H
#define RODIN_VARIATIONAL_NEUMANNBC_H

#include <variant>
#include <functional>

#include "ForwardDecls.h"

#include "ScalarCoefficient.h"
#include "BoundaryCondition.h"

namespace Rodin::Variational
{
   /**
    * @brief Represents a Neumann boundary condition.
    *
    * The usage of this class assumes the proper coefficients are present
    * in the boundary conditions of the equation. In particular, when utilized
    * in a @ref Problem construction, it will add a term of the form:
    * @f[
    *    \int_{\Gamma_N} g v \ dx
    * @f]
    * to the right hand side of the variational formulation.
    *
    * | Detail                | Description                                  |
    * |-----------------------|----------------------------------------------|
    * |  Spaces supported     | L2, H1                                       |
    * |  Dimensions supported | 1D, 2D, 3D                                   |
    *
    */
   template <class T>
   class NeumannBC;

   NeumannBC(int, const ScalarCoefficientBase&)
      -> NeumannBC<ScalarCoefficientBase>;

   template <>
   class NeumannBC<ScalarCoefficientBase>
      : public BoundaryCondition<ScalarCoefficientBase>
   {
      public:
         NeumannBC(int bdrAtr, const ScalarCoefficientBase& v)
            : BoundaryCondition<ScalarCoefficientBase>(bdrAtr, v)
         {}

         void imposeOn(ProblemBase& pb) override;

         NeumannBC* copy() const noexcept override
         {
            return new NeumannBC(*this);
         }
      private:
         mfem::Array<int> m_nbcBdr;
   };

   NeumannBC(int, const VectorCoefficientBase&)
      -> NeumannBC<VectorCoefficientBase>;

   template <>
   class NeumannBC<VectorCoefficientBase>
      : public BoundaryCondition<VectorCoefficientBase>
   {
      public:
         NeumannBC(int bdrAtr, const VectorCoefficientBase& v)
            : BoundaryCondition<VectorCoefficientBase>(bdrAtr, v)
         {}

         void imposeOn(ProblemBase& pb) override;

         NeumannBC* copy() const noexcept override
         {
            return new NeumannBC(*this);
         }
      private:
         mfem::Array<int> m_nbcBdr;
   };
}

#endif
