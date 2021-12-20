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
   class NeumannBC : public BoundaryConditionBase
   {
      public:
         /**
          * @brief Constructs a Neumann boundary condition on the part of the
          * boundary specified by the boundary attribute.
          *
          * @param[in] bdrAttr Attribute corresponding to the segment
          * @f$ \Gamma_D @f$ where the Neumann boundary condition is imposed.
          *
          * @param[in] value Scalar value of the trial function @f$ u @f$ on
          * the boundary @f$ \Gamma_D @f$.
          */
         NeumannBC(int bdrAttr, const ScalarCoefficientBase& value);

         /**
          * @brief Constructs a Neumann boundary condition on the part of the
          * boundary specified by the boundary attribute.
          *
          * @param[in] bdrAttr Attribute corresponding to the segment
          * @f$ \Gamma_D @f$ where the Neumann boundary condition is imposed.
          *
          * @param[in] value Vector value of the trial function @f$ u @f$ on
          * the boundary @f$ \Gamma_D @f$.
          */
         NeumannBC(int bdrAttr, const VectorCoefficientBase& value);

         NeumannBC(const NeumannBC& other);

         int getBoundaryAttribute() const override;

         NeumannBC* copy() const noexcept override
         {
            return new NeumannBC(*this);
         }

         void imposeOn(ProblemBase& pb) const override;

      private:
         int m_bdrAttr;
         std::variant<
            std::unique_ptr<ScalarCoefficientBase>,
            std::unique_ptr<VectorCoefficientBase>
            > m_value;
   };
}

#endif
