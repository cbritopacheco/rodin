/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_DERIVATIVE_H
#define RODIN_VARIATIONAL_DERIVATIVE_H

#include <cstddef>
#include <mfem.hpp>

#include "ForwardDecls.h"
#include "ScalarCoefficient.h"

namespace Rodin::Variational
{
   /**
    * @brief Represents the partial derivative @f$ \frac{\partial u_i}{\partial
    * x_j} @f$ of the function @f$ u @f$.
    *
    * Denote the value of @f$ u \colon \mathbb{R}^n \rightarrow \mathbb{R}^m
    * @f$ at @f$ x = (x_1, \ldots, x_n) @f$ by
    * @f[
    *    u(x) = (u_1(x), \ldots, u_m(x))
    * @f]
    * This class represents the partial derivative with respect to the i-th
    * component and the j-th variable
    * @f[
    *    \dfrac{\partial u_i}{\partial x_j}
    * @f]
    *
    * @tparam DirectionIndex Index of the direction
    * @tparam ComponentIndex Index of the function component
    */
   template <int DirectionIndex, int ComponentIndex>
   class Derivative : public ScalarCoefficientBase
   {
      static_assert(ComponentIndex >= 1 && DirectionIndex >= 1,
            "ComponentIndex and VariableIndex must both be greater than 1.");
      public:
         /**
          * @brief Constructs the Derivative of an @f$ H^1 @f$ function
          * @f$ u @f$.
          * @param[in] u Grid function to be differentiated.
          */
         Derivative(GridFunctionBase& u);

         void buildMFEMCoefficient() override;

         mfem::Coefficient& getMFEMCoefficient() override;

         Derivative* copy() const noexcept override
         {
            return new Derivative(*this);
         }

      private:
         GridFunctionBase& m_u;
         std::optional<mfem::GridFunction> m_mfemGridFunction;
         std::optional<mfem::GridFunctionCoefficient> m_mfemCoefficient;
   };

   /**
    * @brief Represents the partial derivative in the @f$ x @f$ direction of
    * the j-th component.
    * @tparam ComponentIndex Index of the function component
    *
    * In other words, this function returns the Derivative object representing
    * the partial derivative:
    * @f[
    *    \dfrac{\partial u_j}{\partial x_1}
    * @f]
    *
    * @returns Derivative object with respect to the j-th component in the
    * 1st direction.
    */
   template <int ComponentIndex = 1>
   Derivative<1, ComponentIndex> Dx(GridFunctionBase& u)
   {
      return Derivative<1, ComponentIndex>(u);
   }

   /**
    * @brief Represents the partial derivative in the @f$ y @f$ direction of
    * the j-th component.
    * @tparam ComponentIndex Index of the function component
    *
    * In other words, this function returns the Derivative object representing
    * the partial derivative:
    * @f[
    *    \dfrac{\partial u_j}{\partial x_2}
    * @f]
    *
    * @returns Derivative object with respect to the j-th component in the
    * 2nd direction.
    */
   template <int ComponentIndex = 1>
   Derivative<2, ComponentIndex> Dy(GridFunctionBase& u)
   {
      return Derivative<2, ComponentIndex>(u);
   }

   /**
    * @brief Represents the partial derivative in the @f$ x @f$ direction of
    * the j-th component.
    * @tparam ComponentIndex Index of the function component
    *
    * In other words, this function returns the Derivative object representing
    * the partial derivative:
    * @f[
    *    \dfrac{\partial u_j}{\partial x_3}
    * @f]
    *
    * @returns Derivative object with respect to the j-th component in the
    * 3rd direction.
    */
   template <int ComponentIndex = 1>
   Derivative<3, ComponentIndex> Dz(GridFunctionBase& u)
   {
      return Derivative<3, ComponentIndex>(u);
   }
}

#include "Derivative.hpp"

#endif
