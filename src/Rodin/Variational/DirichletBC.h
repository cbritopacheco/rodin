/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_DIRICHLETBC_H
#define RODIN_VARIATIONAL_DIRICHLETBC_H

#include <variant>

#include "ForwardDecls.h"

#include "ScalarCoefficient.h"
#include "VectorCoefficient.h"

#include "BoundaryCondition.h"

namespace Rodin::Variational
{
   template <class T>
   class DirichletBC;

   DirichletBC(int, const ScalarCoefficientBase&)
      -> DirichletBC<ScalarCoefficientBase>;

   /**
    * @brief Represents a Dirichlet boundary condition.
    *
    * When utilized in a Problem construction, it will impose the Dirichlet
    * condition
    * @f[
    *    u = g \text{ on } \Gamma_D
    * @f]
    * on the segment of the boundary @f$ \Gamma_D \subset \partial \Omega @f$
    * specified by the boundary attribute.
    *
    * | Detail               | Description                                   |
    * |----------------------|-----------------------------------------------|
    * | Spaces supported     | L2, H1                                        |
    * | Dimensions supported | 1D, 2D, 3D                                    |
    * | Continuous operator  | @f$ u = g \text{ on } \Gamma_D@f$             |
    * | @f$ g @f$            | ScalarCoefficient                             |
    *
    */
   template <>
   class DirichletBC<ScalarCoefficientBase>
      : public BoundaryCondition<ScalarCoefficientBase>
   {
      public:
         /**
          * @brief Constructs a DirichletBC with a scalar valued coefficient
          * @param[in] bdrAttr Attribute where the condition will be imposed.
          * @param[in] v Derived instance of ScalarCoefficientBase
          */
         DirichletBC(int bdrAtr, const ScalarCoefficientBase& v)
            : BoundaryCondition<ScalarCoefficientBase>(bdrAtr, v)
         {}

         void imposeOn(ProblemBase& pb) override;

         DirichletBC* copy() const noexcept override
         {
            return new DirichletBC(*this);
         }
      private:
         mfem::Array<int> m_essBdr;
   };

   DirichletBC(int, const VectorCoefficientBase&)
      -> DirichletBC<VectorCoefficientBase>;

   /**
    * @brief Represents a Dirichlet boundary condition.
    *
    * When utilized in a Problem construction, it will impose the Dirichlet
    * condition
    * @f[
    *    u = g \text{ on } \Gamma_D
    * @f]
    * on the segment of the boundary @f$ \Gamma_D \subset \partial \Omega @f$
    * specified by the boundary attribute.
    *
    * | Detail               | Description                                   |
    * |----------------------|-----------------------------------------------|
    * | Spaces supported     | L2, H1                                        |
    * | Dimensions supported | 1D, 2D, 3D                                    |
    * | Continuous operator  | @f$ u = g \text{ on } \Gamma_D@f$             |
    * | @f$ g @f$            | VectorCoefficient                             |
    *
    * @see @ref examples-variational-poisson
    *
    */
   template <>
   class DirichletBC<VectorCoefficientBase>
      : public BoundaryCondition<VectorCoefficientBase>
   {
      public:
         /**
          * @brief Constructs a DirichletBC with a vector valued coefficient.
          * @param[in] bdrAttr Attribute where the condition will be imposed.
          * @param[in] v Derived instance of VectorCoefficientBase
          */
         DirichletBC(int bdrAtr, const VectorCoefficientBase& v)
            : BoundaryCondition<VectorCoefficientBase>(bdrAtr, v)
         {}

         void imposeOn(ProblemBase& pb) override;

         DirichletBC* copy() const noexcept override
         {
            return new DirichletBC(*this);
         }
      private:
         mfem::Array<int> m_essBdr;
   };
}

#endif
