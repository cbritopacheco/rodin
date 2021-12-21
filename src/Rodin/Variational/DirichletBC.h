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

#include "BoundaryCondition.h"

namespace Rodin::Variational
{
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
    * | @f$ g @f$            | ScalarCoefficient if @f$ u @f$ is scalar valued, otherwise VectorCoefficient with the same dimension as @f$ u @f$.
    * |                      |                                               |
    *
    * @see @ref examples-variational-poisson
    *
    */
   class DirichletBC : public BoundaryConditionBase
   {
      public:
         /**
          * @brief Constructs a Dirichlet boundary condition on the part of the
          * boundary specified by the boundary attribute.
          *
          * @param[in] bdrAttr Attribute corresponding to the segment
          * @f$ \Gamma_D @f$ where the Dirichlet boundary condition is imposed.
          *
          * @param[in] value Scalar value of the trial function @f$ u @f$ on
          * the boundary @f$ \Gamma_D @f$.
          */
         DirichletBC(int bdrAttr, const ScalarCoefficientBase& value);

         /**
          * @brief Constructs a Dirichlet boundary condition on the part of the
          * boundary specified by the boundary attribute.
          *
          * @param[in] bdrAttr Attribute corresponding to the segment
          * @f$ \Gamma_D @f$ where the Dirichlet boundary condition is imposed.
          *
          * @param[in] value Vector value of the trial function @f$ u @f$ on
          * the boundary @f$ \Gamma_D @f$.
          */
         DirichletBC(int bdrAttr, const VectorCoefficientBase& value);

         /**
          * @internal
          */
         DirichletBC(const DirichletBC& other);

         virtual ~DirichletBC() = default;

         int getBoundaryAttribute() const override;

         void imposeOn(ProblemBase& pb) override;

         DirichletBC* copy() const noexcept override
         {
            return new DirichletBC(*this);
         }

      private:
         int m_bdrAttr;
         mfem::Array<int> m_essBdr;
         std::variant<
            std::unique_ptr<ScalarCoefficientBase>,
            std::unique_ptr<VectorCoefficientBase>
            > m_value;
   };
}

#endif
