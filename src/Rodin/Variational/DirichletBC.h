/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_DIRICHLETBC_H
#define RODIN_VARIATIONAL_DIRICHLETBC_H

#include <set>
#include <variant>

#include "Rodin/Utility.h"

#include "ForwardDecls.h"

#include "ShapeFunction.h"
#include "ScalarCoefficient.h"
#include "VectorCoefficient.h"

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
    * | @f$ g @f$            | ScalarCoefficient                             |
    */
   template <>
   class DirichletBC<TrialFunction<H1>>
   {
      public:
         DirichletBC(const TrialFunction<H1>& u, double v)
            : DirichletBC(u, ScalarCoefficient(v))
         {}

         DirichletBC(const TrialFunction<H1>& u, const ScalarCoefficientBase& v)
            :  m_u(u),
               m_value(std::unique_ptr<ScalarCoefficientBase>(v.copy()))
         {}

         DirichletBC(const TrialFunction<H1>& u, const VectorCoefficientBase& v)
            :  m_u(u),
               m_value(std::unique_ptr<VectorCoefficientBase>(v.copy()))
         {}

         DirichletBC& on(int bdrAtr)
         {
            return on(std::set<int>{bdrAtr});
         }

         DirichletBC& on(const std::set<int>& bdrAtr)
         {
            m_essBdr = bdrAtr;
            return *this;
         }

         DirichletBC(const DirichletBC& other)
            :  m_u(other.m_u),
               m_essBdr(other.m_essBdr),
               m_componentIdx(other.m_componentIdx)
         {
               std::visit(Utility::Overloaded{
                  [this](const std::unique_ptr<ScalarCoefficientBase>& v)
                  { m_value = std::unique_ptr<ScalarCoefficientBase>(v->copy()); },
                  [this](const std::unique_ptr<VectorCoefficientBase>& v)
                  { m_value = std::unique_ptr<VectorCoefficientBase>(v->copy()); },
                  }, other.m_value);
         }

         DirichletBC(DirichletBC&& other)
            :  m_u(other.m_u),
               m_essBdr(std::move(other.m_essBdr)),
               m_value(std::move(other.m_value)),
               m_componentIdx(std::move(other.m_componentIdx))
         {}

         /**
          * @returns Boundary attribute where the boundary condition is
          * imposed.
          */
         const std::set<int>& getBoundaryAttributes() const
         {
            return m_essBdr;
         }

         /**
          * @returns Returns reference to the value of the boundary condition
          * at the boundary
          */
         template <class CoefficientType>
         const CoefficientType& getValue() const
         {
            assert(std::holds_alternative<std::unique_ptr<CoefficientType>>(m_value));
            return *std::get<std::unique_ptr<CoefficientType>>(m_value);
         }

         DirichletBC* copy() const noexcept
         {
            return new DirichletBC(*this);
         }
      private:
         const TrialFunction<H1>& m_u;
         std::set<int> m_essBdr;
         std::variant<
            std::unique_ptr<ScalarCoefficientBase>,
            std::unique_ptr<VectorCoefficientBase>> m_value;
         std::optional<int> m_componentIdx;
   };
   DirichletBC(const TrialFunction<H1>&, double)
      -> DirichletBC<TrialFunction<H1>>;
   DirichletBC(const TrialFunction<H1>&, const ScalarCoefficientBase&)
      -> DirichletBC<TrialFunction<H1>>;
   DirichletBC(const TrialFunction<H1>&, const VectorCoefficientBase&)
      -> DirichletBC<TrialFunction<H1>>;
}

#endif
