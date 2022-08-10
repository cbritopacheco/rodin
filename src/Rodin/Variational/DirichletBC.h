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

#include "Function.h"
#include "ShapeFunction.h"

#include "Exceptions.h"

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
    * | Spaces supported     | H1                                            |
    * | Dimensions supported | 1D, 2D, 3D                                    |
    * | Continuous operator  | @f$ u = g \text{ on } \Gamma_D@f$             |
    * | @f$ g @f$            | ScalarFunction                                |
    */
   template <class Trait>
   class DirichletBC<TrialFunction<H1<Trait>>> : public FormLanguage::Base
   {
      public:
         DirichletBC(const TrialFunction<H1<Trait>>& u, const FunctionBase& v)
            : m_u(u), m_value(v.copy())
         {
            if (v.getRangeType() != RangeType::Scalar && v.getRangeType() != RangeType::Vector)
            {
               UnexpectedRangeTypeException(
                     {RangeType::Scalar, RangeType::Vector}, v.getRangeType()).raise();
            }
         }

         DirichletBC(const DirichletBC& other)
            : FormLanguage::Base(other),
              m_u(other.m_u),
              m_value(other.m_value->copy()),
              m_essBdr(other.m_essBdr)
         {}

         DirichletBC(DirichletBC&& other)
            : FormLanguage::Base(std::move(other)),
              m_u(other.m_u),
              m_value(std::move(other.m_value)),
              m_essBdr(std::move(other.m_essBdr))
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

         /**
          * @returns Returns reference to the value of the boundary condition
          * at the boundary
          */
         const FunctionBase& getValue() const
         {
            assert(m_value);
            return *m_value;
         }

         const TrialFunction<H1<Trait>>& getTrialFunction() const
         {
            return m_u;
         }

         /**
          * @returns Boundary attribute where the boundary condition is
          * imposed.
          */
         const std::set<int>& getBoundaryAttributes() const
         {
            return m_essBdr;
         }

         DirichletBC* copy() const noexcept override
         {
            return new DirichletBC(*this);
         }
      private:
         const TrialFunction<H1<Trait>>& m_u;
         std::unique_ptr<FunctionBase> m_value;
         std::set<int> m_essBdr;
   };
   template <class Trait>
   DirichletBC(const TrialFunction<H1<Trait>>&, const FunctionBase&)
      -> DirichletBC<TrialFunction<H1<Trait>>>;

   template <class Trait>
   class DirichletBC<Component<TrialFunction<H1<Trait>>>>
      : public FormLanguage::Base
   {
      public:
         DirichletBC(const Component<TrialFunction<H1<Trait>>>& ux, const FunctionBase& v)
            : m_ux(ux), m_value(v.copy())
         {
            if (v.getRangeType() != RangeType::Scalar)
               UnexpectedRangeTypeException(RangeType::Scalar, v.getRangeType());
         }

         DirichletBC(const DirichletBC& other)
            :  FormLanguage::Base(other),
               m_ux(other.m_ux),
               m_value(other.m_value->copy()),
               m_essBdr(other.m_essBdr)
         {}

         DirichletBC(DirichletBC&& other)
            :  FormLanguage::Base(std::move(other)),
               m_ux(std::move(other.m_ux)),
               m_value(std::move(other.m_value)),
               m_essBdr(std::move(other.m_essBdr))
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

         /**
          * @returns Returns reference to the value of the boundary condition
          * at the boundary
          */
         const FunctionBase& getValue() const
         {
            assert(m_value);
            return *m_value;
         }

         const Component<TrialFunction<H1<Trait>>>& getComponent() const
         {
            return m_ux;
         }

         const std::set<int>& getBoundaryAttributes() const
         {
            return m_essBdr;
         }

         DirichletBC* copy() const noexcept override
         {
            return new DirichletBC(*this);
         }
      private:
         Component<TrialFunction<H1<Trait>>> m_ux;
         std::unique_ptr<FunctionBase> m_value;
         std::set<int> m_essBdr;
   };
   template <class Trait>
   DirichletBC(const Component<TrialFunction<H1<Trait>>>&, const FunctionBase&)
      -> DirichletBC<Component<TrialFunction<H1<Trait>>>>;
}

#endif
