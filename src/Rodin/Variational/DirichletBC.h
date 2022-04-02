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
#include "ScalarFunction.h"
#include "VectorFunction.h"

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
    * | @f$ g @f$            | ScalarFunction                             |
    */
   template <class FES, class Value>
   class DirichletBC<TrialFunction<FES>, Value> : public FormLanguage::Base
   {
      static_assert(
            std::is_base_of_v<ScalarFunctionBase, Value> ||
            std::is_base_of_v<VectorFunctionBase, Value>,
            "Value must be derived from either ScalarFunctionBase or VectorFunctionBase");
      public:
         DirichletBC(const TrialFunction<FES>& u, const Value& v)
            : m_u(u), m_value(v.copy())
         {
            if constexpr (std::is_base_of_v<ScalarFunctionBase, Value>)
            {
               assert(
                     u.getFiniteElementSpace().getVectorDimension() == 1);
            }
            else if constexpr (std::is_base_of_v<VectorFunctionBase, Value>)
            {
               assert(
                     u.getFiniteElementSpace().getVectorDimension() == v.getDimension());
            }
         }

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
               m_value(other.m_value->copy()),
               m_essBdr(other.m_essBdr)
         {}

         DirichletBC(DirichletBC&& other)
            :  m_u(other.m_u),
               m_value(std::move(other.m_value)),
               m_essBdr(std::move(other.m_essBdr))
         {}

         /**
          * @returns Returns reference to the value of the boundary condition
          * at the boundary
          */
         const Value& getValue() const
         {
            assert(m_value);
            return *m_value;
         }

         const TrialFunction<FES>& getTrialFunction() const
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
         const TrialFunction<FES>& m_u;
         std::unique_ptr<Value> m_value;
         std::set<int> m_essBdr;
   };
   template <class FES>
   DirichletBC(const TrialFunction<FES>&, const ScalarFunctionBase&)
      -> DirichletBC<TrialFunction<FES>, ScalarFunctionBase>;
   template <class FES>
   DirichletBC(const TrialFunction<FES>&, const VectorFunctionBase&)
      -> DirichletBC<TrialFunction<FES>, VectorFunctionBase>;

   template <class FES>
   class DirichletBC<Component<TrialFunction<FES>>, ScalarFunctionBase> : public FormLanguage::Base
   {
      public:
         DirichletBC(const Component<TrialFunction<FES>>& ux, const ScalarFunctionBase& v)
            : m_ux(ux), m_value(v.copy())
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
            :  m_ux(other.m_ux),
               m_value(other.m_value->copy()),
               m_essBdr(other.m_essBdr)
         {}

         DirichletBC(DirichletBC&& other)
            :  m_ux(other.m_u),
               m_value(std::move(other.m_value)),
               m_essBdr(std::move(other.m_essBdr))
         {}

         /**
          * @returns Returns reference to the value of the boundary condition
          * at the boundary
          */
         const ScalarFunctionBase& getValue() const
         {
            assert(m_value);
            return *m_value;
         }

         const Component<TrialFunction<FES>>& getComponent() const
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
         Component<TrialFunction<FES>> m_ux;
         std::unique_ptr<ScalarFunctionBase> m_value;
         std::set<int> m_essBdr;
   };
   template <class FES>
   DirichletBC(const Component<TrialFunction<FES>>&, const ScalarFunctionBase&)
      -> DirichletBC<Component<TrialFunction<FES>>, ScalarFunctionBase>;
}

#endif
