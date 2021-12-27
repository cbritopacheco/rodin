#ifndef RODIN_VARIATIONAL_BOUNDARYCONDITION_H
#define RODIN_VARIATIONAL_BOUNDARYCONDITION_H

#include <variant>

#include "FormLanguage/Base.h"

#include "ForwardDecls.h"

#include "ScalarCoefficient.h"
#include "VectorCoefficient.h"

namespace Rodin::Variational
{
   class BoundaryConditionBase : public FormLanguage::Base
   {
      public:
         virtual int getBoundaryAttribute() const = 0;
         virtual FormLanguage::Base& getValue() = 0;
         virtual void imposeOn(ProblemBase& pb) = 0;
         virtual BoundaryConditionBase* copy() const noexcept override = 0;
   };

   BoundaryCondition(int, const ScalarCoefficientBase&)
      -> BoundaryCondition<ScalarCoefficientBase>;

   template <>
   class BoundaryCondition<ScalarCoefficientBase> : public BoundaryConditionBase
   {
      public:
         BoundaryCondition(int bdrAttr, const ScalarCoefficientBase& value)
            : m_bdrAttr(bdrAttr), m_value(value.copy())
         {}

         BoundaryCondition(const BoundaryCondition& other)
            : m_bdrAttr(other.m_bdrAttr), m_value(other.m_value->copy())
         {}

         virtual ~BoundaryCondition() = default;

         /**
          * @brief Gets the boundary attribute where the boundary condition
          * will be imposed
          */
         virtual int getBoundaryAttribute() const override
         {
            return m_bdrAttr;
         }

         virtual ScalarCoefficientBase& getValue() override
         {
            return *m_value;
         }

         virtual BoundaryCondition* copy() const noexcept override = 0;

      private:
         int m_bdrAttr;
         std::unique_ptr<ScalarCoefficientBase> m_value;
   };

   BoundaryCondition(int, const VectorCoefficientBase&)
      -> BoundaryCondition<VectorCoefficientBase>;

   template <>
   class BoundaryCondition<VectorCoefficientBase> : public BoundaryConditionBase
   {
      public:
         BoundaryCondition(int bdrAttr, const VectorCoefficientBase& value)
            : m_bdrAttr(bdrAttr), m_value(value.copy())
         {}

         BoundaryCondition(const BoundaryCondition& other)
            : m_bdrAttr(other.m_bdrAttr), m_value(other.m_value->copy())
         {}

         virtual ~BoundaryCondition() = default;

         /**
          * @brief Gets the boundary attribute where the boundary condition
          * will be imposed
          */
         virtual int getBoundaryAttribute() const override
         {
            return m_bdrAttr;
         }

         virtual VectorCoefficientBase& getValue() override
         {
            return *m_value;
         }

         virtual BoundaryCondition* copy() const noexcept override = 0;

      private:
         int m_bdrAttr;
         std::unique_ptr<VectorCoefficientBase> m_value;
   };
}

#endif

