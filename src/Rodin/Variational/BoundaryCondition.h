#ifndef RODIN_VARIATIONAL_BOUNDARYCONDITION_H
#define RODIN_VARIATIONAL_BOUNDARYCONDITION_H

#include <variant>

#include "FormLanguage/Base.h"

#include "ForwardDecls.h"

#include "ScalarCoefficient.h"
#include "VectorCoefficient.h"

namespace Rodin::Variational
{
   /**
    * @brief Abstract class for objects representing a boundary condition.
    */
   class BoundaryConditionBase : public FormLanguage::Base
   {
      public:
         /**
          * @returns Boundary attribute where the boundary condition is
          * imposed.
          */
         virtual int getBoundaryAttribute() const = 0;

         /**
          * @returns Returns reference to the value of the boundary condition
          * at the boundary
          */
         virtual const FormLanguage::Base& getValue() const = 0;

         /**
          * @brief Imposes the boundary condition on the given ProblemBase
          * instance.
          */
         virtual void imposeOn(ProblemBase& pb) = 0;

         virtual BoundaryConditionBase* copy() const noexcept override = 0;
   };

   BoundaryCondition(int, const ScalarCoefficientBase&)
      -> BoundaryCondition<ScalarCoefficientBase>;

   template <>
   class BoundaryCondition<ScalarCoefficientBase> : public BoundaryConditionBase
   {
      public:
         /**
          * Constructs a boundary condition
          */
         BoundaryCondition(int bdrAttr, const ScalarCoefficientBase& value)
            : m_bdrAttr(bdrAttr), m_value(value.copy())
         {}

         BoundaryCondition(const BoundaryCondition& other)
            : m_bdrAttr(other.m_bdrAttr), m_value(other.m_value->copy())
         {}

         virtual ~BoundaryCondition() = default;

         virtual int getBoundaryAttribute() const override
         {
            return m_bdrAttr;
         }

         virtual const ScalarCoefficientBase& getValue() const override
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

         virtual VectorCoefficientBase& getValue() const override
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

