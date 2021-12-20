#ifndef RODIN_VARIATIONAL_BOUNDARYCONDITION_H
#define RODIN_VARIATIONAL_BOUNDARYCONDITION_H

#include "FormLanguage/Base.h"

#include "ForwardDecls.h"

namespace Rodin::Variational
{
   class BoundaryConditionBase : public FormLanguage::Base
   {
      public:
         /**
          * @brief Gets the boundary attribute where the boundary condition
          * will be imposed
          */
         virtual int getBoundaryAttribute() const = 0;

         /**
          * @brief Imposes the boundary condition on the given problem.
          * @param[in, out] pb Problem where the boundary condition will be
          * imposed.
          */
         virtual void imposeOn(ProblemBase& pb) const = 0;

         virtual BoundaryConditionBase* copy() const noexcept override = 0;
   };
}

#endif

