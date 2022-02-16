#ifndef RODIN_VARIATIONAL_TRIALFUNCTION_H
#define RODIN_VARIATIONAL_TRIALFUNCTION_H

#include <memory>

#include "Rodin/Alert.h"
#include "ForwardDecls.h"
#include "FiniteElementSpace.h"

namespace Rodin::Variational
{
   class TrialFunctionExpr
   {
      public:
         virtual int getDimension() const = 0;
   };

   template <class FEC>
   class TrialFunction : public TrialFunctionExpr
   {
      public:
         TrialFunction(FiniteElementSpace<FEC>& fes)
            : m_fes(fes)
         {}

         const FiniteElementSpace<FEC>& getFiniteElementSpace() const
         {
            return m_fes;
         }

         int getDimension() const override
         {
            return m_fes.getRangeDimension();
         }

      private:
         FiniteElementSpace<FEC>& m_fes;
         std::unique_ptr<GridFunction<FEC>> m_gf;
   };
}
#endif
