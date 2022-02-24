#ifndef RODIN_VARIATIONAL_TRIALFUNCTION_H
#define RODIN_VARIATIONAL_TRIALFUNCTION_H

#include "ShapeFunction.h"
#include "GridFunction.h"

namespace Rodin::Variational
{
   template <>
   class TrialFunction<H1> : public ShapeFunction<H1, Trial>
   {
      public:
         TrialFunction(H1& fes)
            : ShapeFunction<H1, Trial>(fes)
         {}

         TrialFunction(const TrialFunction& other)
            : ShapeFunction<H1, Trial>(other)
         {}

         TrialFunction(TrialFunction&& other)
            : ShapeFunction<H1, Trial>(std::move(other))
         {}

         TrialFunction& emplaceGridFunction()
         {
            m_gf.emplace(this->getFiniteElementSpace());
            return *this;
         }

         GridFunction<H1>& getGridFunction()
         {
            assert(m_gf);
            return *m_gf;
         }

         const GridFunction<H1>& getGridFunction() const
         {
            assert(m_gf);
            return *m_gf;
         }

         TrialFunction* copy() const noexcept override
         {
            return new TrialFunction(*this);
         }
      private:
         std::optional<GridFunction<H1>> m_gf;
   };
   TrialFunction(H1& fes) -> TrialFunction<H1>;
}
#endif

