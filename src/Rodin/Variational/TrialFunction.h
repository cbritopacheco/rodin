#ifndef RODIN_VARIATIONAL_TRIALFUNCTION_H
#define RODIN_VARIATIONAL_TRIALFUNCTION_H

#include "Component.h"
#include "GridFunction.h"
#include "ShapeFunction.h"

namespace Rodin::Variational
{
   template <class FEC>
   class TrialFunction : public ShapeFunction<FEC, Trial>
   {
      public:
         TrialFunction(FiniteElementSpace<FEC>& fes)
            : ShapeFunction<FEC, Trial>(fes)
         {}

         TrialFunction(const TrialFunction& other)
            : ShapeFunction<FEC, Trial>(other)
         {}

         TrialFunction(TrialFunction&& other)
            : ShapeFunction<FEC, Trial>(std::move(other))
         {}

         void operator=(const TrialFunction&) = delete;

         void operator=(TrialFunction&&) = delete;

         TrialFunction& emplaceGridFunction()
         {
            m_gf.emplace(this->getFiniteElementSpace());
            return *this;
         }

         GridFunction<FEC>& getGridFunction()
         {
            assert(m_gf);
            return *m_gf;
         }

         const GridFunction<FEC>& getGridFunction() const
         {
            assert(m_gf);
            return *m_gf;
         }

         Component<TrialFunction<FEC>> x() const
         {
            assert(this->getFiniteElementSpace().getVectorDimension() >= 1);
            return Component<TrialFunction<FEC>>(*this, 0);
         }

         Component<TrialFunction<FEC>> y() const
         {
            assert(this->getFiniteElementSpace().getVectorDimension() >= 2);
            return Component<TrialFunction<FEC>>(*this, 1);
         }

         Component<TrialFunction<FEC>> z() const
         {
            assert(this->getFiniteElementSpace().getVectorDimension() >= 3);
            return Component<TrialFunction<FEC>>(*this, 2);
         }

         ShapeFunctionBase<Trial>& getRoot()  override
         {
            return *this;
         }

         const ShapeFunctionBase<Trial>& getRoot() const override
         {
            return *this;
         }

         TrialFunction* copy() const noexcept override
         {
            return new TrialFunction(*this);
         }
      private:
         std::optional<GridFunction<FEC>> m_gf;
   };
   template <class FEC>
   TrialFunction(FiniteElementSpace<FEC>&) -> TrialFunction<FEC>;
}
#endif

