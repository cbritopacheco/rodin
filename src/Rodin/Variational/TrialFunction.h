#ifndef RODIN_VARIATIONAL_TRIALFUNCTION_H
#define RODIN_VARIATIONAL_TRIALFUNCTION_H

#include "Component.h"
#include "GridFunction.h"
#include "ShapeFunction.h"

namespace Rodin::Variational
{
   template <class FES>
   class TrialFunction : public ShapeFunction<FES, Trial>
   {
      public:
         TrialFunction(FES& fes)
            : ShapeFunction<FES, Trial>(fes)
         {}

         TrialFunction(const TrialFunction& other)
            : ShapeFunction<FES, Trial>(other)
         {}

         TrialFunction(TrialFunction&& other)
            : ShapeFunction<FES, Trial>(std::move(other))
         {}

         TrialFunction& emplaceGridFunction()
         {
            m_gf.emplace(this->getFiniteElementSpace());
            return *this;
         }

         GridFunction<FES>& getGridFunction()
         {
            assert(m_gf);
            return *m_gf;
         }

         const GridFunction<FES>& getGridFunction() const
         {
            assert(m_gf);
            return *m_gf;
         }

         Component<TrialFunction<FES>> x() const
         {
            assert(this->getFiniteElementSpace().getVectorDimension() >= 1);
            return Component<TrialFunction<FES>>(*this, 0);
         }

         Component<TrialFunction<FES>> y() const
         {
            assert(this->getFiniteElementSpace().getVectorDimension() >= 2);
            return Component<TrialFunction<FES>>(*this, 1);
         }

         Component<TrialFunction<FES>> z() const
         {
            assert(this->getFiniteElementSpace().getVectorDimension() >= 3);
            return Component<TrialFunction<FES>>(*this, 2);
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
         std::optional<GridFunction<FES>> m_gf;
   };
   template <class FES>
   TrialFunction(FES& fes) -> TrialFunction<FES>;
}
#endif

