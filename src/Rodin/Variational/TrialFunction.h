#ifndef RODIN_VARIATIONAL_TRIALFUNCTION_H
#define RODIN_VARIATIONAL_TRIALFUNCTION_H

#include "Component.h"
#include "GridFunction.h"
#include "ShapeFunction.h"

namespace Rodin::Variational
{
   /**
    * @brief Represents a function which belongs to a trial space
    * @tparam FEC Finite Element Collection
    * @tparam Trait Indicates whether the FiniteElementSpace to which the
    * TrialFunction belongs to is parallel or serial
    */
   template <class FEC, class Trait>
   class TrialFunction : public ShapeFunction<FEC, TrialSpace>
   {
      public:
         TrialFunction(FiniteElementSpace<FEC, Trait>& fes)
            : ShapeFunction<FEC, TrialSpace>(fes)
         {}

         TrialFunction(const TrialFunction& other)
            : ShapeFunction<FEC, TrialSpace>(other)
         {}

         TrialFunction(TrialFunction&& other)
            : ShapeFunction<FEC, TrialSpace>(std::move(other))
         {}

         void operator=(const TrialFunction&) = delete;

         void operator=(TrialFunction&&) = delete;

         TrialFunction& emplaceGridFunction()
         {
            m_gf.emplace(this->getFiniteElementSpace());
            return *this;
         }

         GridFunction<FEC, Trait>& getGridFunction()
         {
            assert(m_gf);
            return *m_gf;
         }

         const GridFunction<FEC, Trait>& getGridFunction() const
         {
            assert(m_gf);
            return *m_gf;
         }

         Component<TrialFunction<FEC, Trait>> x() const
         {
            assert(this->getFiniteElementSpace().getVectorDimension() >= 1);
            return Component<TrialFunction<FEC, Trait>>(*this, 0);
         }

         Component<TrialFunction<FEC, Trait>> y() const
         {
            assert(this->getFiniteElementSpace().getVectorDimension() >= 2);
            return Component<TrialFunction<FEC, Trait>>(*this, 1);
         }

         Component<TrialFunction<FEC, Trait>> z() const
         {
            assert(this->getFiniteElementSpace().getVectorDimension() >= 3);
            return Component<TrialFunction<FEC, Trait>>(*this, 2);
         }

         ShapeFunctionBase<TrialSpace>& getLeaf()  override
         {
            return *this;
         }

         const ShapeFunctionBase<TrialSpace>& getLeaf() const override
         {
            return *this;
         }

         TrialFunction* copy() const noexcept override
         {
            return new TrialFunction(*this);
         }
      private:
         std::optional<GridFunction<FEC, Trait>> m_gf;
   };
   template <class FEC, class Trait>
   TrialFunction(FiniteElementSpace<FEC, Trait>&) -> TrialFunction<FEC, Trait>;
}
#endif

