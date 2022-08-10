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
   template <class FES>
   class TrialFunction : public ShapeFunction<FES, TrialSpace>
   {
      public:
         TrialFunction(FES& fes)
            : ShapeFunction<FES, TrialSpace>(fes)
         {}

         TrialFunction(const TrialFunction& other)
            : ShapeFunction<FES, TrialSpace>(other)
         {}

         TrialFunction(TrialFunction&& other)
            : ShapeFunction<FES, TrialSpace>(std::move(other))
         {}

         void operator=(const TrialFunction&) = delete;

         void operator=(TrialFunction&&) = delete;

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

         const TrialFunction& getLeaf() const override
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
   TrialFunction(FES&) -> TrialFunction<FES>;
}
#endif

