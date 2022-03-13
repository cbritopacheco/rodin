#ifndef RODIN_VARIATIONAL_TRIALFUNCTION_H
#define RODIN_VARIATIONAL_TRIALFUNCTION_H

#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>

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
            : ShapeFunction<FES, Trial>(fes),
              m_uuid(boost::uuids::random_generator()())
         {}

         TrialFunction(const TrialFunction& other)
            : ShapeFunction<FES, Trial>(other),
              m_uuid(other.m_uuid)
         {}

         TrialFunction(TrialFunction&& other)
            : ShapeFunction<FES, Trial>(std::move(other)),
              m_uuid(std::move(other.m_uuid))
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
            return Component<TrialFunction<FES>>(*this, 0);
         }

         Component<TrialFunction<FES>> y() const
         {
            return Component<TrialFunction<FES>>(*this, 1);
         }

         Component<TrialFunction<FES>> z() const
         {
            return Component<TrialFunction<FES>>(*this, 2);
         }

         boost::uuids::uuid getUUID() const
         {
            return m_uuid;
         }

         TrialFunction* copy() const noexcept override
         {
            return new TrialFunction(*this);
         }
      private:
         std::optional<GridFunction<FES>> m_gf;
         const boost::uuids::uuid m_uuid;
   };
   template <class FES>
   TrialFunction(FES& fes) -> TrialFunction<FES>;
}
#endif

