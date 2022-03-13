#ifndef RODIN_VARIATIONAL_COMPONENT_H
#define RODIN_VARIATIONAL_COMPONENT_H

#include "ForwardDecls.h"
#include "TrialFunction.h"

namespace Rodin::Variational
{
   template <class FES>
   class Component<TrialFunction<FES>>
   {
      public:
         Component(const TrialFunction<FES>& u, int component)
            : m_u(u),
              m_component(component)
         {}

         Component(const Component& other)
            : m_u(other.m_u),
              m_component(other.m_component)
         {}

         Component(Component&& other)
            : m_u(other.m_u),
              m_component(other.m_component)
         {}

         const TrialFunction<FES>& getTrialFunction() const
         {
            return m_u;
         }

         int getComponent() const
         {
            return m_component;
         }

      private:
         const int m_component;
         const TrialFunction<FES>& m_u;
   };
   template <class FES>
   Component(TrialFunction<FES>&, int) -> Component<TrialFunction<FES>>;


   template <class FES>
   class Component<GridFunction<FES>>
   {
      public:
         Component(GridFunction<FES>& u, int component)
            : m_u(u),
              m_component(component)
         {}

         Component(const Component& other)
            : m_u(other.m_u),
              m_component(other.m_component)
         {}

         Component(Component&& other)
            : m_u(other.m_u),
              m_component(other.m_component)
         {}

         GridFunction<FES>& getGridFunction()
         {
            return m_u;
         }

         const GridFunction<FES>& getGridFunction() const
         {
            return m_u;
         }

         int getComponent() const
         {
            return m_component;
         }

         GridFunction<FES>& operator=(const ScalarCoefficientBase& v)
         {
            assert(false);
         }

      private:
         const int m_component;
         GridFunction<FES>& m_u;
   };
   template <class FES>
   Component(GridFunction<FES>&, int) -> Component<GridFunction<FES>>;
}

#endif

