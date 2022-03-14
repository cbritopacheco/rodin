#ifndef RODIN_VARIATIONAL_COMPONENT_H
#define RODIN_VARIATIONAL_COMPONENT_H

#include "Rodin/Utility.h"

#include "ForwardDecls.h"
#include "GridFunction.h"
#include "TrialFunction.h"

namespace Rodin::Variational
{
   template <class FES>
   class Component<TrialFunction<FES>>
   {
      public:
         Component(const TrialFunction<FES>& u, int component)
            : m_u(u),
              m_idx(component)
         {}

         Component(const Component& other)
            : m_u(other.m_u),
              m_idx(other.m_idx)
         {}

         Component(Component&& other)
            : m_u(other.m_u),
              m_idx(other.m_idx)
         {}

         const TrialFunction<FES>& getTrialFunction() const
         {
            return m_u;
         }

         int getIndex() const
         {
            return m_idx;
         }

      private:
         const int m_idx;
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
              m_idx(component)
         {}

         Component(const Component& other)
            : m_u(other.m_u),
              m_idx(other.m_component)
         {}

         Component(Component&& other)
            : m_u(other.m_u),
              m_idx(other.m_component)
         {}

         const GridFunction<FES>& getGridFunction() const
         {
            return m_u;
         }

         int getIndex() const
         {
            return m_idx;
         }

         Component& projectOnBoundary(const ScalarCoefficientBase& s, int attr)
         {
            return projectOnBoundary(s, std::set<int>{attr});
         }

         Component& projectOnBoundary(const ScalarCoefficientBase& s, const std::set<int>& attrs = {})
         {
            int maxAttr = *m_u.getFiniteElementSpace()
                              .getMesh()
                              .getBoundaryAttributes().rbegin();
            std::vector<mfem::Coefficient*> mfemCoeffs(
                  m_u.getFiniteElementSpace().getVectorDimension(), nullptr);
            mfemCoeffs[getIndex()] = s.build().release();
            if (attrs.size() == 0)
            {
               mfem::Array<int> marker(maxAttr);
               marker = 1;
               m_u.getHandle().ProjectBdrCoefficient(mfemCoeffs.data(), marker);
            }
            else
            {
               assert(mfemCoeffs[getIndex()] != nullptr);
               mfem::Array<int> marker = Utility::set2marker(attrs, maxAttr);
               m_u.getHandle().ProjectBdrCoefficient(mfemCoeffs.data(), marker);
            }
            delete mfemCoeffs[getIndex()];
            return *this;
         }

      private:
         const int m_idx;
         GridFunction<FES>& m_u;
   };
   template <class FES>
   Component(GridFunction<FES>&, int) -> Component<GridFunction<FES>>;
}

#endif

