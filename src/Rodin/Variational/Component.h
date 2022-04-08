#ifndef RODIN_VARIATIONAL_COMPONENT_H
#define RODIN_VARIATIONAL_COMPONENT_H

#include "Rodin/Utility.h"

#include "ForwardDecls.h"
#include "GridFunction.h"
#include "TrialFunction.h"

namespace Rodin::Variational
{
   /**
    * @brief Represents the component (or entry) of a vectorial TrialFunction.
    */
   template <class FEC, class Trait>
   class Component<TrialFunction<FEC, Trait>>
   {
      public:
         Component(const TrialFunction<FEC, Trait>& u, int component)
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

         const TrialFunction<FEC, Trait>& getTrialFunction() const
         {
            return m_u;
         }

         int getIndex() const
         {
            return m_idx;
         }

      private:
         const int m_idx;
         const TrialFunction<FEC, Trait>& m_u;
   };
   template <class FEC, class Trait>
   Component(TrialFunction<FEC, Trait>&, int) -> Component<TrialFunction<FEC, Trait>>;

   /**
    * @brief Represents the component (or entry) of a vectorial GridFunction.
    */
   template <class FEC, class Trait>
   class Component<GridFunction<FEC, Trait>>
   {
      public:
         Component(GridFunction<FEC, Trait>& u, int component)
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

         const GridFunction<FEC, Trait>& getGridFunction() const
         {
            return m_u;
         }

         int getIndex() const
         {
            return m_idx;
         }

         Component& projectOnBoundary(const ScalarFunctionBase& s, int attr)
         {
            return projectOnBoundary(s, std::set<int>{attr});
         }

         Component& projectOnBoundary(const ScalarFunctionBase& s, const std::set<int>& attrs = {})
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
         GridFunction<FEC, Trait>& m_u;
   };
   template <class FEC, class Trait>
   Component(GridFunction<FEC, Trait>&, int) -> Component<GridFunction<FEC, Trait>>;

   /**
    * @brief Represents the component (or entry) of a VectorFunctionBase
    * instance.
    */
   template <>
   class Component<VectorFunctionBase> : public ScalarFunctionBase
   {
      public:
         Component(const VectorFunctionBase& v, int component)
            :  m_v(v.copy()),
               m_idx(component)
         {}

         Component(const Component& other)
            :  ScalarFunctionBase(other),
               m_v(other.m_v->copy()),
               m_idx(other.m_idx)
         {}

         Component(Component&& other)
            :  ScalarFunctionBase(std::move(other)),
               m_v(std::move(other.m_v)),
               m_idx(other.m_idx)
         {}

         int getIndex() const
         {
            return m_idx;
         }

         double getValue(
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) const override
         {
            mfem::Vector v;
            m_v->getValue(v, trans, ip);
            assert(m_idx < v.Size());
            return v(m_idx);
         }

         Component* copy() const noexcept override
         {
            return new Component(*this);
         }
      private:
         std::unique_ptr<VectorFunctionBase> m_v;
         const int m_idx;
   };
   Component(const VectorFunctionBase&, int) -> Component<VectorFunctionBase>;
}

#endif

