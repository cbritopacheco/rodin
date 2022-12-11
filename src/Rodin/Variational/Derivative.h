#ifndef RODIN_VARIATIONAL_DERIVATIVE_H
#define RODIN_VARIATIONAL_DERIVATIVE_H

#include <cassert>
#include <cstdlib>

#include "ForwardDecls.h"
#include "H1.h"
#include "FiniteElementSpace.h"
#include "GridFunction.h"
#include "ScalarFunction.h"

namespace Rodin::Variational
{
   /**
    * @brief Function representing the Derivative of a GridFunction in H1 space
    *
    * Given a GridFunction @f$ u : \mathbb{R}^s \rightarrow \mathbb{R}^d @f$,
    * this class represents the derivative in the @f$ i @f$-th direction of the
    * @f$ j @f$-th component:
    * @f[
    *    \dfrac{\partial u_j}{ \partial x_i } ,
    * @f]
    * where @f$ \quad 0 \leq i < s @f$ and @f$ \ 0 \leq j < d @f$.
    */
   class Derivative : public ScalarFunctionBase
   {
      public:
         /**
          * @brief Constructs the derivative in the i-th direction of the j-th
          * component
          * @param[in] direction Spatial direction @f$ x_i @f$
          * @param[in] component Component @f$ u_j @f$ to differentiate
          * @param[in] u GridFunction in H1 space
          */
         template <class Trait>
         Derivative(int direction, int component, const GridFunction<H1<Trait>>& u)
            :  m_direction(direction),
               m_component(component),
               m_u(u)
         {
            if (u.getRangeType() != RangeType::Scalar)
               UnexpectedRangeTypeException(RangeType::Scalar, u.getRangeType());
         }

         Derivative(const Derivative& other)
            : ScalarFunctionBase(other),
              m_direction(other.m_direction),
              m_component(other.m_component),
              m_u(other.m_u)
         {}

         Derivative(Derivative&& other)
            : ScalarFunctionBase(std::move(other)),
              m_direction(other.m_direction),
              m_component(other.m_component),
              m_u(other.m_u)
         {}

         FunctionValue getValue(const Geometry::Point& p) const override
         {
            mfem::Vector grad;
            switch (p.getElement().getRegion())
            {
               case Geometry::Region::Domain:
               {
                  assert(dynamic_cast<const Geometry::Element*>(&p.getElement()));
                  const auto& element = p.getElement();
                  auto& trans = element.getTransformation();
                  m_u.getHandle().GetGradient(trans, grad);
                  break;
               }
               case Geometry::Region::Boundary:
               {
                  assert(dynamic_cast<const Geometry::Boundary*>(&p.getElement()));
                  const auto& boundary = static_cast<const Geometry::Boundary&>(p.getElement());
                  auto& ft = boundary.getFaceTransformations();
                  assert(ft.Elem1);
                  auto& trans1 = ft.GetElement1Transformation();
                  assert(ft.Elem2);
                  auto& trans2 = ft.GetElement2Transformation();
                  const auto& mesh = p.getElement().getMesh();
                  if (mesh.isSubMesh())
                  {
                     const auto& submesh = static_cast<const Geometry::SubMesh<Context::Serial>&>(mesh);
                     const auto& parent = submesh.getParent();

                     if (getTraceDomain().count(trans1.Attribute))
                     {
                        int parentIdx = submesh.getElementMap().left.at(trans1.ElementNo);
                        trans1.ElementNo = parentIdx;
                        m_u.getHandle().GetGradient(trans1, grad);
                     }
                     else if (getTraceDomain().count(trans2.Attribute))
                     {
                        int parentIdx = submesh.getElementMap().left.at(trans2.ElementNo);
                        const auto& parentElement = parent.getElement(parentIdx);
                        assert(parentElement);
                        m_u.getHandle().GetGradient(parentElement->getTransformation(), grad);
                     }
                     else
                        assert(false);
                  }
                  else
                  {
                     if (getTraceDomain().count(trans1.Attribute))
                        m_u.getHandle().GetGradient(trans1, grad);
                     else if (getTraceDomain().count(trans2.Attribute))
                        m_u.getHandle().GetGradient(trans2, grad);
                     else
                        assert(false);
                  }
                  break;
               }
               case Geometry::Region::Interface:
               {
                  assert(false);
                  break;
               }
            }
            return grad(m_direction);
         }

         Derivative* copy() const noexcept override
         {
            return new Derivative(*this);
         }
      private:
         const int m_direction;
         const int m_component;
         const GridFunctionBase& m_u;
   };

   /**
    * @brief %Utility function for computing @f$ \partial_x u @f$
    * @param[in] u GridFunction in H1 space
    *
    * Given a scalar function @f$ u : \mathbb{R}^s \rightarrow \mathbb{R} @f$,
    * this function constructs the derivative in the @f$ x @f$ direction
    * @f$
    *    \dfrac{\partial u}{\partial x}
    * @f$
    */
   template <class Trait>
   Derivative Dx(GridFunction<H1<Trait>>& u)
   {
      assert(u.getRangeType() == RangeType::Scalar);
      return Derivative(0, 0, u);
   }

   /**
    * @brief %Utility function for computing @f$ \partial_y u @f$
    * @param[in] u GridFunction in H1 space
    *
    * Given a scalar function @f$ u : \mathbb{R}^s \rightarrow \mathbb{R} @f$,
    * this function constructs the derivative in the @f$ y @f$ direction
    * @f$
    *    \dfrac{\partial u}{\partial y}
    * @f$
    */
   template <class Trait>
   Derivative Dy(GridFunction<H1<Trait>>& u)
   {
      assert(u.getRangeType() == RangeType::Scalar);
      return Derivative(1, 0, u);
   }

   /**
    * @brief %Utility function for computing @f$ \partial_z u @f$
    * @param[in] u GridFunction in H1 space
    *
    * Given a scalar function @f$ u : \mathbb{R}^s \rightarrow \mathbb{R} @f$,
    * this function constructs the derivative in the @f$ y @f$ direction
    * @f$
    *    \dfrac{\partial u}{\partial z}
    * @f$
    */
   template <class Trait>
   Derivative Dz(GridFunction<H1<Trait>>& u)
   {
      assert(u.getRangeType() == RangeType::Scalar);
      return Derivative(2, 0, u);
   }
}

#endif
