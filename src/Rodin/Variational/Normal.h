#ifndef RODIN_VARIATIONAL_NORMAL_H
#define RODIN_VARIATIONAL_NORMAL_H

#include "Rodin/Geometry/Mesh.h"

#include "ForwardDecls.h"
#include "VectorFunction.h"

namespace Rodin::Variational
{
   /**
    * @brief Outward unit normal.
    */
   class Normal : public VectorFunctionBase
   {
      public:
         /**
          * @brief Constructs the outward unit normal.
          */
         Normal(int dimension)
            : m_dimension(dimension)
         {
            assert(dimension > 0);
         }

         Normal(const Normal& other)
            :  VectorFunctionBase(other),
               m_dimension(other.m_dimension)
         {}

         Normal(Normal&& other)
            :  VectorFunctionBase(std::move(other)),
               m_dimension(other.m_dimension)
         {}

         int getDimension() const override
         {
            return m_dimension;
         }

         FunctionValue getValue(const Geometry::Point& p) const override
         {
            FunctionValue::Vector value;
            const auto& element = p.getSimplex();
            auto& trans = element.getTransformation();
            const auto& mesh = element.getMesh();
            assert(
               // We are on a boundary element of a d-mesh in d-space
               (
                  mesh.getDimension() == mesh.getSpaceDimension() &&
                  trans.ElementType == mfem::ElementTransformation::BDR_ELEMENT
               ) ||
               // Or we are on an element of a d-mesh in (d + 1)-space.
               (
                  mesh.getDimension() == (mesh.getSpaceDimension() - 1) &&
                  trans.ElementType == mfem::ElementTransformation::ELEMENT
               )
            );
            value.SetSize(m_dimension);
            mfem::CalcOrtho(trans.Jacobian(), value);
            const double norm = value.Norml2();
            assert(norm > 0.0);
            switch (trans.ElementType)
            {
               case mfem::ElementTransformation::BDR_ELEMENT:
               {
                  value /= norm * (
                        1.0 - 2.0 * mesh.getHandle().FaceIsInterior(
                           mesh.getHandle().GetBdrFace(trans.ElementNo)));
                  break;
               }
               default:
               {
                  value /= norm;
               }
            }
            return value;
         }

         Normal* copy() const noexcept override
         {
            return new Normal(*this);
         }
      private:
         const int m_dimension;
   };
}

#endif
