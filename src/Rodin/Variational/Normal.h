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
         Normal(size_t dimension)
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
            const auto& simplex = p.getSimplex();
            auto& trans = simplex.getTransformation();
            const auto& mesh = simplex.getMesh();
            assert(
               // We are on a face of a d-mesh in d-space
               (
                  (mesh.getDimension() == mesh.getSpaceDimension()) &&
                  simplex.getDimension() == mesh.getDimension() - 1
               ) ||
               // Or we are on an element of a d-mesh in (d + 1)-space.
               (
                  (mesh.getDimension() + 1 == mesh.getSpaceDimension()) &&
                  simplex.getDimension() == mesh.getDimension()
               )
            );
            value.SetSize(m_dimension);
            mfem::CalcOrtho(trans.Jacobian(), value);
            const double norm = value.Norml2();
            assert(norm > 0.0);
            if (simplex.getDimension() == mesh.getDimension() - 1)
            {
               assert(dynamic_cast<const Geometry::Face*>(&simplex));
               const auto& face = static_cast<const Geometry::Face&>(simplex);
               value /= norm * (1.0 - 2.0 * face.isInterface());
            }
            else
            {
               value /= norm;
            }
            return value;
         }

         Normal* copy() const noexcept override
         {
            return new Normal(*this);
         }
      private:
         const size_t m_dimension;
   };
}

#endif
