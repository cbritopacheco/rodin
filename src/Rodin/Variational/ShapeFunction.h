#ifndef RODIN_VARIATIONAL_SHAPEFUNCTION_H
#define RODIN_VARIATIONAL_SHAPEFUNCTION_H

#include <mfem.hpp>

#include "FormLanguage/Base.h"

#include "ForwardDecls.h"

namespace Rodin::Variational
{
   class ShapeFunction : public FormLanguage::Base
   {
      public:
         using ScalarShape = mfem::Vector;
         using VectorShape = mfem::DenseMatrix;
         using MatrixShape = std::vector<mfem::DenseMatrix>;

         enum ValueType
         {
            Scalar,
            Vector,
            Matrix
         };

         enum SpaceType
         {
            Trial,
            Test
         };

         virtual SpaceType getSpaceType() const = 0;

         virtual ValueType getValueType() const = 0;

         virtual void getValue(
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans,
               ScalarShape& values
               ) const;

         virtual void getValue(
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans,
               VectorShape& values
               ) const;

         virtual void getValue(
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans,
               MatrixShape& values
               ) const;

         virtual const FiniteElementSpaceBase& getFiniteElementSpace() const = 0;
   };


}

#endif
