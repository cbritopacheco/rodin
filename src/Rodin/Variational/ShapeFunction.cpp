#include "ShapeFunction.h"

namespace Rodin::Variational
{
   void ShapeFunction::getValue(const mfem::FiniteElement&,
         mfem::ElementTransformation&, ScalarShape&) const
   {
      assert(false);
   }

   void ShapeFunction::getValue(const mfem::FiniteElement&,
         mfem::ElementTransformation&, VectorShape&) const
   {
      assert(false);
   }

   void ShapeFunction::getValue(const mfem::FiniteElement&,
         mfem::ElementTransformation&, MatrixShape&) const
   {
      assert(false);
   }
}
