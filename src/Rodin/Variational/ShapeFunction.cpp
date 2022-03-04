#include <cmath>

#include "ShapeFunction.h"

namespace Rodin::Variational
{}

namespace Rodin::Variational::Internal
{
   int ScalarShapeR3O::GetRows() const
   {
      return m_vdim;
   }

   int ScalarShapeR3O::GetColumns() const
   {
      return 1;
   }

   int ScalarShapeR3O::GetDOFs() const
   {
      return m_shape.Size() * m_vdim;
   }

   double ScalarShapeR3O::operator()(int row, int col, int dof) const
   {
      assert(0 <= row && row < GetRows());
      assert(col == 0);
      assert(0 <= dof && dof < GetDOFs());
      return (std::floor(dof / m_shape.Size()) == row) * m_shape(dof % m_shape.Size());
   }

   ScalarShapeR3O&
   ScalarShapeR3O::operator*=(double s)
   {
      m_shape *= s;
      return *this;
   }

   ScalarShapeR3O&
   ScalarShapeR3O::operator=(double s)
   {
      m_shape = s;
      return *this;
   }

   std::unique_ptr<Rank3Operator>
   ScalarShapeR3O::VectorDot(const mfem::Vector& rhs) const
   {
      assert(GetRows() == rhs.Size());
      auto result = new DenseRank3Operator(1, 1, GetDOFs());
      for (int i = 0; i < m_vdim; i++)
      {
         for (int k = 0; k < m_shape.Size(); k++)
         {
            (*result)(0, 0, k + i * m_shape.Size()) = m_shape(k) * rhs(i);
         }
      }
      return std::unique_ptr<Rank3Operator>(result);
   }
}
