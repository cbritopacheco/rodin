
#include "Jacobian.h"

namespace Rodin::Variational::Internal
{
   int JacobianShapeR3O::GetRows() const
   {
      return m_sdim;
   }

   int JacobianShapeR3O::GetColumns() const
   {
      return m_vdim;
   }

   int JacobianShapeR3O::GetDOFs() const
   {
      return m_dshape.NumRows() * m_vdim;
   }

   double JacobianShapeR3O::operator()(int row, int col, int dof) const
   {
      assert(0 <= row && row < GetRows());
      assert(0 <= col && col < GetColumns());
      assert(0 <= dof && dof < GetDOFs());
      return (std::floor(dof / m_dshape.NumRows()) == col) * m_dshape(dof % m_dshape.NumRows(), row);
   }

   JacobianShapeR3O& JacobianShapeR3O::operator*=(double s)
   {
      m_dshape *= s;
      return *this;
   }

   JacobianShapeR3O& JacobianShapeR3O::operator=(double s)
   {
      m_dshape = s;
      return *this;
   }

   std::unique_ptr<Rank3Operator> JacobianShapeR3O::Trace() const
   {
      assert(GetRows() == GetColumns());
      auto result = new DenseRank3Operator(1, 1, GetDOFs());
      for (int k = 0; k < GetDOFs(); k++)
         (*result)(0, 0, k) = m_dshape.GetData()[k];
      return std::unique_ptr<Rank3Operator>(result);
   }
}
