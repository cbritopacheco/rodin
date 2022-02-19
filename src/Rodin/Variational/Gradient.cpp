#include "Gradient.h"

#include "GridFunction.h"
#include "H1.h"

namespace Rodin::Variational
{
   // ---- GridFunction<H1> --------------------------------------------------
   Gradient<GridFunction<H1>>::Gradient(const GridFunction<H1>& u)
      : m_u(u),
        m_mfemVectorCoefficient(&m_u.getHandle())
   {}

   size_t Gradient<GridFunction<H1>>::getDimension() const
   {
      return m_u.getFiniteElementSpace().getMesh().getDimension();
   }

   // ---- TrialFunction<H1> -------------------------------------------------
   Gradient<TrialFunction<H1>>::Gradient(const TrialFunction<H1>& u)
      : m_u(u)
   {}

   ShapeFunction::ValueType Gradient<TrialFunction<H1>>::getValueType() const
   {
      return Vector;
   }

   void Gradient<TrialFunction<H1>>::getValue(
         const mfem::FiniteElement& fe,
         mfem::ElementTransformation& trans,
         VectorShape& values) const
   {
      fe.CalcPhysDShape(trans, values);
   }

   const FiniteElementSpaceBase& Gradient<TrialFunction<H1>>::getFiniteElementSpace() const
   {
      return m_u.getFiniteElementSpace();
   }

   // ---- TestFunction<H1> -------------------------------------------------
   Gradient<TestFunction<H1>>::Gradient(const TestFunction<H1>& v)
      : m_v(v)
   {}

   ShapeFunction::ValueType Gradient<TestFunction<H1>>::getValueType() const
   {
      return Vector;
   }

   void Gradient<TestFunction<H1>>::getValue(
         const mfem::FiniteElement& fe,
         mfem::ElementTransformation& trans,
         VectorShape& values) const
   {
      fe.CalcPhysDShape(trans, values);
   }

   const FiniteElementSpaceBase& Gradient<TestFunction<H1>>::getFiniteElementSpace() const
   {
      return m_v.getFiniteElementSpace();
   }
}

