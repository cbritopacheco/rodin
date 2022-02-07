#include "GridFunction.h"
#include "GridFunctionIndex.h"

#include "GridFunctionView.h"

namespace Rodin::Variational
{
   GridFunctionView::GridFunctionView(GridFunctionBase& gf)
      : m_gf(gf)
   {}

   GridFunctionView& GridFunctionView::operator=(double v)
   {
      assert(m_idx);
      for (const auto& i : m_idx->getIndices())
         m_gf.getHandle()[i] = v;
      return *this;
   }
}
