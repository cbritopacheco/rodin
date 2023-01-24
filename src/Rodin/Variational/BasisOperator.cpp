#include "BasisOperator.h"

namespace Rodin::Variational
{
  void DenseBasisOperator::addToVector(mfem::Vector& vec) const
  {
    assert(getRows() == 1 && getColumns() == 1);
    assert(getDOFs() == vec.Size());
    double* vdata = vec.GetData();
    for (int i = 0; i < getDOFs(); i++)
      vdata[i] += operator()(0, 0, i);
  }
}
