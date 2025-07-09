#ifndef RODIN_PETSC_VARIATIONAL_TESTFUNCTION_H
#define RODIN_PETSC_VARIATIONAL_TESTFUNCTION_H

#include <petsc.h>

#include "Rodin/Variational/TestFunction.h"

namespace Rodin::PETSc
{
  template <class FES>
  class TestFunction : public Variational::TestFunction<FES>
  {
    public:
      using Parent = Variational::TestFunction<FES>;
      using FESType = FES;

      TestFunction(const FESType& fes)
        : Parent(fes)
      {}
  };

  template <class FES>
  TestFunction(const FES& fes) -> TestFunction<FES>;
}

namespace Rodin::FormLanguage
{
  template <class FES>
  struct Traits<PETSc::TestFunction<FES>>
  {
    using FESType = FES;
  };
}

namespace Rodin::Variational
{
  template <class FES>
  struct IsTestFunction<PETSc::TestFunction<FES>>
  {
    static constexpr Boolean Value = true;
  };
}


#endif

