#ifndef RODIN_PETSC_VARIATIONAL_TESTFUNCTION_H
#define RODIN_PETSC_VARIATIONAL_TESTFUNCTION_H

#include <petsc.h>

#include "Rodin/Variational/TestFunction.h"

namespace Rodin::PETSc::Variational
{
  template <class FES>
  class TestFunction : public Rodin::Variational::TestFunction<FES>
  {
    public:
      using FESType =
        FES;

      using Parent =
        Rodin::Variational::TestFunction<FESType>;

      using Parent::Parent;
  };

  template <class FES>
  TestFunction(const FES& fes) -> TestFunction<FES>;
}

namespace Rodin::Variational
{
  template <class FES>
  struct IsTestFunction<PETSc::Variational::TestFunction<FES>>
  {
    static constexpr Boolean Value = true;
  };
}

namespace Rodin::FormLanguage
{
  template <class FES>
  struct Traits<PETSc::Variational::TestFunction<FES>>
  {
    using FESType = FES;
  };
}

#endif

