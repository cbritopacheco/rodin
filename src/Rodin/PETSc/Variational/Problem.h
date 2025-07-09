#ifndef RODIN_PETSC_VARIATIONAL_PROBLEM_H
#define RODIN_PETSC_VARIATIONAL_PROBLEM_H

#include <petsc.h>

#include "Rodin/Variational/Problem.h"

#include "Rodin/PETSc/Math/LinearSystem.h"

namespace Rodin::PETSc
{
  template <class ... Ps>
  class Problem;

  template <class U1, class U2, class ... Us>
  class Problem<Tuple<U1, U2, Us ...>, LinearSystem> : public Variational::Problem<Tuple<U1, U2, Us ...>, LinearSystem>
  {
    using Parent = Variational::Problem<Tuple<U1, U2, Us ...>, LinearSystem>;
    using Parent::operator=;

    public:
      Problem(U1& u1, U2& u2, Us& ... us)
        : Parent(u1, u2, us..., LinearSystem(u1.getFiniteElementSpace().getMesh().getContext().getCommunicator()))
      {
        assert(
            u1.getFiniteElementSpace().getMesh().getContext().getCommunicator()
            ==
            u2.getFiniteElementSpace().getMesh().getContext().getCommunicator());

        assert((... && (
            u2.getFiniteElementSpace().getMesh().getContext().getCommunicator()
            ==
            us.getFiniteElementSpace().getMesh().getContext().getCommunicator())));
      }

      template <class AXB>
      Problem(U1& u1, U2& u2, Us& ... us, AXB&& axb)
        : Parent(u1, u2, us..., std::forward<AXB>(axb))
      {}
  };

  template <class U1, class U2, class ... Us>
  Problem(U1&, U2&, Us&...) -> Problem<Tuple<U1, U2, Us...>, LinearSystem>;

  template <class U1, class U2, class ... Us, class AXB>
  Problem(U1& u1, U2& u2, Us& ... us, AXB&& axb) -> Problem<Tuple<U1, U2, Us...>, AXB>;
}

#endif
