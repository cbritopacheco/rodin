#include "UMFPack.h"

namespace Rodin::Solver
{
   void UMFPack::solve(Variational::ProblemBase& problem)
   {
      problem.update().assemble();

      m_umfpack.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
      m_umfpack.SetOperator(problem.getStiffnessMatrix());
      m_umfpack.Mult(problem.getMassVector(), problem.getInitialGuess());

      problem.recoverSolution();
   }

   UMFPack& UMFPack::useLongInts(bool v)
   {
      m_useLongInts = v;
      return *this;
   }
}
