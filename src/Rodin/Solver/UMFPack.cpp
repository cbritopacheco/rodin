#include "UMFPack.h"

namespace Rodin::Solver
{
   void UMFPack::solve(Variational::ProblemBase& problem)
   {
      // Update and assemble problem
      problem.update();
      problem.assemble();

      // Form the linear system
      mfem::SparseMatrix A;
      mfem::Vector B, X;
      problem.getLinearSystem(A, B, X);

      // Solve
      m_umfpack.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
      m_umfpack.SetOperator(A);
      m_umfpack.Mult(B, X);

      // Recover solution
      problem.getSolution(X);
   }

   UMFPack& UMFPack::useLongInts(bool v)
   {
      m_useLongInts = v;
      return *this;
   }
}
