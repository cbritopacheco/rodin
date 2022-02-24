#include "UMFPack.h"

namespace Rodin::Solver
{
   void UMFPack::solve(Variational::ProblemBase& problem)
   {
      auto& a = problem.getBilinearForm();
      auto& b = problem.getLinearForm();
      auto& u = problem.getSolution();

      problem.update();
      problem.assemble();

      // Compute essential true degrees of freedom
      int maxBdrAttr = u.getFiniteElementSpace().getMesh().getHandle().bdr_attributes.Max();
      mfem::Array<int> essTrueDofList;
      u.getFiniteElementSpace()
       .getFES()
       .GetEssentialTrueDofs(
                   Utility::set2marker(problem.getEssentialBoundary(), maxBdrAttr), essTrueDofList);

      // Form the linear system
      mfem::SparseMatrix A;
      mfem::Vector B, X;
      a.getHandle()
       .FormLinearSystem(essTrueDofList, u.getHandle(), b.getHandle(), A, X, B);

      // Solve
      m_umfpack.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
      m_umfpack.SetOperator(A);
      m_umfpack.Mult(B, X);

      a.getHandle()
       .RecoverFEMSolution(X, b.getHandle(), u.getHandle());
   }

   UMFPack& UMFPack::useLongInts(bool v)
   {
      m_useLongInts = v;
      return *this;
   }
}
