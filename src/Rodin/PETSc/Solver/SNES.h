#ifndef RODIN_SOLVER_PETSC_SNES_H
#define RODIN_SOLVER_PETSC_SNES_H

#include <petscsnes.h>

#include "Rodin/PETSc/Object.h"
#include "KSP.h"

namespace Rodin::Solver
{
  class SNES : public PETSc::Object<::SNES>
  {
    public:
      using HandleType = ::SNES;
      using KSPType = ::KSP;
      using VectorType = ::Vec;
      using Parent = PETSc::Object<HandleType>;

      explicit SNES(MPI_Comm comm = PETSC_COMM_WORLD);

      virtual ~SNES() override;

      SNES& setType(::SNESType type) noexcept;

      SNES& setTolerances(PetscReal abstol,
                          PetscReal rtol,
                          PetscReal stol,
                          PetscInt maxIt,
                          PetscInt maxF) noexcept;

      SNES& setKSP(KSPType ksp) noexcept;

      void solve(VectorType b, VectorType x);

      HandleType& getHandle() noexcept override;

      const HandleType& getHandle() const noexcept override;

      virtual SNES* copy() const noexcept
      {
        return new SNES(*this);
      }

    private:
      HandleType m_snes;
      ::SNESType m_type;
      PetscReal m_abstol, m_rtol, m_stol;
      PetscInt m_maxIt, m_maxF;
      KSPType m_kspHandle;
  };
}

namespace Rodin::PETSc::Solver
{
  using SNES = Rodin::Solver::SNES;
}

#endif // RODIN_SOLVER_PETSC_SNES_H
