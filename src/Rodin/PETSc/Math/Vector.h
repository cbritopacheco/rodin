#ifndef RODIN_PETSC_MATH_VECTOR_H
#define RODIN_PETSC_MATH_VECTOR_H

#include <boost/mpi/communicator.hpp>

#include <mpi.h>

#include <petsc.h>
#include <petscsys.h>
#include <petscsystypes.h>

#include "Rodin/FormLanguage/Traits.h"

#include "Rodin/PETSc/Object.h"

namespace Rodin::PETSc
{
  class Vector : public Object
  {
    public:
      using Parent = Object;

      class MPI
      {
        public:
          constexpr
          MPI(Vector& vec)
            : m_vec(vec)
          {}

          Vector& setGhost(PetscInt nghost, const PetscInt ghosts[])
          {
            PetscErrorCode ierr;
            ierr = VecMPISetGhost(m_vec.m_vec, nghost, ghosts);
            assert(ierr == PETSC_SUCCESS);
            return m_vec;
          }

        private:
          Vector& m_vec; ///< Reference to the parent Vector object.
      } MPI;

      Vector(const boost::mpi::communicator& comm);

      Vector(const Vector& other);

      Vector(Vector&& other) noexcept;

      virtual ~Vector() override;

      Vector& operator=(const Vector& other);

      Vector& operator=(Vector&& other) noexcept;

      Vector& operator+=(const Vector& other);

      Vector& operator-=(const Vector& other);

      Vector& operator*=(const Vector& rhs);

      Vector& operator/=(const Vector& rhs);

      Vector& operator+=(const PetscScalar& other);

      Vector& operator-=(const PetscScalar& other);

      Vector& operator*=(const PetscScalar& rhs);

      Vector& operator/=(const PetscScalar& rhs);

      Vector& setType(const ::VecType& type);

      Vector& setFromOptions();

      Vector& setSizes(PetscInt localSize, PetscInt globalSize);

      Vector& zeroEntries();

      Vector& setValue(PetscInt idx, const PetscScalar& value, ::InsertMode mode);

      Vector& assemblyBegin();

      Vector& assemblyEnd();

      Vector& axpy(const PetscScalar& alpha, const Vector& x);

      Vector& aypx(const PetscScalar& beta, const Vector& x);

      Vector& axpby(const PetscScalar& alpha, const PetscScalar& beta, const Vector& x);

      Vector& waxpy(const PetscScalar& alpha, const Vector& x, const Vector& y);

      Vector& axpbypcz(
          const PetscScalar& alpha, const Vector& x,
          const PetscScalar& beta, const Vector& y,
          const PetscScalar& gamma);

      const Vector& getArrayRead(const PetscScalar *a[]) const;

      const Vector& restoreArrayRead(const PetscScalar *a[]) const;

      Vector& getArrayWrite(PetscScalar *a[]);

      Vector& restoreArrayWrite(PetscScalar *a[]);

      Vector& ghostUpdateBegin(::InsertMode mode, ::ScatterMode scatterMode);

      Vector& ghostUpdateEnd(::InsertMode mode, ::ScatterMode scatterMode);

      const Vector& getLocalSize(PetscInt* size) const;

      const Vector& getComm(MPI_Comm* comm) const noexcept;

      ::PetscObject& getHandle() noexcept override;

      const ::PetscObject& getHandle() const noexcept;

    private:
      ::Vec m_vec;
  };
}

namespace Rodin::FormLanguage
{
  template <>
  struct Traits<PETSc::Vector>
  {
    using ScalarType = PetscScalar;
  };
}

#endif
