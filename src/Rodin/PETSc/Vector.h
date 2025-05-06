#ifndef RODIN_PETSC_VECTOR_H
#define RODIN_PETSC_VECTOR_H

#include <petsc.h>
#include <boost/mpi/communicator.hpp>

namespace Rodin::PETSc
{
  class Vector
  {
    public:
      /*
       * @brief Constructor for creating a PETSc vector.
       * @param comm The MPI communicator for the vector.
       *
       * Delegates to VecCreate(comm, &m_vec) to create a new vector object.
       */
      Vector() = default;

      /*
       * @brief Move constructor for transferring ownership of a vector.
       * @param other The vector to move from.
       *
       * Transfers ownership of the PETSc Vec handle from other to this object.
       */
      Vector(Vector&& other) noexcept;

      Vector(const Vec& other) = delete;

      Vector& operator=(const Vector& other) = delete;

      /*
       * @brief Move assignment operator for transferring ownership of a vector.
       * @param other The vector to move from.
       *
       * Transfers ownership of the PETSc Vec handle from other to this object.
       */
      Vector& operator=(Vector&& other) noexcept;

      /*
       * @brief Destructor, destroys the PETSc Vec.
       */
      ~Vector();

      /*
       * @brief Access the underlying PETSc Vec.
       * @return A reference to the PETSc Vec handle.
       */
      ::Vec& getVec();

      /*
       * @brief Access the underlying PETSc Vec (const).
       * @return A const reference to the PETSc Vec handle.
       */
      const ::Vec& getVec() const;

      PetscErrorCode create(MPI_Comm comm);

      /*
       * @brief Sets the local and global sizes of the vector.
       * @param localSize Number of local entries.
       * @param globalSize Total number of entries across all processes.
       */
      PetscErrorCode setSizes(PetscInt localSize, PetscInt globalSize);

      /*
       * @brief Sets Vec options from the command-line/database.
       */
      PetscErrorCode setFromOptions();

      PetscErrorCode zeroEntries();

      PetscErrorCode setValue(PetscInt idx, PetscScalar value, InsertMode mode);

      PetscErrorCode assemblyBegin();

      PetscErrorCode assemblyEnd();

    private:
      ::Vec          m_vec;
  };
}

#endif // RODIN_PETSC_VEC_H

