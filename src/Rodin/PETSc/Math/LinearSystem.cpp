#include <mpi.h>
#include <petscsys.h>
#include <utility>

#include "LinearSystem.h"

namespace Rodin::Math
{
  LinearSystem<::Mat, ::Vec>::LinearSystem(MPI_Comm comm)
    : m_comm(comm)
  {
    PetscErrorCode ierr;
    ierr = MatCreate(comm, &m_operator);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecCreate(comm, &m_solution);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecCreate(comm, &m_vector);
    assert(ierr == PETSC_SUCCESS);
    (void) ierr;
  }

  LinearSystem<::Mat, ::Vec>::LinearSystem(const LinearSystem& other)
    : Parent(other),
      m_comm(other.m_comm),
      m_operator(other.m_operator),
      m_solution(other.m_solution),
      m_vector(other.m_vector),
      m_fieldSplits(other.m_fieldSplits)
  {}

  LinearSystem<::Mat, ::Vec>::LinearSystem(LinearSystem&& other) noexcept
    : Parent(std::move(other)),
      m_comm(std::exchange(other.m_comm, MPI_COMM_NULL)),
      m_operator(std::exchange(other.m_operator, PETSC_NULLPTR)),
      m_solution(std::exchange(other.m_solution, PETSC_NULLPTR)),
      m_vector(std::exchange(other.m_vector, PETSC_NULLPTR)),
      m_fieldSplits(std::move(other.m_fieldSplits))
  {}

  LinearSystem<::Mat, ::Vec>::~LinearSystem()
  {
    m_comm = MPI_COMM_NULL;
    PetscErrorCode ierr;

    ierr = MatDestroy(&m_operator);
    assert(ierr == PETSC_SUCCESS);

    ierr = VecDestroy(&m_solution);
    assert(ierr == PETSC_SUCCESS);

    ierr = VecDestroy(&m_vector);
    assert(ierr == PETSC_SUCCESS);

    (void) ierr;
  }

  LinearSystem<::Mat, ::Vec>&
  LinearSystem<::Mat, ::Vec>::operator=(const LinearSystem& other)
  {
    if (this != &other)
    {
      Parent::operator=(other);
      m_comm = other.m_comm;
      m_operator = other.m_operator;
      m_solution = other.m_solution;
      m_vector = other.m_vector;
      m_fieldSplits = other.m_fieldSplits;
    }
    return *this;
  }

  LinearSystem<::Mat, ::Vec>&
  LinearSystem<::Mat, ::Vec>::operator=(LinearSystem&& other) noexcept
  {
    if (this != &other)
    {
      Parent::operator=(std::move(other));
      m_comm = std::exchange(other.m_comm, MPI_COMM_NULL);
      m_operator = std::exchange(other.m_operator, PETSC_NULLPTR);
      m_solution = std::exchange(other.m_solution, PETSC_NULLPTR);
      m_vector = std::exchange(other.m_vector, PETSC_NULLPTR);
      m_fieldSplits = std::move(other.m_fieldSplits);
    }
    return *this;
  }
}
