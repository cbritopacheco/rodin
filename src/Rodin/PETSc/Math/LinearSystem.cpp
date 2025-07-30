#include <mpi.h>
#include <petscsys.h>
#include <utility>

#include "LinearSystem.h"

namespace Rodin::Math
{
  LinearSystem<PETSc::Math::Matrix, PETSc::Math::Vector>::LinearSystem()
    : m_comm(MPI_COMM_NULL)
  {}

  LinearSystem<PETSc::Math::Matrix, PETSc::Math::Vector>::LinearSystem(MPI_Comm comm)
    : m_comm(comm)
  {
    MatCreate(comm, &m_operator);
    VecCreate(comm, &m_solution);
    VecCreate(comm, &m_vector);
  }

  LinearSystem<PETSc::Math::Matrix, PETSc::Math::Vector>::LinearSystem(const LinearSystem& other)
    : Parent(other),
      m_comm(other.m_comm),
      m_operator(other.m_operator),
      m_solution(other.m_solution),
      m_vector(other.m_vector)
  {}

  LinearSystem<PETSc::Math::Matrix, PETSc::Math::Vector>::LinearSystem(LinearSystem&& other) noexcept
    : Parent(std::move(other)),
      m_comm(std::exchange(other.m_comm, MPI_COMM_NULL)),
      m_operator(std::exchange(other.m_operator, PETSC_NULLPTR)),
      m_solution(std::exchange(other.m_solution, PETSC_NULLPTR)),
      m_vector(std::exchange(other.m_vector, PETSC_NULLPTR))
  {}

  LinearSystem<PETSc::Math::Matrix, PETSc::Math::Vector>&
  LinearSystem<PETSc::Math::Matrix, PETSc::Math::Vector>::operator=(const LinearSystem& other)
  {
    if (this != &other)
    {
      Parent::operator=(other);
      m_comm = other.m_comm;
      m_operator = other.m_operator;
      m_solution = other.m_solution;
      m_vector = other.m_vector;
    }
    return *this;
  }

  LinearSystem<PETSc::Math::Matrix, PETSc::Math::Vector>&
  LinearSystem<PETSc::Math::Matrix, PETSc::Math::Vector>::operator=(LinearSystem&& other) noexcept
  {
    if (this != &other)
    {
      Parent::operator=(std::move(other));
      m_comm = std::exchange(other.m_comm, MPI_COMM_NULL);
      m_operator = std::exchange(other.m_operator, PETSC_NULLPTR);
      m_solution = std::exchange(other.m_solution, PETSC_NULLPTR);
      m_vector = std::exchange(other.m_vector, PETSC_NULLPTR);
    }
    return *this;
  }

  LinearSystem<PETSc::Math::Matrix, PETSc::Math::Vector>&
  LinearSystem<PETSc::Math::Matrix, PETSc::Math::Vector>::create(MPI_Comm comm)
  {
    m_comm = comm;
    MatCreate(comm, &m_operator);
    VecCreate(comm, &m_solution);
    VecCreate(comm, &m_vector);
    return *this;
  }

  LinearSystem<PETSc::Math::Matrix, PETSc::Math::Vector>&
  LinearSystem<PETSc::Math::Matrix, PETSc::Math::Vector>::destroy()
  {
    m_comm = MPI_COMM_NULL;
    MatDestroy(&m_operator);
    VecDestroy(&m_solution);
    VecDestroy(&m_vector);
    return *this;
  }
}
