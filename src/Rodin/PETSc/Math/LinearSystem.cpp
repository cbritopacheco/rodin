#include <mpi.h>
#include <petscsys.h>

#include "LinearSystem.h"

namespace Rodin::Math
{
  LinearSystem<PETSc::Matrix, PETSc::Vector>::LinearSystem()
    : m_comm(MPI_COMM_NULL)
  {}

  LinearSystem<PETSc::Matrix, PETSc::Vector>::LinearSystem(MPI_Comm comm)
    : m_comm(comm),
      m_operator(comm),
      m_vector(comm),
      m_solution(comm)
  {}

  LinearSystem<PETSc::Matrix, PETSc::Vector>::LinearSystem(const LinearSystem& other)
    : Parent(other),
      m_comm(other.m_comm),
      m_operator(other.m_operator),
      m_vector(other.m_vector),
      m_solution(other.m_solution)
  {}

  LinearSystem<PETSc::Matrix, PETSc::Vector>::LinearSystem(LinearSystem&& other) noexcept
    : Parent(std::move(other)),
      m_comm(std::exchange(other.m_comm, MPI_COMM_NULL)),
      m_operator(std::move(other.m_operator)),
      m_vector(std::move(other.m_vector)),
      m_solution(std::move(other.m_solution))
  {}

  LinearSystem<PETSc::Matrix, PETSc::Vector>&
  LinearSystem<PETSc::Matrix, PETSc::Vector>::operator=(const LinearSystem& other)
  {
    if (this != &other)
    {
      Parent::operator=(other);
      m_comm = other.m_comm;
      m_operator = other.m_operator;
      m_vector = other.m_vector;
      m_solution = other.m_solution;
    }
    return *this;
  }

  LinearSystem<PETSc::Matrix, PETSc::Vector>&
  LinearSystem<PETSc::Matrix, PETSc::Vector>::operator=(LinearSystem&& other) noexcept
  {
    if (this != &other)
    {
      Parent::operator=(std::move(other));
      m_comm = std::exchange(other.m_comm, MPI_COMM_NULL);
      m_operator = std::move(other.m_operator);
      m_vector = std::move(other.m_vector);
      m_solution = std::move(other.m_solution);
    }
    return *this;
  }

  LinearSystem<PETSc::Matrix, PETSc::Vector>&
  LinearSystem<PETSc::Matrix, PETSc::Vector>::create(MPI_Comm comm)
  {
    m_comm = comm;
    m_operator.create(comm);
    m_vector.create(comm);
    m_solution.create(comm);
    return *this;
  }

  LinearSystem<PETSc::Matrix, PETSc::Vector>&
  LinearSystem<PETSc::Matrix, PETSc::Vector>::destroy()
  {
    m_comm = MPI_COMM_NULL;
    m_operator.destroy();
    m_vector.destroy();
    m_solution.destroy();
    return *this;
  }
}
