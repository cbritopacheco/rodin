#include <mpi.h>
#include <petscsys.h>

#include "LinearSystem.h"

namespace Rodin::Math
{
  LinearSystem<PETSc::Math::Matrix, PETSc::Math::Vector>::LinearSystem()
    : m_comm(MPI_COMM_NULL)
  {}

  LinearSystem<PETSc::Math::Matrix, PETSc::Math::Vector>::LinearSystem(MPI_Comm comm)
    : m_comm(comm),
      m_operator(comm),
      m_solution(comm),
      m_vector(comm)
  {}

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
      m_operator(std::move(other.m_operator)),
      m_solution(std::move(other.m_solution)),
      m_vector(std::move(other.m_vector))
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
      m_operator = std::move(other.m_operator);
      m_solution = std::move(other.m_solution);
      m_vector = std::move(other.m_vector);
    }
    return *this;
  }

  LinearSystem<PETSc::Math::Matrix, PETSc::Math::Vector>&
  LinearSystem<PETSc::Math::Matrix, PETSc::Math::Vector>::create(MPI_Comm comm)
  {
    m_comm = comm;
    m_operator.create(comm);
    m_solution.create(comm);
    m_vector.create(comm);
    return *this;
  }

  LinearSystem<PETSc::Math::Matrix, PETSc::Math::Vector>&
  LinearSystem<PETSc::Math::Matrix, PETSc::Math::Vector>::destroy()
  {
    m_comm = MPI_COMM_NULL;
    m_operator.destroy();
    m_solution.destroy();
    m_vector.destroy();
    return *this;
  }
}
