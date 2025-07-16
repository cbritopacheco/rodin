#include <petscsys.h>

#include "LinearSystem.h"

namespace Rodin::Math
{
  LinearSystem<PETSc::Matrix, PETSc::Vector>::LinearSystem(const boost::mpi::communicator& comm)
    : m_operator(comm),
      m_vector(comm),
      m_solution(comm)
  {}

  LinearSystem<PETSc::Matrix, PETSc::Vector>::LinearSystem(const LinearSystem& other)
    : Parent(other),
      m_operator(other.m_operator),
      m_vector(other.m_vector),
      m_solution(other.m_solution)
  {}

  LinearSystem<PETSc::Matrix, PETSc::Vector>::LinearSystem(LinearSystem&& other) noexcept
    : Parent(std::move(other)),
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
      m_operator = std::move(other.m_operator);
      m_vector = std::move(other.m_vector);
      m_solution = std::move(other.m_solution);
    }
    return *this;
  }
}
