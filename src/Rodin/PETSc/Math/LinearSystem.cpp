#include "LinearSystem.h"

namespace Rodin::Math
{
  LinearSystem<::Mat, ::Vec>::LinearSystem(const boost::mpi::communicator& comm)
    : Parent(::Mat(), ::Vec(), ::Vec())
  {
    PetscErrorCode ierr;
    ierr = MatCreate(comm, &this->getOperator());
    assert(ierr == PETSC_SUCCESS);
    ierr = VecCreate(comm, &this->getVector());
    assert(ierr == PETSC_SUCCESS);
    ierr = VecCreate(comm, &this->getVector());
    assert(ierr == PETSC_SUCCESS);
  }

  LinearSystem<::Mat, ::Vec>::LinearSystem(::Mat&& a, ::Vec&& x, ::Vec&& b) noexcept
    : Parent(::Mat(), ::Vec(), ::Vec())
  {
    this->getOperator() = std::exchange(a, nullptr);
    this->getVector() = std::exchange(x, nullptr);
    this->getSolution() = std::exchange(b, nullptr);
  }

  LinearSystem<::Mat, ::Vec>::LinearSystem(::Mat& a, ::Vec& x, ::Vec& b)
    : Parent(a, x, b)
  {}
}
