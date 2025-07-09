#include "LinearSystem.h"
#include <petscsys.h>

namespace Rodin::Math
{
  LinearSystem<::Mat, ::Vec>::LinearSystem(const boost::mpi::communicator& comm)
  {
    PetscErrorCode ierr;
    ierr = MatCreate(comm, &this->getOperator());
    assert(ierr == PETSC_SUCCESS);
    ierr = VecCreate(comm, &this->getVector());
    assert(ierr == PETSC_SUCCESS);
    ierr = VecCreate(comm, &this->getVector());
    assert(ierr == PETSC_SUCCESS);
  }

  LinearSystem<::Mat, ::Vec>::LinearSystem(::Mat& a, ::Vec& x, ::Vec& b)
    : Parent(a, x, b)
  {}

  LinearSystem<::Mat, ::Vec>::LinearSystem(::Mat&& a, ::Vec&& x, ::Vec&& b) noexcept
    : Parent(std::exchange(a, nullptr), std::exchange(x, nullptr), std::exchange(b, nullptr))
  {}

  LinearSystem<::Mat, ::Vec>::LinearSystem(const LinearSystem& other)
    : Parent(other)
  {}

  LinearSystem<::Mat, ::Vec>::LinearSystem(LinearSystem&& other) noexcept
    : Parent(std::move(other))
  {}

  LinearSystem<::Mat, ::Vec>::~LinearSystem()
  {
    PetscErrorCode ierr;
    if (this->getOperator())
    {
      ierr = MatDestroy(&this->getOperator());
      assert(ierr == PETSC_SUCCESS);
    }
    if (this->getVector())
    {
      ierr = VecDestroy(&this->getVector());
      assert(ierr == PETSC_SUCCESS);
    }
    if (this->getSolution())
    {
      ierr = VecDestroy(&this->getSolution());
      assert(ierr == PETSC_SUCCESS);
    }
  }
}
