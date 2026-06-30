/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_CHECK_H
#define RODIN_PETSC_CHECK_H

#include <petscsys.h>

#include "Rodin/Alert/Exception.h"
#include "Rodin/Alert/Raise.h"

/**
 * @file Check.h
 * @brief Live error checking for PETSc calls.
 *
 * Historically the PETSc backends guarded every call with
 * @c assert(ierr == PETSC_SUCCESS). Release builds define @c NDEBUG, which
 * compiles those asserts out, so a failing PETSc call returned a nonzero error
 * code that was silently dropped and execution continued into corrupt state
 * (e.g. a failed @c MatSetSizes followed by an out-of-range @c MatSetValues).
 *
 * @c RODIN_PETSC_CHECK_OK surfaces the error in every build: it inspects the
 * supplied error code and, when nonzero, raises a Rodin exception carrying the
 * PETSc error number, its message, and the source location before aborting.
 */
#define RODIN_PETSC_CHECK_OK(code)                                            \
  do                                                                          \
  {                                                                           \
    const PetscErrorCode rodinPetscErr_ = (code);                             \
    if (rodinPetscErr_ != PETSC_SUCCESS)                                      \
    {                                                                         \
      const char* rodinPetscMsg_ = nullptr;                                   \
      PetscErrorMessage(rodinPetscErr_, &rodinPetscMsg_, nullptr);            \
      Rodin::Alert::Exception()                                               \
        << "PETSc error " << static_cast<int>(rodinPetscErr_) << " ("         \
        << (rodinPetscMsg_ ? rodinPetscMsg_ : "unknown error") << ") at "     \
        << __FILE__ << ":" << __LINE__                                        \
        << Rodin::Alert::Raise;                                               \
    }                                                                         \
  } while (false)

#endif
