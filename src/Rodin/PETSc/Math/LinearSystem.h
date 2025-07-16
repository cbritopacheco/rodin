/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_MATH_LINEARSYSTEM_H
#define RODIN_PETSC_MATH_LINEARSYSTEM_H

#include <boost/mpi/communicator.hpp>
#include <petsc.h>
#include "Rodin/Math/LinearSystem.h"
#include "Rodin/PETSc/Math/Matrix.h"
#include "Rodin/PETSc/Math/Vector.h"

namespace Rodin::Math
{
  template <>
  class LinearSystem<PETSc::Matrix, PETSc::Vector>
    : public LinearSystemBase<PETSc::Matrix, PETSc::Vector, LinearSystem<PETSc::Matrix, PETSc::Vector>>
  {
    public:
      using MatrixType =
        PETSc::Matrix;

      using VectorType =
        PETSc::Vector;

      using Parent =
        LinearSystemBase<MatrixType, VectorType, LinearSystem<MatrixType, VectorType>>;

      LinearSystem() = default;

      LinearSystem(const boost::mpi::communicator& comm);

      LinearSystem(const LinearSystem& other);

      LinearSystem(LinearSystem&& other) noexcept;

      virtual ~LinearSystem() = default;

      LinearSystem& operator=(const LinearSystem& other);

      LinearSystem& operator=(LinearSystem&& other) noexcept;

      template <class DOFScalar>
      LinearSystem& eliminate(const IndexMap<DOFScalar>& dofs, size_t offset = 0)
      {
        auto& A = this->getOperator();
        auto& b = this->getVector();
        auto& x = this->getSolution();
        std::vector<PetscInt> rows;
        rows.reserve(dofs.size());
        for (auto const& kv : dofs)
          rows.push_back(PetscInt(kv.first) + PetscInt(offset));
        for (auto const& kv : dofs)
        {
          const PetscInt  i   = PetscInt(kv.first) + PetscInt(offset);
          const auto& ui  = kv.second;
          x.setValue(i, ui, INSERT_VALUES);
        }
        x.assemblyBegin();
        x.assemblyEnd();
        A.zeroRowsColumns(rows.size(), rows.data(), 1.0, x, b);
        return *this;
      }

      template <class DOFScalar>
      LinearSystemBase& merge(
          const IndexMap<std::pair<IndexArray, Math::Vector<DOFScalar>>>& dofs, size_t offset = 0)
      {
        throw "Unimplemented.";
        return *this;
      }

      constexpr
      MatrixType& getOperator()
      {
        return m_operator;
      }

      constexpr
      const MatrixType& getOperator() const
      {
        return m_operator;
      }

      constexpr
      VectorType& getVector()
      {
        return m_vector;
      }

      constexpr
      const VectorType& getVector() const
      {
        return m_vector;
      }

      constexpr
      VectorType& getSolution()
      {
        return m_solution;
      }

      constexpr
      const VectorType& getSolution() const
      {
        return m_solution;
      }

    private:
      MatrixType m_operator; ///< The operator of the linear system.
      VectorType m_vector;   ///< The vector of the linear system.
      VectorType m_solution; ///< The solution vector of the linear system.
  };
}

namespace Rodin::PETSc
{
  using LinearSystem = Math::LinearSystem<PETSc::Matrix, PETSc::Vector>;
}

#endif
