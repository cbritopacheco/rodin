/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_VARIATIONAL_GRIDFUNCTION_H
#define RODIN_PETSC_VARIATIONAL_GRIDFUNCTION_H

#include <petsc.h>
#include <petscmacros.h>
#include <petscsys.h>
#include <petscsystypes.h>
#include <petscvec.h>
#include <utility>

#include "Rodin/Types.h"
#include "Rodin/Context/Local.h"
#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Variational/GridFunction.h"

#include "Rodin/PETSc/Math/Vector.h"

#ifdef RODIN_USE_MPI
#include "Rodin/MPI/Context.h"
#endif

namespace Rodin::Variational
{
  template <class FES>
  class GridFunction<FES, PETSc::Math::Vector>
    : public GridFunctionBase<GridFunction<FES, PETSc::Math::Vector>>
  {
    enum class State
    {
      Acquired,
      Unacquired
    };

    struct ArrayWrite
    {
      State state = State::Unacquired;
      PetscScalar* raw = PETSC_NULLPTR;
    };

    struct ArrayRead
    {
      State state = State::Unacquired;
      const PetscScalar* raw = PETSC_NULLPTR;
    };

    public:
      using FESType =
        FES;

      using DataType =
        PETSc::Math::Vector;

      using ScalarType =
        PetscScalar;

      using RangeType =
        typename FormLanguage::Traits<FESType>::RangeType;

      using FESMeshType =
        typename FormLanguage::Traits<FESType>::MeshType;

      using FESMeshContextType =
        typename FormLanguage::Traits<FESMeshType>::ContextType;

      using Parent =
        GridFunctionBase<GridFunction<FESType, DataType>>;

      using Parent::operator=;

      using Parent::min;

      using Parent::max;

      static_assert(
          std::is_same_v<FESMeshContextType, Context::Local> ||
          std::is_same_v<FESMeshContextType, Context::MPI>);

      GridFunction(const FESType& fes)
        : Parent(fes),
          m_read{.state = State::Unacquired, .raw = PETSC_NULLPTR},
          m_write{.state = State::Unacquired, .raw = PETSC_NULLPTR}
      {
        auto& data = this->getData();
        if constexpr (std::is_same_v<FESMeshContextType, Context::Local>)
        {
          VecCreate(PETSC_COMM_SELF, &data);
          VecSetSizes(data, fes.getSize(), PETSC_DECIDE);
          VecSetFromOptions(data);
        }
        else if constexpr (std::is_same_v<FESMeshContextType, Context::MPI>)
        {
          const auto& mesh = fes.getMesh();
          const auto& ctx = mesh.getContext();
          PetscErrorCode ierr;
          data.create(ctx.getCommunicator());
          const size_t globalSize = fes.getSize();
          const size_t localSize = fes.getShard().getSize();
          fes.getOwnershipRange(m_begin, m_end);
          const size_t ownedSize = m_end - m_begin;
          ierr = VecSetSizes(data, localSize, globalSize);
          assert(ierr == PETSC_SUCCESS);
          ierr = VecSetFromOptions(data);
          assert(ierr == PETSC_SUCCESS);
          size_t ghostSize = localSize - ownedSize;
          m_ghosts.resize(ghostSize);
          for (size_t i = 0; i < ghostSize; ++i)
            m_ghosts[i] = fes.getGlobalIndex(ownedSize + i);
          ierr = VecMPISetGhost(data, ghostSize, m_ghosts.data());
          assert(ierr == PETSC_SUCCESS);
        }
        else
        {
          assert(false);
        }
        VecZeroEntries(data);
      }

      GridFunction(const GridFunction& other)
        : Parent(other.getFiniteElementSpace()),
          m_data(other.getData()),
          m_begin(other.m_begin),
          m_end(other.m_end),
          m_ghosts(other.m_ghosts),
          m_read{.state = State::Unacquired, .raw = PETSC_NULLPTR},
          m_write{.state = State::Unacquired, .raw = PETSC_NULLPTR}
      {}

      GridFunction(GridFunction&& other) noexcept
        : Parent(other.getFiniteElementSpace(), std::exchange(other.getData(), nullptr)),
          m_data(std::move(other.getData())),
          m_begin(std::move(other.m_begin)),
          m_end(std::move(other.m_end)),
          m_ghosts(std::move(other.m_ghosts))
      {
        m_read.state = std::exchange(other.m_read, State::Unacquired);
        m_read.raw = std::exchange(other.m_read.raw, PETSC_NULLPTR);
        m_write.state = std::exchange(other.m_write, State::Unacquired);
        m_write.raw = std::exchange(other.m_write.raw, PETSC_NULLPTR);
      }

      GridFunction& operator=(const GridFunction& other) = delete;

      GridFunction& operator=(GridFunction&& other) noexcept
      {
        if (this != &other)
        {
          Parent::operator=(std::move(other));
          m_data = std::move(other.m_data);
          m_begin = std::move(other.m_begin);
          m_end = std::move(other.m_end);
          m_ghosts = std::move(other.m_ghosts);
          m_read.state = std::exchange(other.m_read, State::Unacquired);
          m_read.raw = std::exchange(other.m_read.raw, PETSC_NULLPTR);
          m_write.state = std::exchange(other.m_write, State::Unacquired);
          m_write.raw = std::exchange(other.m_write.raw, PETSC_NULLPTR);
        }
        return *this;
      }

      virtual ~GridFunction()
      {
        PetscErrorCode ierr;
        if (m_read.state == State::Acquired)
        {
          assert(m_read.raw);
          ierr = VecRestoreArrayRead(m_data, &m_read.raw);
          assert(ierr == PETSC_SUCCESS);
        }
        if (m_write.state == State::Acquired)
        {
          assert(m_write.raw);
          ierr = VecRestoreArrayWrite(m_data, &m_write.raw);
          assert(ierr == PETSC_SUCCESS);
        }
        ierr = VecDestroy(&m_data);
        assert(ierr == PETSC_SUCCESS);
      }

      constexpr
      ScalarType min(Index& idx) const
      {
        flush();
        PetscErrorCode ierr;
        auto& data = this->getData();
        ScalarType res;
        ierr = VecMin(data, idx, &res);
        assert(ierr == PETSC_SUCCESS);
        return res;
      }

      constexpr
      ScalarType max(Index& idx) const
      {
        flush();
        PetscErrorCode ierr;
        auto& data = this->getData();
        ScalarType res;
        ierr = VecMax(data, idx, &res);
        assert(ierr == PETSC_SUCCESS);
        return res;
      }

      ScalarType& operator[](Index global)
      {
        acquire();

        if constexpr (std::is_same_v<FESMeshContextType, Context::Local>)
        {
          return m_write.raw[static_cast<PetscInt>(global)];
        }
        else if constexpr (std::is_same_v<FESMeshContextType, Context::MPI>)
        {
          PetscInt local;
          if (m_begin <= global && global < m_end)
          {
            local = static_cast<PetscInt>(global - m_begin);
          }
          else
          {
            const auto& fes = this->getFiniteElementSpace();
            const Optional<Index> localIdx = fes.getLocalIndex(global);
            assert(localIdx);
            local = *localIdx;
          }
          return m_write.raw[local];
        }
        else
        {
          assert(false);
        }
      }

      const ScalarType& operator[](Index global) const
      {
        acquire();

        if constexpr (std::is_same_v<FESMeshContextType, Context::Local>)
        {
          return m_read.raw[static_cast<PetscInt>(global)];
        }
        else if constexpr (std::is_same_v<FESMeshContextType, Context::MPI>)
        {
          PetscInt local;
          if (m_begin <= global && global < m_end)
          {
            local = static_cast<PetscInt>(global - m_begin);
          }
          else
          {
            const auto& fes = this->getFiniteElementSpace();
            const Optional<Index> localIdx = fes.getLocalIndex(global);
            assert(localIdx);
            local = *localIdx;
          }
          return m_read.raw[local];
        }
        else
        {
          assert(false);
        }
      }

      GridFunction& operator+=(const ScalarType& rhs)
      {
        static_assert(std::is_same_v<RangeType, ScalarType>);
        flush();
        PetscErrorCode ierr;
        auto& data = this->getData();
        ierr = VecShift(data, rhs);
        assert(ierr == PETSC_SUCCESS);
        return *this;
      }

      GridFunction& operator-=(const ScalarType& rhs)
      {
        static_assert(std::is_same_v<RangeType, ScalarType>);
        flush();
        PetscErrorCode ierr;
        auto& data = this->getData();
        ierr = VecShift(data, -rhs);
        assert(ierr == PETSC_SUCCESS);
        return *this;
      }

      GridFunction& operator*=(const ScalarType& rhs)
      {
        flush();
        PetscErrorCode ierr;
        auto& data = this->getData();
        ierr = VecScale(data, rhs);
        assert(ierr == PETSC_SUCCESS);
        return *this;
      }

      GridFunction& operator/=(const ScalarType& rhs)
      {
        flush();
        PetscErrorCode ierr;
        auto& data = this->getData();
        ierr = VecScale(data, 1.0 / rhs);
        assert(ierr == PETSC_SUCCESS);
        return static_cast<GridFunction&>(*this);
      }

      GridFunction& operator+=(const GridFunction& rhs)
      {
        assert(&this->getFiniteElementSpace() == &rhs.getFiniteElementSpace());
        flush();
        PetscErrorCode ierr;
        auto& data = this->getData();
        ierr = VecAXPY(data, 1.0, rhs.getData());
        assert(ierr == PETSC_SUCCESS);
        return *this;
      }

      GridFunction& operator-=(const GridFunction& rhs)
      {
        assert(&this->getFiniteElementSpace() == &rhs.getFiniteElementSpace());
        flush();
        PetscErrorCode ierr;
        auto& data = this->getData();
        ierr = VecAXPY(data, -1.0, rhs.getData());
        assert(ierr == PETSC_SUCCESS);
        return *this;
      }

      GridFunction& operator*=(const GridFunction& rhs)
      {
        flush();
        PetscErrorCode ierr;
        auto& data = this->getData();
        ierr = VecPointwiseMult(data, data, rhs.getData());
        assert(ierr == PETSC_SUCCESS);
        return *this;
      }

      GridFunction& operator/=(const GridFunction& rhs)
      {
        flush();
        PetscErrorCode ierr;
        auto& data = this->getData();
        ierr = VecPointwiseDivide(data, data, rhs.getData());
        assert(ierr == PETSC_SUCCESS);
        return *this;
      }

      GridFunction& setData(const DataType& other, size_t offset = 0)
      {
        PetscErrorCode ierr;
        if constexpr (std::is_same_v<FESMeshContextType, Context::Local>)
        {
          assert(offset == 0);
          ierr = VecCopy(other, m_data);
          assert(ierr == PETSC_SUCCESS);
        }
        else
        {
          // flush();
          // auto& data = this->getData();
          // PetscInt localSize;
          // data.getLocalSize(&localSize);
          // const PetscScalar* src = nullptr;
          // other.getArrayRead(&src);
          // PetscScalar* dst = nullptr;
          // data.getArrayWrite(&dst);
          // std::memcpy(dst, src + static_cast<PetscInt>(offset), localSize * sizeof(PetscScalar));
          // data.restoreArrayWrite(&dst);
          // other.restoreArrayRead(&src);
          // data.ghostUpdateBegin(INSERT_VALUES, SCATTER_FORWARD);
          // data.ghostUpdateEnd(INSERT_VALUES, SCATTER_FORWARD);
        }
        return *this;
      }

      PetscScalar* getRaw()
      {
        return m_write.raw;
      }

      const PetscScalar* getRaw() const
      {
        return m_read.raw;
      }

      void acquire()
      {
        if (m_write.state == State::Unacquired)
        {
          VecGetArrayWrite(m_data, &m_write.raw);
          m_write.state = State::Acquired;
        }
      }

      void acquire() const
      {
        if (m_read.state == State::Unacquired)
        {
          VecGetArrayRead(m_data, &m_read.raw);
          m_read.state = State::Acquired;
        }
      }

      void flush()
      {
        if (m_write.state == State::Acquired)
        {
          VecRestoreArrayWrite(m_data, &m_write.raw);
          m_write.state = State::Unacquired;
        }
      }

      void flush() const
      {
        if (m_read.state == State::Acquired)
        {
          VecRestoreArrayRead(m_data, &m_read.raw);
          m_read.state = State::Unacquired;
        }
      }

      auto& getData()
      {
        return m_data;
      }

      const DataType& getData() const
      {
        return m_data;
      }

      const size_t& getBegin() const
      {
        return m_begin;
      }

      const size_t& getEnd() const
      {
        return m_end;
      }

      const auto& getGhosts() const
      {
        return m_ghosts;
      }

      const ArrayRead& getArrayRead() const
      {
        return m_read;
      }

      const ArrayWrite& getArrayWrite() const
      {
        return m_write;
      }

    private:
      DataType m_data;
      size_t m_begin, m_end;
      std::vector<PetscInt> m_ghosts;
      mutable ArrayRead m_read;
      mutable ArrayWrite m_write;
  };
}

namespace Rodin::PETSc::Variational
{
  template <class FES>
  class GridFunction
    : public Rodin::Variational::GridFunction<FES, PETSc::Math::Vector>
  {
    public:
      using Parent = Rodin::Variational::GridFunction<FES, PETSc::Math::Vector>;
      using Parent::Parent;
      using Parent::operator[];
      using Parent::operator=;
      using Parent::operator+=;
      using Parent::operator-=;
      using Parent::operator*=;
      using Parent::operator/=;
  };

  template <class FES>
  GridFunction(const FES&) -> GridFunction<FES>;
}

#endif
