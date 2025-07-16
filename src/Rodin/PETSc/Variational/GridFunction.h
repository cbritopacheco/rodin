/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_VARIATIONAL_GRIDFUNCTION_H
#define RODIN_PETSC_VARIATIONAL_GRIDFUNCTION_H

#include <petsc.h>
#include <petscsys.h>
#include <petscsystypes.h>
#include <petscvec.h>

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
  class GridFunction<FES, PETSc::Vector> final
    : public GridFunctionBase<GridFunction<FES, PETSc::Vector>>
  {
    enum class State
    {
      Acquired,
      Unacquired
    };

    public:
      using FESType = FES;

      using DataType = PETSc::Vector;

      using ScalarType = PetscScalar;

      using RangeType = typename FormLanguage::Traits<FESType>::RangeType;

      using ContextType = typename FormLanguage::Traits<FESType>::ContextType;

      using Parent = GridFunctionBase<GridFunction<FESType, DataType>>;
      using Parent::operator=;
      using Parent::min;
      using Parent::max;

      static_assert(
          std::is_same_v<ContextType, Context::Local> ||
          std::is_same_v<ContextType, Context::MPI>);

      GridFunction(const FESType& fes)
        : Parent(fes),
          m_owned(true),
          m_write(State::Unacquired),
          m_rawWrite(nullptr),
          m_read(State::Unacquired),
          m_rawRead(nullptr)
      {
        auto& data = this->getData();
        if constexpr (std::is_same_v<ContextType, Context::Local>)
        {
          data.create(PETSC_COMM_SELF);
          data.setSizes(fes.getSize(), PETSC_DECIDE).setFromOptions();
        }
        else if constexpr (std::is_same_v<ContextType, Context::MPI>)
        {
          const auto& mesh = fes.getMesh();
          const auto& ctx = mesh.getContext();
          data.create(ctx.getCommunicator());
          const size_t globalSize = fes.getSize();
          const size_t localSize = fes.getShard().getSize();
          fes.getOwnershipRange(m_begin, m_end);
          const size_t ownedSize = m_end - m_begin;
          data.setSizes(localSize, globalSize).setFromOptions();
          size_t ghostSize = localSize - ownedSize;
          m_ghosts.resize(ghostSize);
          for (size_t i = 0; i < ghostSize; ++i)
            m_ghosts[i] = fes.getGlobalIndex(m_end + i);
          data.MPI.setGhost(ghostSize, m_ghosts.data());
        }
        else
        {
          assert(false);
        }
        data.zeroEntries();
      }

      GridFunction(const GridFunction& other)
        : Parent(other.getFiniteElementSpace()),
          m_data(other.getData()),
          m_owned(true),
          m_begin(other.m_begin),
          m_end(other.m_end),
          m_ghosts(other.m_ghosts),
          m_read(State::Unacquired),
          m_rawRead(nullptr),
          m_write(State::Unacquired),
          m_rawWrite(nullptr)
      {}

      GridFunction(GridFunction&& other) noexcept
        : Parent(other.getFiniteElementSpace(), std::exchange(other.getData(), nullptr)),
          m_data(std::move(other.getData())),
          m_owned(std::exchange(other.m_owned, false)),
          m_begin(std::move(other.m_begin)),
          m_end(std::move(other.m_end)),
          m_ghosts(std::move(other.m_ghosts)),
          m_read(std::move(other.m_read)),
          m_rawRead(std::exchange(other.m_rawRead, nullptr)),
          m_write(std::move(other.m_write)),
          m_rawWrite(std::exchange(other.m_rawWrite, nullptr))
      {}

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
          m_read = std::exchange(other.m_read, State::Unacquired);
          m_rawRead = std::exchange(other.m_rawRead, nullptr);
          m_write = std::exchange(other.m_write, State::Unacquired);
          m_rawWrite = std::exchange(other.m_rawWrite, nullptr);
        }
        return *this;
      }

      virtual ~GridFunction()
      {
        if (m_rawRead && m_read == State::Acquired)
          m_data.restoreArrayRead(&m_rawRead);
        if (m_rawWrite && m_write == State::Acquired)
          m_data.restoreArrayWrite(&m_rawWrite);
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

        if constexpr (std::is_same_v<ContextType, Context::Local>)
        {
          return m_rawWrite[static_cast<PetscInt>(global)];
        }
        else if constexpr (std::is_same_v<ContextType, Context::MPI>)
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
          return m_rawWrite[local];
        }
        else
        {
          assert(false);
        }
      }

      const ScalarType& operator[](Index global) const
      {
        acquire();

        if constexpr (std::is_same_v<ContextType, Context::Local>)
        {
          return m_rawRead[static_cast<PetscInt>(global)];
        }
        else if constexpr (std::is_same_v<ContextType, Context::MPI>)
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
          return m_rawRead[local];
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
        flush();
        auto& data = this->getData();
        PetscInt localSize;
        data.getLocalSize(&localSize);
        const PetscScalar* src = nullptr;
        other.getArrayRead(&src);
        PetscScalar* dst = nullptr;
        data.getArrayWrite(&dst);
        std::memcpy(
            dst, src + static_cast<PetscInt>(offset), localSize * sizeof(PetscScalar));
        data.restoreArrayWrite(&dst);
        other.restoreArrayRead(&src);
        data.ghostUpdateBegin(INSERT_VALUES, SCATTER_FORWARD);
        data.ghostUpdateEnd(INSERT_VALUES, SCATTER_FORWARD);
        return *this;
      }

      PetscScalar* getRaw()
      {
        return m_rawWrite;
      }

      const PetscScalar* getRaw() const
      {
        return m_rawRead;
      }

      void acquire()
      {
        if (m_write == State::Unacquired)
        {
          m_data.getArrayWrite(&m_rawWrite);
          m_write = State::Acquired;
        }
      }

      void acquire() const
      {
        if (m_read == State::Unacquired)
        {
          m_data.getArrayRead(&m_rawRead);
          m_read = State::Acquired;
        }
      }

      void flush()
      {
        if (m_write == State::Acquired)
        {
          m_data.restoreArrayWrite(&m_rawWrite);
          m_write = State::Unacquired;
        }
      }

      void flush() const
      {
        if (m_read == State::Acquired)
        {
          m_data.restoreArrayRead(&m_rawRead);
          m_read = State::Unacquired;
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

    private:
      DataType m_data;

      bool m_owned;
      size_t m_begin, m_end;
      std::vector<PetscInt> m_ghosts;

      mutable State m_read;
      mutable const PetscScalar* m_rawRead;

      mutable State m_write;
      mutable PetscScalar* m_rawWrite;
  };
}

namespace Rodin::PETSc
{
  template <class FES>
  using GridFunction = Variational::GridFunction<FES, PETSc::Vector>;
}

#endif
