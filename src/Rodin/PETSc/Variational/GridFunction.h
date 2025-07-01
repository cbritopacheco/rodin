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
#include <utility>

#include "Rodin/Context/Local.h"
#include "Rodin/Types.h"
#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Variational/GridFunction.h"

#ifdef RODIN_USE_MPI
#include "Rodin/MPI/Context.h"
#endif

namespace Rodin::Variational
{
  template <class FES>
  class GridFunction<FES, ::Vec> final
    : public GridFunctionBase<GridFunction<FES, ::Vec>>
  {
    public:
      using FESType = FES;

      using DataType = ::Vec;

      using ScalarType = PetscScalar;

      using RangeType = typename FormLanguage::Traits<FESType>::RangeType;

      using ContextType = typename FormLanguage::Traits<FESType>::ContextType;

      using Parent = GridFunctionBase<GridFunction<FESType, DataType>>;
      using Parent::operator=;
      using Parent::min;
      using Parent::max;

      GridFunction(const FESType& fes)
        : Parent(fes)
      {
        init(this->getData());
      }

      GridFunction(const FESType& fes, DataType& data)
        : Parent(fes, data)
      {
        init(data);
      }

      GridFunction(const GridFunction& other)
        : Parent(other.getFiniteElementSpace()),
          m_begin(other.m_begin),
          m_end(other.m_end),
          m_ghosts(other.m_ghosts)
      {
        PetscErrorCode ierr;

        auto& data = this->getData();

        ierr = VecDuplicate(other.getData(), &data);
        assert(ierr == PETSC_SUCCESS);

        ierr = VecCopy(other.getData(), data);
        assert(ierr == PETSC_SUCCESS);

        ierr = VecGetArray(data, &m_raw);
        assert(ierr == PETSC_SUCCESS);
      }

      GridFunction(GridFunction&& other) noexcept
        : Parent(other.getFiniteElementSpace(), std::exchange(other.getData(), nullptr)),
          m_begin(std::move(other.m_begin)),
          m_end(std::move(other.m_end)),
          m_ghosts(std::move(other.m_ghosts)),
          m_raw(std::exchange(other.m_raw, nullptr))
      {}

      GridFunction& operator=(const GridFunction& other) = delete;

      GridFunction& operator=(GridFunction&& other) noexcept
      {
        if (this != &other)
        {
          PetscErrorCode ierr;
          if (m_raw)
          {
            ierr = VecRestoreArray(this->getData(), &m_raw);
            assert(ierr == PETSC_SUCCESS);
          }

          auto& data = this->getData();
          if (data)
          {
            ierr = VecDestroy(&data);
            assert(ierr == PETSC_SUCCESS);
          }

          data = std::exchange(other.getData(), nullptr);
          m_begin = std::move(other.m_begin);
          m_end = std::move(other.m_end);
          m_ghosts = std::move(other.m_ghosts);
          m_raw = std::exchange(other.m_raw, nullptr);
        }
        return *this;
      }

      virtual ~GridFunction()
      {
        PetscErrorCode ierr;
        if (m_raw)
        {
          ierr = VecRestoreArray(this->getData(), &m_raw);
          assert(ierr == PETSC_SUCCESS);
        }

        auto& data = this->getData();
        if (data)
        {
          ierr = VecDestroy(&data);
          assert(ierr == PETSC_SUCCESS);
        }
      }

      void init(DataType& data)
      {
        PetscErrorCode ierr;
        const auto& fes = this->getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        if constexpr (std::is_same_v<ContextType, Context::Local>)
        {
          VecCreate(PETSC_COMM_SELF, &data);
          VecSetType(data, VECSEQ);
          VecSetSizes(data, fes.getSize(), PETSC_DECIDE);
        }
        else if constexpr (std::is_same_v<ContextType, Context::MPI>)
        {
          const auto& ctx = mesh.getContext();
          const auto& comm = ctx.getCommunicator();
          const size_t globalSize = fes.getSize();
          const size_t localSize = fes.getShard().getSize();
          fes.getOwnershipRange(m_begin, m_end);
          const size_t ownedSize = m_end - m_begin;
          size_t ghostSize = localSize - ownedSize;
          m_ghosts.resize(ghostSize);
          for (size_t i = 0; i < ghostSize; ++i)
            m_ghosts[i] = fes.getGlobalIndex(m_end + i);
          ierr = VecCreateGhost(comm, localSize, globalSize, ghostSize, m_ghosts.data(), &data);
          assert(ierr == PETSC_SUCCESS);
        }
        else
        {
          static_assert(std::is_same_v<ContextType, Context::Local> || std::is_same_v<ContextType, Context::MPI>);
        }

        ierr = VecZeroEntries(data);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecGetArray(data, &m_raw);
        assert(ierr == PETSC_SUCCESS);
      }

      constexpr
      ScalarType min(Index& idx) const
      {
        auto& data = this->getData();
        ScalarType res;
        VecMin(data, idx, &res);
        return res;
      }

      constexpr
      ScalarType max(Index& idx) const
      {
        auto& data = this->getData();
        ScalarType res;
        VecMax(data, idx, &res);
        return res;
      }

      ScalarType& operator[](Index global)
      {
        if constexpr (std::is_same_v<ContextType, Context::Local>)
        {
          return reinterpret_cast<ScalarType*>(m_raw)[static_cast<PetscInt>(global)];
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
          return reinterpret_cast<ScalarType*>(m_raw)[local];
        }
        else
        {
          static_assert(
              std::is_same_v<ContextType, Context::Local>
              || std::is_same_v<ContextType, Context::MPI>);
        }
      }

      const ScalarType& operator[](Index global) const
      {
        if constexpr (std::is_same_v<ContextType, Context::Local>)
        {
          return reinterpret_cast<const ScalarType*>(m_raw)[static_cast<PetscInt>(global)];
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
          return reinterpret_cast<const ScalarType*>(m_raw)[local];
        }
        else
        {
          static_assert(
              std::is_same_v<ContextType, Context::Local>
              || std::is_same_v<ContextType, Context::MPI>);
        }
      }

      GridFunction& operator+=(const ScalarType& rhs)
      {
        static_assert(std::is_same_v<RangeType, ScalarType>);
        PetscErrorCode ierr;
        auto& data = this->getData();
        ierr = VecShift(data, rhs);
        assert(ierr == PETSC_SUCCESS);
        return *this;
      }

      GridFunction& operator-=(const ScalarType& rhs)
      {
        static_assert(std::is_same_v<RangeType, ScalarType>);
        PetscErrorCode ierr;
        auto& data = this->getData();
        ierr = VecShift(data, -rhs);
        assert(ierr == PETSC_SUCCESS);
        return *this;
      }

      GridFunction& operator*=(const ScalarType& rhs)
      {
        PetscErrorCode ierr;
        auto& data = this->getData();
        ierr = VecScale(data, rhs);
        assert(ierr == PETSC_SUCCESS);
        return *this;
      }

      GridFunction& operator/=(const ScalarType& rhs)
      {
        PetscErrorCode ierr;
        auto& data = this->getData();
        ierr = VecScale(data, 1.0 / rhs);
        assert(ierr == PETSC_SUCCESS);
        return static_cast<GridFunction&>(*this);
      }

      GridFunction& operator+=(const GridFunction& rhs)
      {
        assert(&this->getFiniteElementSpace() == &rhs.getFiniteElementSpace());
        PetscErrorCode ierr;
        auto& data = this->getData();
        ierr = VecAXPY(data, 1.0, rhs.getData());
        assert(ierr == PETSC_SUCCESS);
        return *this;
      }

      GridFunction& operator-=(const GridFunction& rhs)
      {
        assert(&this->getFiniteElementSpace() == &rhs.getFiniteElementSpace());
        PetscErrorCode ierr;
        auto& data = this->getData();
        ierr = VecAXPY(data, -1.0, rhs.getData());
        assert(ierr == PETSC_SUCCESS);
        return *this;
      }

      GridFunction& operator*=(const GridFunction& rhs)
      {
        PetscErrorCode ierr;
        auto& data = this->getData();
        ierr = VecPointwiseMult(data, data, rhs.getData());
        assert(ierr == PETSC_SUCCESS);
        return *this;
      }

      GridFunction& operator/=(const GridFunction& rhs)
      {
        PetscErrorCode ierr;
        auto& data = this->getData();
        ierr = VecPointwiseDivide(data, data, rhs.getData());
        assert(ierr == PETSC_SUCCESS);
        return *this;
      }

      GridFunction& setData(const DataType& _data, size_t offset = 0)
      {
        auto& data = this->getData();

        PetscErrorCode ierr;
        ierr = VecRestoreArray(data, &m_raw);
        assert(ierr == PETSC_SUCCESS);

        PetscScalar* _raw = nullptr;
        ierr = VecGetArray(_data, &_raw);
        assert(ierr == PETSC_SUCCESS);

        PetscInt localSize;
        ierr = VecGetLocalSize(data, &localSize);
        assert(ierr == PETSC_SUCCESS);

        std::memcpy(m_raw, _raw + static_cast<PetscInt>(offset), localSize * sizeof(PetscScalar));
        ierr = VecRestoreArray(_data, &_raw);
        assert(ierr == PETSC_SUCCESS);

        ierr = VecGhostUpdateBegin(data, INSERT_VALUES, SCATTER_FORWARD);
        assert(ierr == PETSC_SUCCESS);

        ierr = VecGhostUpdateEnd(data, INSERT_VALUES, SCATTER_FORWARD);
        assert(ierr == PETSC_SUCCESS);

        ierr = VecGetArray(data, &m_raw);
        assert(ierr == PETSC_SUCCESS);

        return *this;
      }

    private:
      size_t m_begin, m_end;
      std::vector<PetscInt> m_ghosts;
      PetscScalar* m_raw;
  };
}

#endif
