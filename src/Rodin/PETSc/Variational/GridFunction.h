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
    struct ArrayWrite
    {
      Boolean acquired = false;
      ::Vec ghost = PETSC_NULLPTR;
      PetscScalar* raw = PETSC_NULLPTR;
    };

    struct ArrayRead
    {
      Boolean acquired = false;
      ::Vec ghost = PETSC_NULLPTR;
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

      using Parent::projectOnCells;

      using Parent::operator=;

      using Parent::min;

      using Parent::max;

      static_assert(
          std::is_same_v<FESMeshContextType, Context::Local> ||
          std::is_same_v<FESMeshContextType, Context::MPI>);

      GridFunction(const FESType& fes)
        : Parent(fes),
          m_read{.acquired = false, .raw = PETSC_NULLPTR},
          m_write{.acquired = false, .raw = PETSC_NULLPTR}
      {
        PetscErrorCode ierr;
        auto& data = this->getData();
        if constexpr (std::is_same_v<FESMeshContextType, Context::Local>)
        {
          ierr = VecCreate(PETSC_COMM_SELF, &data);
          assert(ierr == PETSC_SUCCESS);

          ierr = VecSetSizes(data, fes.getSize(), PETSC_DECIDE);
          assert(ierr == PETSC_SUCCESS);

          ierr = VecSetFromOptions(data);
          assert(ierr == PETSC_SUCCESS);

          ierr = VecZeroEntries(data);
          assert(ierr == PETSC_SUCCESS);
        }
        else if constexpr (std::is_same_v<FESMeshContextType, Context::MPI>)
        {
          const auto& mesh = fes.getMesh();
          const auto& ctx = mesh.getContext();
          const auto& comm = ctx.getCommunicator();
          const size_t globalSize = fes.getSize();
          const size_t localSize = fes.getShard().getSize();
          fes.getOwnershipRange(m_begin, m_end);
          const size_t ownedSize = m_end - m_begin;

          ierr = VecCreate(comm, &data);
          assert(ierr == PETSC_SUCCESS);

          ierr = VecSetSizes(data, localSize, globalSize);
          assert(ierr == PETSC_SUCCESS);

          ierr = VecSetFromOptions(data);
          assert(ierr == PETSC_SUCCESS);

          size_t ghostSize = localSize - ownedSize;
          m_ghosts.resize(ghostSize);
          for (size_t i = 0; i < ghostSize; ++i)
            m_ghosts[i] = fes.getGlobalIndex(ownedSize + i);

for (int r = 0; r < comm.size(); ++r) {
  if (comm.rank() == r) {
    std::cout << "=== rank " << r << " ===\n";
    for (size_t i = 0; i < m_ghosts.size(); ++i) {
      std::cout << " ghost[" << i << "] = " << m_ghosts[i] << "\n";
    }
    std::cout << std::flush;
  }
  comm.barrier();    // wait so next rank prints cleanly
}

          ierr = VecMPISetGhost(data, ghostSize, m_ghosts.data());
          assert(ierr == PETSC_SUCCESS);

          ierr = VecZeroEntries(data);
          assert(ierr == PETSC_SUCCESS);
        }
        else
        {
          assert(false);
        }
      }

      GridFunction(const GridFunction& other)
        : Parent(other.getFiniteElementSpace()),
          m_data(other.getData()),
          m_begin(other.m_begin),
          m_end(other.m_end),
          m_ghosts(other.m_ghosts),
          m_read{.acquired = false, .raw = PETSC_NULLPTR},
          m_write{.acquired = false, .raw = PETSC_NULLPTR}
      {}

      GridFunction(GridFunction&& other) noexcept
        : Parent(other.getFiniteElementSpace(), std::exchange(other.getData(), nullptr)),
          m_data(std::move(other.getData())),
          m_begin(std::move(other.m_begin)),
          m_end(std::move(other.m_end)),
          m_ghosts(std::move(other.m_ghosts))
      {
        m_read.acquired = std::exchange(other.m_read, false);
        m_read.raw = std::exchange(other.m_read.raw, PETSC_NULLPTR);
        m_write.acquired = std::exchange(other.m_write, false);
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
          m_read.acquired = std::exchange(other.m_read, false);
          m_read.raw = std::exchange(other.m_read.raw, PETSC_NULLPTR);
          m_write.acquired = std::exchange(other.m_write, false);
          m_write.raw = std::exchange(other.m_write.raw, PETSC_NULLPTR);
        }
        return *this;
      }

      virtual ~GridFunction()
      {
        PetscErrorCode ierr;
        if (m_read.acquired)
        {
          assert(m_read.raw);
          ierr = VecRestoreArrayRead(m_data, &m_read.raw);
          assert(ierr == PETSC_SUCCESS);
        }
        if (m_write.acquired)
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
          return m_write.raw[global];
        }
        else if constexpr (std::is_same_v<FESMeshContextType, Context::MPI>)
        {
          const auto& fes = this->getFiniteElementSpace();
          const Optional<Index> localIdx = fes.getLocalIndex(global);
          assert(localIdx);
          return m_write.raw[*localIdx];
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
          return m_read.raw[global];
        }
        else if constexpr (std::is_same_v<FESMeshContextType, Context::MPI>)
        {
          PetscInt local;
          if (m_begin <= global && global < m_end)
          {
            local = global - m_begin;
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
        else if constexpr (std::is_same_v<FESMeshContextType, Context::MPI>)
        {
          flush();
          auto& data = this->getData();
          PetscInt localSize;
          ierr = VecGetLocalSize(data, &localSize);
          assert(ierr == PETSC_SUCCESS);
          const PetscScalar* src = nullptr;
          ierr = VecGetArrayRead(other, &src);
          assert(ierr == PETSC_SUCCESS);
          PetscScalar* dst = nullptr;
          ierr = VecGetArrayWrite(data, &dst);
          assert(ierr == PETSC_SUCCESS);
          std::memcpy(dst, src + offset, localSize * sizeof(PetscScalar));
          ierr = VecRestoreArrayWrite(data, &dst);
          assert(ierr == PETSC_SUCCESS);
          ierr = VecRestoreArrayRead(other, &src);
          assert(ierr == PETSC_SUCCESS);
        }
        else
        {
          assert(false);
        }
        return *this;
      }

      constexpr
      void interpolate(RangeType& res, const Geometry::Point& p) const
      {
        const auto& fes = this->getFiniteElementSpace();
        const auto& polytope = p.getPolytope();
        const size_t d = polytope.getDimension();
        const Index  i = polytope.getIndex();
        const auto& fe = fes.getFiniteElement(d, i);
        const size_t count = fe.getCount();
        RangeType v;
        for (Index local = 0; local < count; ++local)
        {
          const auto mapping = fes.getInverseMapping({ d, i }, fe.getBasis(local));
          mapping(v, p);
          const auto k = this->operator[](fes.getGlobalIndex({ d, i }, local)) * v;
          if (local == 0)
            res = k; // Initializes the result (resizes)
          else
            res += k; // Accumulates the result (does not resize)
        }
      }

      template <class NestedDerived>
      GridFunction& projectOnCells(
          const FunctionBase<NestedDerived>& fn, const FlatSet<Geometry::Attribute>& attrs)
      {
        if constexpr (std::is_same_v<FESMeshContextType, Context::Local>)
        {
          const auto& fes = this->getFiniteElementSpace();
          const auto& mesh = fes.getMesh();
          for (auto it = mesh.getCell(); it; ++it)
          {
            const auto& polytope = *it;
            if (attrs.size() == 0 || attrs.count(polytope.getAttribute()))
              project(fn, { polytope.getDimension(), polytope.getIndex() });
          }
        }
        else if constexpr (std::is_same_v<FESMeshContextType, Context::MPI>)
        {
          const auto& fes = this->getFiniteElementSpace();
          const auto& shard = fes.getShard();
          const auto& mesh = shard.getMesh();
          for (Index i = 0; i < mesh.getCellCount(); ++i)
          {
            const auto it = mesh.getCell(i);
            const auto& polytope = *it;
            if (attrs.size() == 0 || attrs.count(polytope.getAttribute()))
              project(fn, { polytope.getDimension(), polytope.getIndex() });
          }
          this->flush();
        }
        else
        {
          assert(false);
        }
        return *this;
      }

      template <class NestedDerived>
      void project(const FunctionBase<NestedDerived>& fn, const std::pair<size_t, Index>& p)
      {
        const auto& fes = this->getFiniteElementSpace();
        const auto& [d, i] = p;
        if constexpr (std::is_same_v<FESMeshContextType, Context::Local>)
        {
          const auto& fe = fes.getFiniteElement(d, i);
          const auto mapping =
            fes.getMapping({ d, i }, fn.template cast<RangeType>());
          for (Index local = 0; local < fe.getCount(); local++)
          {
            const Index global = fes.getGlobalIndex({ d, i }, local);
            this->operator[](global) = fe.getLinearForm(local)(mapping);
          }
        }
        else if constexpr (std::is_same_v<FESMeshContextType, Context::MPI>)
        {
          const auto& shard = fes.getShard();
          const auto& fe = shard.getFiniteElement(d, i);
          const auto mapping =
            shard.getMapping({ d, i }, fn.template cast<RangeType>());
          for (Index local = 0; local < fe.getCount(); local++)
          {
            const Index global =
              fes.getGlobalIndex(shard.getGlobalIndex({ d, i }, local));
            this->operator[](global) = fe.getLinearForm(local)(mapping);
          }
        }
        else
        {
          assert(false);
        }
      }

      GridFunction& acquire()
      {
        if constexpr (std::is_same_v<FESMeshContextType, Context::Local>)
        {
          if (!m_write.acquired)
          {
            PetscErrorCode ierr;

            ierr = VecGetArrayWrite(m_data, &m_write.raw);
            assert(ierr == PETSC_SUCCESS);

            m_write.acquired = true;
          }
        }
        else if constexpr (std::is_same_v<FESMeshContextType, Context::MPI>)
        {
          if (!m_write.acquired)
          {
            PetscErrorCode ierr;

            ierr = VecGhostUpdateBegin(m_data, INSERT_VALUES, SCATTER_FORWARD);
            assert(ierr == PETSC_SUCCESS);

            ierr = VecGhostUpdateEnd(m_data, INSERT_VALUES, SCATTER_FORWARD);
            assert(ierr == PETSC_SUCCESS);

            ierr = VecGhostGetLocalForm(m_data, &m_write.ghost);
            assert(ierr == PETSC_SUCCESS);

            ierr = VecGetArrayWrite(m_write.ghost, &m_write.raw);
            assert(ierr == PETSC_SUCCESS);

            m_write.acquired = true;
          }
        }
        else
        {
          assert(false);
        }
        return *this;
      }

      const GridFunction& acquire() const
      {
        if constexpr (std::is_same_v<FESMeshContextType, Context::Local>)
        {
          if (!m_read.acquired)
          {
            PetscErrorCode ierr;

            ierr = VecGetArrayRead(m_data, &m_read.raw);
            assert(ierr == PETSC_SUCCESS);

            m_read.acquired = true;
          }
        }
        else if constexpr (std::is_same_v<FESMeshContextType, Context::MPI>)
        {
          if (!m_read.acquired)
          {
            PetscErrorCode ierr;

            ierr = VecGhostUpdateBegin(m_data, INSERT_VALUES, SCATTER_FORWARD);
            assert(ierr == PETSC_SUCCESS);

            ierr = VecGhostUpdateEnd(m_data, INSERT_VALUES, SCATTER_FORWARD);
            assert(ierr == PETSC_SUCCESS);

            ierr = VecGhostGetLocalForm(m_data, &m_read.ghost);
            assert(ierr == PETSC_SUCCESS);

            ierr = VecGetArrayRead(m_read.ghost, &m_read.raw);
            assert(ierr == PETSC_SUCCESS);

            m_read.acquired = true;
          }
        }
        else
        {
          assert(false);
        }
        return *this;
      }

      GridFunction& flush()
      {
        if constexpr (std::is_same_v<FESMeshContextType, Context::Local>)
        {
          if (m_write.acquired)
          {
            PetscErrorCode ierr;

            ierr = VecRestoreArrayWrite(m_write.ghost, &m_write.raw);
            assert(ierr == PETSC_SUCCESS);

            m_write.acquired = false;
          }
        }
        else if constexpr (std::is_same_v<FESMeshContextType, Context::MPI>)
        {
          if (m_write.acquired)
          {
            PetscErrorCode ierr;

            ierr = VecRestoreArrayWrite(m_write.ghost, &m_write.raw);
            assert(ierr == PETSC_SUCCESS);

            ierr = VecGhostRestoreLocalForm(m_data, &m_write.ghost);
            assert(ierr == PETSC_SUCCESS);

            ierr = VecGhostUpdateBegin(m_data, ADD_VALUES, SCATTER_REVERSE);
            assert(ierr == PETSC_SUCCESS);

            ierr = VecGhostUpdateEnd(m_data, ADD_VALUES, SCATTER_REVERSE);
            assert(ierr == PETSC_SUCCESS);

            m_write.acquired = false;
          }
        }
        else
        {
          assert(false);
        }
        return *this;
      }

      const GridFunction& flush() const
      {
        if constexpr (std::is_same_v<FESMeshContextType, Context::Local>)
        {
          if (m_read.acquired)
          {
            PetscErrorCode ierr;

            ierr = VecRestoreArrayRead(m_data, &m_read.raw);
            assert(ierr == PETSC_SUCCESS);

            m_read.acquired = false;
          }
        }
        else if constexpr (std::is_same_v<FESMeshContextType, Context::MPI>)
        {
          if (m_read.acquired)
          {
            PetscErrorCode ierr;

            ierr = VecRestoreArrayRead(m_read.ghost, &m_read.raw);
            assert(ierr == PETSC_SUCCESS);

            ierr = VecGhostRestoreLocalForm(m_data, &m_read.ghost);
            assert(ierr == PETSC_SUCCESS);

            ierr = VecGhostUpdateBegin(m_data, ADD_VALUES, SCATTER_REVERSE);
            assert(ierr == PETSC_SUCCESS);

            ierr = VecGhostUpdateEnd(m_data, ADD_VALUES, SCATTER_REVERSE);
            assert(ierr == PETSC_SUCCESS);

            m_read.acquired = false;
          }
        }
        else
        {
          assert(false);
        }
        return *this;
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
