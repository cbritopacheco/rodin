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
  class GridFunction<FES, ::Vec>
    : public GridFunctionBase<GridFunction<FES, ::Vec>>
  {
    public:
      struct GhostBimap
      {
        std::vector<PetscInt> left;
        FlatMap<PetscInt, PetscInt> right;
      };

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

      using FESType =
        FES;

      using DataType =
        ::Vec;

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

      using Parent::project;

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

          ierr = VecSetUp(data);
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
          const auto& shard = fes.getShard();
          fes.getOwnershipRange(m_begin, m_end);
          const size_t ownedSize = m_end - m_begin;

          Index offset = 0;
          for (size_t i = 0; i < shard.getSize(); ++i)
          {
            const Index global = fes.getGlobalIndex(i);
            if (global < m_begin || global >= m_end)
            {
              auto [it, inserted] = m_ghosts.right.emplace(global, offset);
              assert(inserted);
              offset++;
            }
          }

          m_ghosts.left.resize(m_ghosts.right.size());
          for (const auto& [global, offset] : m_ghosts.right)
          {
            assert(offset >= 0);
            assert(offset < m_ghosts.left.size());
            m_ghosts.left[offset] = global;
          }

          ierr = VecCreate(comm, &data);
          assert(ierr == PETSC_SUCCESS);

          ierr = VecSetSizes(data, ownedSize, globalSize);
          assert(ierr == PETSC_SUCCESS);

          ierr = VecSetFromOptions(data);
          assert(ierr == PETSC_SUCCESS);

          ierr = VecMPISetGhost(data, m_ghosts.left.size(), m_ghosts.left.data());
          assert(ierr == PETSC_SUCCESS);

          ierr = VecZeroEntries(data);
          assert(ierr == PETSC_SUCCESS);
        }
        else
        {
          assert(false);
        }
        (void) ierr;
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
          m_begin(std::exchange(other.m_begin, 0)),
          m_end(std::exchange(other.m_end, 0)),
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
          m_begin = std::exchange(other.m_begin, 0);
          m_end = std::exchange(other.m_end, 0);
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
          m_read.acquired = false;
        }
        if (m_write.acquired)
        {
          assert(m_write.raw);
          ierr = VecRestoreArrayWrite(m_data, &m_write.raw);
          assert(ierr == PETSC_SUCCESS);
          m_write.acquired = false;
        }
        ierr = VecDestroy(&m_data);
        assert(ierr == PETSC_SUCCESS);
        m_data = PETSC_NULLPTR;
        (void) ierr;
      }

      constexpr
      ScalarType min(Index& idx) const
      {
        this->flush();
        PetscErrorCode ierr;
        auto& data = this->getData();
        ScalarType res;
        ierr = VecMin(data, idx, &res);
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;
        return res;
      }

      constexpr
      ScalarType max(Index& idx) const
      {
        this->flush();
        PetscErrorCode ierr;
        auto& data = this->getData();
        ScalarType res;
        ierr = VecMax(data, idx, &res);
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;
        return res;
      }

      ScalarType& operator[](Index global)
      {
        this->acquire();

        if constexpr (std::is_same_v<FESMeshContextType, Context::Local>)
        {
          return m_write.raw[global];
        }
        else if constexpr (std::is_same_v<FESMeshContextType, Context::MPI>)
        {
          PetscInt local;
          if (m_begin <= global && global < m_end)
            local = global - m_begin;
          else
            local = m_end - m_begin + m_ghosts.right.at(global);
          return m_write.raw[local];
        }
        else
        {
          assert(false);
        }
      }

      const ScalarType& operator[](Index global) const
      {
        this->acquire();

        if constexpr (std::is_same_v<FESMeshContextType, Context::Local>)
        {
          return m_read.raw[global];
        }
        else if constexpr (std::is_same_v<FESMeshContextType, Context::MPI>)
        {
          PetscInt local;
          if (m_begin <= global && global < m_end)
            local = global - m_begin;
          else
            local = m_end - m_begin + m_ghosts.right.at(global);
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
        this->flush();
        PetscErrorCode ierr;
        auto& data = this->getData();
        ierr = VecShift(data, rhs);
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;
        return *this;
      }

      GridFunction& operator-=(const ScalarType& rhs)
      {
        static_assert(std::is_same_v<RangeType, ScalarType>);
        this->flush();
        PetscErrorCode ierr;
        auto& data = this->getData();
        ierr = VecShift(data, -rhs);
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;
        return *this;
      }

      GridFunction& operator*=(const ScalarType& rhs)
      {
        this->flush();
        PetscErrorCode ierr;
        auto& data = this->getData();
        ierr = VecScale(data, rhs);
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;
        return *this;
      }

      GridFunction& operator/=(const ScalarType& rhs)
      {
        this->flush();
        PetscErrorCode ierr;
        auto& data = this->getData();
        ierr = VecScale(data, 1.0 / rhs);
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;
        return static_cast<GridFunction&>(*this);
      }

      GridFunction& operator+=(const GridFunction& rhs)
      {
        assert(&this->getFiniteElementSpace() == &rhs.getFiniteElementSpace());
        this->flush();
        PetscErrorCode ierr;
        auto& data = this->getData();
        ierr = VecAXPY(data, 1.0, rhs.getData());
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;
        return *this;
      }

      GridFunction& operator-=(const GridFunction& rhs)
      {
        assert(&this->getFiniteElementSpace() == &rhs.getFiniteElementSpace());
        this->flush();
        PetscErrorCode ierr;
        auto& data = this->getData();
        ierr = VecAXPY(data, -1.0, rhs.getData());
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;
        return *this;
      }

      GridFunction& operator*=(const GridFunction& rhs)
      {
        this->flush();
        PetscErrorCode ierr;
        auto& data = this->getData();
        ierr = VecPointwiseMult(data, data, rhs.getData());
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;
        return *this;
      }

      GridFunction& operator/=(const GridFunction& rhs)
      {
        this->flush();
        PetscErrorCode ierr;
        auto& data = this->getData();
        ierr = VecPointwiseDivide(data, data, rhs.getData());
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;
        return *this;
      }

      GridFunction& setData(const DataType& data, size_t offset = 0)
      {
        this->flush();

        PetscErrorCode ierr;

        // Global size checks
        PetscInt nDst = 0, nSrc = 0;
        ierr = VecGetSize(m_data, &nDst); assert(ierr == PETSC_SUCCESS);
        ierr = VecGetSize(data,   &nSrc); assert(ierr == PETSC_SUCCESS);

        const PetscInt off = static_cast<PetscInt>(offset);
        assert(off >= 0);
        assert(off + nDst <= nSrc);

        if constexpr (std::is_same_v<FESMeshContextType, Context::Local>)
        {
          // Sequential: simple local slice
          IS is = nullptr;
          ierr = ISCreateStride(PETSC_COMM_SELF, nDst, off, 1, &is);
          assert(ierr == PETSC_SUCCESS);

          Vec sub = nullptr;
          ierr = VecGetSubVector(data, is, &sub);
          assert(ierr == PETSC_SUCCESS);

          ierr = VecCopy(sub, m_data);
          assert(ierr == PETSC_SUCCESS);

          ierr = VecRestoreSubVector(data, is, &sub);
          assert(ierr == PETSC_SUCCESS);

          ierr = ISDestroy(&is);
          assert(ierr == PETSC_SUCCESS);
        }
        else
        {
          static_assert(std::is_same_v<FESMeshContextType, Context::MPI>);

          // MPI: copy only the owned range of m_data on this rank
          PetscInt rb = 0, re = 0;
          ierr = VecGetOwnershipRange(m_data, &rb, &re);
          assert(ierr == PETSC_SUCCESS);

          const PetscInt nLocal = re - rb;
          assert(nLocal >= 0);

          // Indices in 'data' we want are shifted by 'off'
          const PetscInt start = off + rb;

          // Important: IS local size must be nLocal, not nDst (global)
          MPI_Comm comm;
          PetscObjectGetComm((PetscObject) m_data, &comm);

          IS is = nullptr;
          ierr = ISCreateStride(comm, nLocal, start, 1, &is);
          assert(ierr == PETSC_SUCCESS);

          Vec sub = nullptr;
          ierr = VecGetSubVector(data, is, &sub);
          assert(ierr == PETSC_SUCCESS);

          // Now sub and m_data have identical local layouts -> VecCopy is valid
          ierr = VecCopy(sub, m_data);
          assert(ierr == PETSC_SUCCESS);

          ierr = VecRestoreSubVector(data, is, &sub);
          assert(ierr == PETSC_SUCCESS);

          ierr = ISDestroy(&is);
          assert(ierr == PETSC_SUCCESS);

          // Refresh ghosts from owned values
          ierr = VecGhostUpdateBegin(m_data, INSERT_VALUES, SCATTER_FORWARD);
          assert(ierr == PETSC_SUCCESS);
          ierr = VecGhostUpdateEnd(m_data, INSERT_VALUES, SCATTER_FORWARD);
          assert(ierr == PETSC_SUCCESS);
        }

        (void)ierr;
        return *this;
      }

      constexpr
      void interpolate(RangeType& res, const Geometry::Point& p) const
      {
        const auto& fes = this->getFiniteElementSpace();
        const auto& polytope = p.getPolytope();
        const size_t d = polytope.getDimension();
        const Index  i = polytope.getIndex();

        if constexpr (std::is_same_v<FESMeshContextType, Context::Local>)
        {
          const auto& fe = fes.getFiniteElement(d, i);
          const size_t count = fe.getCount();
          RangeType v;
          for (Index local = 0; local < count; ++local)
          {
            const auto mapping = fes.getPushforward({ d, i }, fe.getBasis(local));
            const auto k = this->operator[](fes.getGlobalIndex({ d, i }, local)) * mapping(p);
            if (local == 0)
              res = k; // Initializes the result (resizes)
            else
              res += k; // Accumulates the result (does not resize)
          }
        }
        else if (std::is_same_v<FESMeshContextType, Context::MPI>)
        {
          const auto& shard = fes.getShard();
          const auto& fe = shard.getFiniteElement(d, i);
          const size_t count = fe.getCount();
          RangeType v;
          for (Index local = 0; local < count; ++local)
          {
            const auto mapping = fes.getPushforward({ d, i }, fe.getBasis(local));
            const auto k = this->operator[](fes.getGlobalIndex({ d, i }, local)) * mapping(p);
            if (local == 0)
              res = k; // Initializes the result (resizes)
            else
              res += k; // Accumulates the result (does not resize)
          }
        }
        else
        {
          assert(false);
        }
      }

      template <class NestedDerived, class Predicate>
      GridFunction& projectOnCells(const FunctionBase<NestedDerived>& fn, const Predicate& pred)
      {
        this->acquire();

        if constexpr (std::is_same_v<FESMeshContextType, Context::Local>)
        {
          const auto& fes = this->getFiniteElementSpace();
          const auto& mesh = fes.getMesh();
          for (auto it = mesh.getCell(); it; ++it)
          {
            const auto& polytope = *it;
            if (pred(polytope))
              this->project(fn, { polytope.getDimension(), polytope.getIndex() });
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
            if (pred(polytope))
              this->project(fn, { polytope.getDimension(), polytope.getIndex() });
          }
        }
        else
        {
          assert(false);
        }

        this->flush();

        return *this;
      }

      template <class Function, class Pred>
      GridFunction& project(const Geometry::Region& region, const Function& v, const Pred& pred)
      {
        const auto& fes = this->getFiniteElementSpace();
        const auto& mesh = fes.getMesh();

        Geometry::PolytopeIterator it;
        switch (region)
        {
          case Geometry::Region::Cells:
          {
            if constexpr (std::is_same_v<FESMeshContextType, Context::Local>)
            {
              it = mesh.getCell();
            }
            else
            {
              static_assert(std::is_same_v<FESMeshContextType, Context::MPI>);
              it = mesh.getShard().getCell();
            }
            break;
          }
          case Geometry::Region::Faces:
          {
            if constexpr (std::is_same_v<FESMeshContextType, Context::Local>)
            {
              it = mesh.getFace();
            }
            else
            {
              static_assert(std::is_same_v<FESMeshContextType, Context::MPI>);
              it = mesh.getShard().getFace();
            }
            break;
          }
          case Geometry::Region::Boundary:
          {
            if constexpr (std::is_same_v<FESMeshContextType, Context::Local>)
            {
              it = mesh.getBoundary();
            }
            else
            {
              static_assert(std::is_same_v<FESMeshContextType, Context::MPI>);
              it = mesh.getShard().getBoundary();
            }
            break;
          }
          case Geometry::Region::Interface:
          {
            if constexpr (std::is_same_v<FESMeshContextType, Context::Local>)
            {
              it = mesh.getInterface();
            }
            else
            {
              static_assert(std::is_same_v<FESMeshContextType, Context::MPI>);
              it = mesh.getShard().getInterface();
            }
            break;
          }
        }

        this->acquire();

        while (it)
        {
          const auto& polytope = *it;
          if (pred(polytope))
            this->project({ polytope.getDimension(), polytope.getIndex() }, v);
          ++it;
        }

        this->flush();

        return *this;
      }

      template <class Function>
      void project(const std::pair<size_t, Index>& p, const Function& fn)
      {
        const auto& fes = this->getFiniteElementSpace();
        const auto& [d, i] = p;
        if constexpr (std::is_same_v<FESMeshContextType, Context::Local>)
        {
          const auto& fe = fes.getFiniteElement(d, i);
          const auto mapping = fes.getPullback({ d, i }, fn);
          for (Index local = 0; local < fe.getCount(); local++)
          {
            const Index global = fes.getGlobalIndex({ d, i }, local);
            this->operator[](global) = fe.getLinearForm(local)(mapping);
          }
        }
        else
        {
          static_assert(std::is_same_v<FESMeshContextType, Context::MPI>);
          const auto& shard = fes.getShard();
          const auto& fe = shard.getFiniteElement(d, i);
          const auto mapping = shard.getPullback({ d, i }, fn);
          for (Index local = 0; local < fe.getCount(); local++)
          {
            const Index global =
              fes.getGlobalIndex(
                  shard.getGlobalIndex({ d, i }, local));
            this->operator[](global) = fe.getLinearForm(local)(mapping);
          }
        }
      }

      GridFunction& acquire()
      {
        PetscErrorCode ierr;
        if constexpr (std::is_same_v<FESMeshContextType, Context::Local>)
        {
          if (!m_write.acquired)
          {
            ierr = VecGetArrayWrite(m_data, &m_write.raw);
            assert(ierr == PETSC_SUCCESS);

            m_write.acquired = true;
          }
        }
        else if constexpr (std::is_same_v<FESMeshContextType, Context::MPI>)
        {
          if (!m_write.acquired)
          {
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
        (void) ierr;
        return *this;
      }

      const GridFunction& acquire() const
      {
        PetscErrorCode ierr;
        if constexpr (std::is_same_v<FESMeshContextType, Context::Local>)
        {
          if (!m_read.acquired)
          {
            ierr = VecGetArrayRead(m_data, &m_read.raw);
            assert(ierr == PETSC_SUCCESS);

            m_read.acquired = true;
          }
        }
        else if constexpr (std::is_same_v<FESMeshContextType, Context::MPI>)
        {
          if (!m_read.acquired)
          {
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
        (void) ierr;
        return *this;
      }

      GridFunction& flush()
      {
        PetscErrorCode ierr;
        if constexpr (std::is_same_v<FESMeshContextType, Context::Local>)
        {
          if (m_write.acquired)
          {
            ierr = VecRestoreArrayWrite(m_data, &m_write.raw);
            assert(ierr == PETSC_SUCCESS);

            m_write.acquired = false;
          }
        }
        else if constexpr (std::is_same_v<FESMeshContextType, Context::MPI>)
        {
          if (m_write.acquired)
          {
            ierr = VecRestoreArrayWrite(m_write.ghost, &m_write.raw);
            assert(ierr == PETSC_SUCCESS);

            ierr = VecGhostRestoreLocalForm(m_data, &m_write.ghost);
            assert(ierr == PETSC_SUCCESS);

            ierr = VecGhostUpdateBegin(m_data, INSERT_VALUES, SCATTER_REVERSE);
            assert(ierr == PETSC_SUCCESS);

            ierr = VecGhostUpdateEnd(m_data, INSERT_VALUES, SCATTER_REVERSE);
            assert(ierr == PETSC_SUCCESS);

            m_write.acquired = false;
          }
        }
        else
        {
          assert(false);
        }
        (void) ierr;
        return *this;
      }

      const GridFunction& flush() const
      {
        PetscErrorCode ierr;
        if constexpr (std::is_same_v<FESMeshContextType, Context::Local>)
        {
          if (m_read.acquired)
          {
            ierr = VecRestoreArrayRead(m_data, &m_read.raw);
            assert(ierr == PETSC_SUCCESS);

            m_read.acquired = false;
          }
        }
        else if constexpr (std::is_same_v<FESMeshContextType, Context::MPI>)
        {
          if (m_read.acquired)
          {
            ierr = VecRestoreArrayRead(m_read.ghost, &m_read.raw);
            assert(ierr == PETSC_SUCCESS);

            ierr = VecGhostRestoreLocalForm(m_data, &m_read.ghost);
            assert(ierr == PETSC_SUCCESS);

            ierr = VecGhostUpdateBegin(m_data, INSERT_VALUES, SCATTER_REVERSE);
            assert(ierr == PETSC_SUCCESS);

            ierr = VecGhostUpdateEnd(m_data, INSERT_VALUES, SCATTER_REVERSE);
            assert(ierr == PETSC_SUCCESS);

            m_read.acquired = false;
          }
        }
        else
        {
          assert(false);
        }
        (void) ierr;
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

      const ArrayRead& getArrayRead() const
      {
        return m_read;
      }

      const ArrayWrite& getArrayWrite() const
      {
        return m_write;
      }

      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope&) const
      {
        return std::nullopt;
      }

    private:
      DataType m_data;
      size_t m_begin, m_end;
      GhostBimap m_ghosts;

      mutable ArrayRead m_read;
      mutable ArrayWrite m_write;
  };
}

namespace Rodin::PETSc::Variational
{
  template <class FES>
  class GridFunction
    : public Rodin::Variational::GridFunction<FES, ::Vec>
  {
    public:
      using Parent = Rodin::Variational::GridFunction<FES, ::Vec>;
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
