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
          m_begin(other.m_begin),
          m_end(other.m_end),
          m_ghosts(other.m_ghosts),
          m_read{.acquired = false, .ghost = PETSC_NULLPTR, .raw = PETSC_NULLPTR},
          m_write{.acquired = false, .ghost = PETSC_NULLPTR, .raw = PETSC_NULLPTR}
      {
        PetscErrorCode ierr = VecDuplicate(other.m_data, &m_data);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecCopy(other.m_data, m_data);
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;
      }

      GridFunction(GridFunction&& other) noexcept
        : Parent(std::move(other)),
          m_data(std::exchange(other.m_data, PETSC_NULLPTR)),
          m_begin(std::exchange(other.m_begin, 0)),
          m_end(std::exchange(other.m_end, 0)),
          m_ghosts(std::move(other.m_ghosts)),
          m_read{.acquired = false, .ghost = PETSC_NULLPTR, .raw = PETSC_NULLPTR},
          m_write{.acquired = false, .ghost = PETSC_NULLPTR, .raw = PETSC_NULLPTR}
      {
        other.m_read = {.acquired = false, .ghost = PETSC_NULLPTR, .raw = PETSC_NULLPTR};
        other.m_write = {.acquired = false, .ghost = PETSC_NULLPTR, .raw = PETSC_NULLPTR};
      }

      GridFunction& operator=(const GridFunction& other)
      {
        if (this == &other)
          return *this;

        assert(&this->getFiniteElementSpace() == &other.getFiniteElementSpace());

        this->release();

        m_begin = other.m_begin;
        m_end = other.m_end;
        m_ghosts = other.m_ghosts;

        m_read = {.acquired = false, .ghost = PETSC_NULLPTR, .raw = PETSC_NULLPTR};
        m_write = {.acquired = false, .ghost = PETSC_NULLPTR, .raw = PETSC_NULLPTR};

        PetscErrorCode ierr = VecDuplicate(other.m_data, &m_data);
        assert(ierr == PETSC_SUCCESS);

        ierr = VecCopy(other.m_data, m_data);
        assert(ierr == PETSC_SUCCESS);

        (void) ierr;
        return *this;
      }

      GridFunction& operator=(GridFunction&& other) noexcept
      {
        if (this == &other)
          return *this;

        assert(&this->getFiniteElementSpace() == &other.getFiniteElementSpace());

        this->release();

        m_data = std::exchange(other.m_data, PETSC_NULLPTR);
        m_begin = std::exchange(other.m_begin, 0);
        m_end = std::exchange(other.m_end, 0);
        m_ghosts = std::move(other.m_ghosts);

        m_read = {.acquired = false, .ghost = PETSC_NULLPTR, .raw = PETSC_NULLPTR};
        m_write = {.acquired = false, .ghost = PETSC_NULLPTR, .raw = PETSC_NULLPTR};

        other.m_read = {.acquired = false, .ghost = PETSC_NULLPTR, .raw = PETSC_NULLPTR};
        other.m_write = {.acquired = false, .ghost = PETSC_NULLPTR, .raw = PETSC_NULLPTR};

        return *this;
      }

      virtual ~GridFunction()
      {
        this->release();
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

      template <class Function, class Pred>
      GridFunction& project(const Geometry::Region& region, const Function& v, const Pred& pred)
      {
        const auto& fes = this->getFiniteElementSpace();
        const auto& mesh = fes.getMesh();

        Geometry::PolytopeIterator it;
        switch (region)
        {
          case Geometry::Region::Cells:
            it = mesh.getCell();
            break;
          case Geometry::Region::Faces:
            it = mesh.getFace();
            break;
          case Geometry::Region::Boundary:
            it = mesh.getBoundary();
            break;
          case Geometry::Region::Interface:
            it = mesh.getInterface();
            break;
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

      decltype(auto) operator()(const Geometry::Point& p) const
      {
        return this->getValue(p);
      }

      decltype(auto) getValue(const Geometry::Point& p) const
      {
        if constexpr (std::is_same_v<FESMeshContextType, Context::Local>)
        {
          return Parent::getValue(p);
        }
        else if constexpr (std::is_same_v<FESMeshContextType, Context::MPI>)
        {
          static thread_local RangeType s_out;

          const auto& fes = this->getFiniteElementSpace();
          const auto& mesh = fes.getMesh();
          const auto& shard = mesh.getShard();
          const auto& polytope = p.getPolytope();
          const auto& polytopeMesh = polytope.getMesh();

          if (polytopeMesh == shard)
          {
            this->interpolate(s_out, p);
            return s_out;
          }

          if (const auto inclusion = shard.inclusion(p))
          {
            this->interpolate(s_out, *inclusion);
            return s_out;
          }

          if (shard.isSubMesh())
          {
            const auto& submesh = shard.asSubMesh();
            if (const auto restriction = submesh.restriction(p))
            {
              this->interpolate(s_out, *restriction);
              return s_out;
            }
          }

          Alert::MemberFunctionException(*this, __func__)
            << "Point does not belong to the PETSc GridFunction shard."
            << Alert::Raise;

          assert(false);
          return s_out;
        }
        else
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Unsupported mesh context type for PETSc GridFunction."
            << Alert::Raise;
          return Parent::getValue(p);
        }
      }

      void release()
      {
        PetscErrorCode ierr;

        if (m_read.acquired)
        {
          assert(m_read.raw);
          if constexpr (std::is_same_v<FESMeshContextType, Context::Local>)
          {
            ierr = VecRestoreArrayRead(m_data, &m_read.raw);
            assert(ierr == PETSC_SUCCESS);
          }
          else if constexpr (std::is_same_v<FESMeshContextType, Context::MPI>)
          {
            ierr = VecRestoreArrayRead(m_read.ghost, &m_read.raw);
            assert(ierr == PETSC_SUCCESS);

            ierr = VecGhostRestoreLocalForm(m_data, &m_read.ghost);
            assert(ierr == PETSC_SUCCESS);
          }

          m_read.acquired = false;
          m_read.raw = PETSC_NULLPTR;
          m_read.ghost = PETSC_NULLPTR;
        }

        if (m_write.acquired)
        {
          assert(m_write.raw);
          if constexpr (std::is_same_v<FESMeshContextType, Context::Local>)
          {
            ierr = VecRestoreArrayWrite(m_data, &m_write.raw);
            assert(ierr == PETSC_SUCCESS);
          }
          else if constexpr (std::is_same_v<FESMeshContextType, Context::MPI>)
          {
            ierr = VecRestoreArrayWrite(m_write.ghost, &m_write.raw);
            assert(ierr == PETSC_SUCCESS);

            ierr = VecGhostRestoreLocalForm(m_data, &m_write.ghost);
            assert(ierr == PETSC_SUCCESS);

            ierr = VecGhostUpdateBegin(m_data, INSERT_VALUES, SCATTER_REVERSE);
            assert(ierr == PETSC_SUCCESS);

            ierr = VecGhostUpdateEnd(m_data, INSERT_VALUES, SCATTER_REVERSE);
            assert(ierr == PETSC_SUCCESS);
          }

          m_write.acquired = false;
          m_write.raw = PETSC_NULLPTR;
          m_write.ghost = PETSC_NULLPTR;
        }

        if (m_data)
        {
          ierr = VecDestroy(&m_data);
          assert(ierr == PETSC_SUCCESS);
          m_data = PETSC_NULLPTR;
        }

        (void) ierr;
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
      using Parent::operator();
      using Parent::getValue;
      using Parent::operator=;
      using Parent::operator+=;
      using Parent::operator-=;
      using Parent::operator*=;
      using Parent::operator/=;
  };

  template <class FES>
  GridFunction(const FES&) -> GridFunction<FES>;
}

namespace Rodin::FormLanguage
{
  template <class FES>
  struct Traits<PETSc::Variational::GridFunction<FES>>
  {
    using FESType = FES;
    using DataType = ::Vec;
  };
}

#endif
