#ifndef RODIN_PETSC_VARIATIONAL_P1_P1_H
#define RODIN_PETSC_VARIATIONAL_P1_P1_H

#include <petsc.h>

#include "Rodin/Math/Vector.h"
#include "Rodin/MPI/Geometry/Mesh.h"
#include "Rodin/MPI/Variational/P1/P1.h"

namespace Rodin::PETSc::Variational
{
  template <class Range, class Mesh>
  class P1;

  template <class Range>
  class P1<Range, Geometry::Mesh<Context::Local>>
    : public Rodin::Variational::P1<Range, Geometry::Mesh<Context::Local>>
  {
    public:
      using Parent = Rodin::Variational::P1<Range, Geometry::Mesh<Context::Local>>;

      using Parent::Parent;
  };

  template <class Range>
  class P1<Range, Geometry::Mesh<Context::MPI>>
    : public Rodin::Variational::P1<Range, Geometry::Mesh<Context::MPI>>
  {
    public:
      using MeshType = Geometry::Mesh<Context::MPI>;

      struct GhostBimap
      {
        std::vector<PetscInt> left;
        FlatMap<PetscInt, PetscInt> right;
      };

      P1(const MeshType& mesh)
      {
        const auto& shard = this->getShard();
        size_t begin, end;
        this->getOwnershipRange(begin, end);
        Index offset = 0;
        for (size_t i = 0; i < shard.getSize(); ++i)
        {
          const Index global = this->getGlobalIndex(i);
          if (global < begin || end <= global)
          {
            auto [it, inserted] = m_ghosts.right.insert({ global, offset });
            assert(inserted);
            m_ghosts.left.push_back(offset);
            offset++;
          }
        }
      }

      P1(const MeshType& mesh, size_t vdim)
      {
        assert(false);
      }

      virtual ~P1() = default;

      const auto& getGhosts() const
      {
        return m_ghosts;
      }

    private:
      GhostBimap m_ghosts;
  };

  template <class Mesh>
  P1(const Mesh&) -> P1<Real, Mesh>;

  template <class Mesh>
  P1(const Mesh&, size_t) -> P1<Rodin::Math::Vector<Real>, Mesh>;

  // template <class Mesh>
  // P1(const Mesh&, size_t) -> P1<Math::Vector<Real>, Mesh>;
}

#endif
