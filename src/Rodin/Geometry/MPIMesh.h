/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_MPI_MESH_H
#define RODIN_GEOMETRY_MPI_MESH_H

#ifdef RODIN_USE_MPI

#include "Rodin/Configure.h"
#include "Rodin/Context/MPI.h"

#include "Mesh.h"

namespace Rodin::Geometry
{
  using MPIMesh = Mesh<Context::MPI>;

  template <>
  class Mesh<Context::MPI> : public MeshBase
  {
    public:
      class Fragment : public Mesh<Context::Sequential>
      {
        public:
          using Parent = Mesh<Context::Sequential>;
          using Parent::Parent;
      };

      Index getGlobalIndex(const std::pair<size_t, Index>& p, Index fragmentId);

      const Fragment& getFragment(Index i) const
      {
        return m_fragments[i];
      }

    private:
      std::vector<Fragment> m_fragments;
  };
}

#endif
#endif
