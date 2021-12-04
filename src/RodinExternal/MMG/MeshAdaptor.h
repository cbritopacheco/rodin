/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MESH_MESHADAPTOR_H
#define RODIN_MESH_MESHADAPTOR_H

namespace Rodin::External::MMG
{
  template <class MeshType, class Derived>
  class MeshAdaptor
  {
    public:
      void adapt(MeshType& mesh)
      {
        Derived::optimize(mesh);
      }
  };
}

#endif
