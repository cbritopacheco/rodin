/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_RODININTEGRATION_MMG_FORWARDDECLS_H
#define RODIN_RODININTEGRATION_MMG_FORWARDDECLS_H

namespace Rodin::External::MMG
{
   template <int Dimension, class Derived>
   class Mesh;

   class Mesh2D;

   template <int Dimension, class Derived>
   class ScalarSolution;

   class ScalarSolution2D;

   template <class MeshType, class Derived>
   class MeshAdaptor;

   class MeshAdaptor2D;

   template <class MeshType, class Derived>
   class MeshOptimizer;

   class MeshOptimizer2D;
}

#endif
