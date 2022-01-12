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
   template <int Dimension>
   class Mesh;

   class Mesh2D;

   class MeshAdaptor2D;

   class MeshOptimizer2D;

   template <int Dimension>
   class ScalarSolution;

   template <bool HasMesh = true>
   class ScalarSolution2D;

   template <>
   class ScalarSolution2D<true>;

   template <>
   class ScalarSolution2D<false>;

   class Distancer2D;
   class Redistancer2D;
}

#endif
