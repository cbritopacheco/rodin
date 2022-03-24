/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_RODININTEGRATION_MMG_CAST_H
#define RODIN_RODININTEGRATION_MMG_CAST_H

#include "Rodin/Cast.h"

#include "Rodin/Mesh/ForwardDecls.h"
#include "Rodin/Variational/ForwardDecls.h"
#include "ForwardDecls.h"

namespace Rodin
{
   // ---- mmg2d -------------------------------------------------------------
   /**
    * MMG::Mesh2D -> Rodin::Mesh
    */
   template <>
   template <>
   Rodin::Mesh
   Cast<External::MMG::Mesh2D>::to<Rodin::Mesh>() const;

   /**
    * Rodin::Mesh -> MMG::Mesh2D
    */
   template <>
   template <>
   External::MMG::Mesh2D
   Cast<Rodin::Mesh>::to<External::MMG::Mesh2D>() const;

   /**
    * MMG::ScalarSolution2D -> Rodin::Variational::IncompleteGridFunction
    */
   template <>
   template <>
   Variational::IncompleteGridFunction
   Cast<External::MMG::IncompleteScalarSolution2D>
   ::to<Variational::IncompleteGridFunction>() const;

   /**
    * Rodin::Variational::GridFunction<H1> -> MMG::IncompleteScalarSolution2D
    */
   template <>
   template <>
   External::MMG::IncompleteScalarSolution2D
   Cast<Variational::GridFunction<Variational::H1>>
   ::to<External::MMG::IncompleteScalarSolution2D>() const;

   /**
    * Rodin::Variational::GridFunction<H1> -> MMG::IncompleteVectorSolution2D
    */
   template <>
   template <>
   External::MMG::IncompleteVectorSolution2D
   Cast<Variational::GridFunction<Variational::H1>>
   ::to<External::MMG::IncompleteVectorSolution2D>() const;

   /**
    * Rodin::Variational::IncompleteGridFunction -> MMG::ScalarSolution2D
    */
   template <>
   template <>
   Variational::IncompleteGridFunction
   Cast<External::MMG::ScalarSolution2D>
   ::to<Variational::IncompleteGridFunction>() const;

   /**
    * Rodin::Variational::IncompleteGridFunction -> MMG::IncompleteScalarSolution2D
    */
   template <>
   template <>
   External::MMG::IncompleteScalarSolution2D
   Cast<Variational::IncompleteGridFunction>
   ::to<External::MMG::IncompleteScalarSolution2D>() const;

   // ---- mmg3d -------------------------------------------------------------
   /**
    * MMG::Mesh3D -> Rodin::Mesh
    */
   template <>
   template <>
   Rodin::Mesh
   Cast<External::MMG::Mesh3D>::to<Rodin::Mesh>() const;

   /**
    * Rodin::Mesh -> MMG::Mesh3D
    */
   template <>
   template <>
   External::MMG::Mesh3D
   Cast<Rodin::Mesh>::to<External::MMG::Mesh3D>() const;

   /**
    * Rodin::Variational::GridFunction<H1> -> MMG::IncompleteScalarSolution3D
    */
   template <>
   template <>
   External::MMG::IncompleteScalarSolution3D
   Cast<Variational::GridFunction<Variational::H1>>
   ::to<External::MMG::IncompleteScalarSolution3D>() const;

   /**
    * Rodin::Variational::IncompleteGridFunction -> MMG::IncompleteScalarSolution3D
    */
   template <>
   template <>
   External::MMG::IncompleteScalarSolution3D
   Cast<Variational::IncompleteGridFunction>
   ::to<External::MMG::IncompleteScalarSolution3D>() const;

   // ---- mmgs --------------------------------------------------------------
   /**
    * MMG::SurfaceMesh -> Rodin::Mesh.
    */
   template <>
   template <>
   Rodin::Mesh
   Cast<External::MMG::MeshS>::to<Rodin::Mesh>() const;

   /**
    * Rodin::Mesh -> MMG::SurfaceMesh
    */
   template <>
   template <>
   External::MMG::MeshS
   Cast<Rodin::Mesh>::to<External::MMG::MeshS>() const;
}


#endif
