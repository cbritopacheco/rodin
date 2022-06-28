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
#include "Rodin/Variational/H1.h"
#include "Rodin/Variational/GridFunction.h"
#include "Rodin/Variational/FiniteElementSpace.h"

#include "Utility.h"
#include "ScalarSolution2D.h"
#include "ScalarSolution3D.h"
#include "ScalarSolutionS.h"

#include "VectorSolution2D.h"
#include "VectorSolution3D.h"
#include "VectorSolutionS.h"

#include "Common.h"
#include "ForwardDecls.h"

namespace Rodin
{
   // ---- mmg2d -------------------------------------------------------------
   /**
    * MMG::Mesh2D -> Rodin::Mesh<>
    */
   template <>
   template <>
   Rodin::Mesh<Traits::Serial>
   Cast<External::MMG::Mesh2D>::to<Rodin::Mesh<>>() const;

   /**
    * Rodin::Mesh<Traits::Serial> -> MMG::Mesh2D
    */
   template <>
   template <>
   External::MMG::Mesh2D
   Cast<Rodin::Mesh<Traits::Serial>>::to<External::MMG::Mesh2D>() const;

   template <>
   template <>
   External::MMG::Mesh2D
   Cast<Rodin::SubMesh<Traits::Serial>>::to<External::MMG::Mesh2D>() const;

   /**
    * Rodin::Variational::GridFunction<FEC, Traits::Serial> -> MMG::IncompleteScalarSolution2D
    */
   template <class FEC>
   class ADLCaster<
      Variational::GridFunction<FEC, Traits::Serial>, External::MMG::ScalarSolution2D>
   {
      public:
         ADLCaster(External::MMG::Mesh2D& mesh)
            : m_mesh(mesh)
         {}

         External::MMG::ScalarSolution2D cast(
               const Variational::GridFunction<FEC, Traits::Serial>& src)
         {
            assert(src.getFiniteElementSpace().getVectorDimension() == 1);
            External::MMG::ScalarSolution2D res(m_mesh);
            External::MMG::Rodin_GridFunction_To_MMG5_Sol(src, res.getHandle());
            return res;
         }

      private:
         External::MMG::Mesh2D& m_mesh;
   };

   /**
    * Rodin::Variational::GridFunction<FEC, Traits::Serial> -> MMG::VectorSolution2D
    */
   template <class FEC>
   class ADLCaster<
      Variational::GridFunction<FEC, Traits::Serial>, External::MMG::VectorSolution2D>
   {
      public:
         ADLCaster(External::MMG::Mesh2D& mesh)
            : m_mesh(mesh)
         {}

         External::MMG::VectorSolution2D cast(
               const Variational::GridFunction<FEC, Traits::Serial>& src)
         {
            assert(src.getFiniteElementSpace().getVectorDimension() == 2);
            External::MMG::VectorSolution2D res(m_mesh);
            External::MMG::Rodin_GridFunction_To_MMG5_Sol(src, res.getHandle());
            return res;
         }

      private:
         External::MMG::Mesh2D& m_mesh;
   };

   /**
    * MMG::ScalarSolution2D -> Rodin::Variational::GridFunction<Variational::H1, Traits::Serial>
    */
   template <class FEC>
   class ADLCaster<External::MMG::ScalarSolution2D, Variational::GridFunction<FEC, Traits::Serial>>
   {
      public:
         ADLCaster(Variational::FiniteElementSpace<FEC>& mesh)
            : m_fes(mesh)
         {}

         Variational::GridFunction<FEC, Traits::Serial> cast(
               const External::MMG::ScalarSolution2D& src)
         {
            Variational::GridFunction<FEC> res(m_fes);
            External::MMG::MMG5_Sol_To_Rodin_GridFunction(src.getHandle(), res);
            return res;
         }
      private:
         Variational::FiniteElementSpace<FEC>& m_fes;
   };

   // ---- mmg3d -------------------------------------------------------------
   /**
    * MMG::Mesh3D -> Rodin::Mesh<>
    */
   template <>
   template <>
   Rodin::Mesh<>
   Cast<External::MMG::Mesh3D>::to<Rodin::Mesh<Traits::Serial>>() const;

   /**
    * Rodin::Mesh<> -> MMG::Mesh3D
    */
   template <>
   template <>
   External::MMG::Mesh3D
   Cast<Rodin::Mesh<Traits::Serial>>::to<External::MMG::Mesh3D>() const;

   template <>
   template <>
   External::MMG::Mesh3D
   Cast<Rodin::SubMesh<Traits::Serial>>::to<External::MMG::Mesh3D>() const;

   /**
    * MMG::ScalarSolution3D -> Rodin::Variational::GridFunction<Variational::H1, Traits::Serial>
    */
   template <class FEC>
   class ADLCaster<External::MMG::ScalarSolution3D, Variational::GridFunction<FEC, Traits::Serial>>
   {
      public:
         ADLCaster(Variational::FiniteElementSpace<FEC>& mesh)
            : m_fes(mesh)
         {}

         Variational::GridFunction<FEC, Traits::Serial> cast(
               const External::MMG::ScalarSolution3D& src)
         {
            Variational::GridFunction<FEC> res(m_fes);
            External::MMG::MMG5_Sol_To_Rodin_GridFunction(src.getHandle(), res);
            return res;
         }
      private:
         Variational::FiniteElementSpace<FEC>& m_fes;
   };

   template <class FEC>
   class ADLCaster<
      Variational::GridFunction<FEC, Traits::Serial>, External::MMG::ScalarSolution3D>
   {
      public:
         ADLCaster(External::MMG::Mesh3D& mesh)
            : m_mesh(mesh)
         {}

         External::MMG::ScalarSolution3D cast(
               const Variational::GridFunction<FEC, Traits::Serial>& src)
         {
            assert(src.getFiniteElementSpace().getVectorDimension() == 1);
            External::MMG::ScalarSolution3D res(m_mesh);
            External::MMG::Rodin_GridFunction_To_MMG5_Sol(src, res.getHandle());
            return res;
         }

      private:
         External::MMG::Mesh3D& m_mesh;
   };

   template <class FEC>
   class ADLCaster<
      Variational::GridFunction<FEC, Traits::Serial>, External::MMG::VectorSolution3D>
   {
      public:
         ADLCaster(External::MMG::Mesh3D& mesh)
            : m_mesh(mesh)
         {}

         External::MMG::VectorSolution3D cast(
               const Variational::GridFunction<FEC, Traits::Serial>& src)
         {
            assert(src.getFiniteElementSpace().getVectorDimension() == 3);
            External::MMG::VectorSolution3D res(m_mesh);
            External::MMG::Rodin_GridFunction_To_MMG5_Sol(src, res.getHandle());
            return res;
         }

      private:
         External::MMG::Mesh3D& m_mesh;
   };

   // ---- mmgs --------------------------------------------------------------
   /**
    * MMG::SurfaceMesh -> Rodin::Mesh<>.
    */
   template <>
   template <>
   Rodin::Mesh<>
   Cast<External::MMG::MeshS>::to<Rodin::Mesh<>>() const;

   /**
    * Rodin::Mesh<> -> MMG::SurfaceMesh
    */
   template <>
   template <>
   External::MMG::MeshS
   Cast<Rodin::Mesh<Traits::Serial>>::to<External::MMG::MeshS>() const;

   template <>
   template <>
   External::MMG::MeshS
   Cast<Rodin::SubMesh<Traits::Serial>>::to<External::MMG::MeshS>() const;

   /**
    * MMG::ScalarSolutionS -> Rodin::Variational::GridFunction<Variational::H1, Traits::Serial>
    */
   template <class FEC>
   class ADLCaster<External::MMG::ScalarSolutionS, Variational::GridFunction<FEC, Traits::Serial>>
   {
      public:
         ADLCaster(Variational::FiniteElementSpace<FEC>& mesh)
            : m_fes(mesh)
         {}

         Variational::GridFunction<FEC, Traits::Serial> cast(
               const External::MMG::ScalarSolutionS& src)
         {
            Variational::GridFunction<FEC> res(m_fes);
            External::MMG::MMG5_Sol_To_Rodin_GridFunction(src.getHandle(), res);
            return res;
         }
      private:
         Variational::FiniteElementSpace<FEC>& m_fes;
   };

   template <class FEC>
   class ADLCaster<
      Variational::GridFunction<FEC, Traits::Serial>, External::MMG::ScalarSolutionS>
   {
      public:
         ADLCaster(External::MMG::MeshS& mesh)
            : m_mesh(mesh)
         {}

         External::MMG::ScalarSolutionS cast(
               const Variational::GridFunction<FEC, Traits::Serial>& src)
         {
            assert(src.getFiniteElementSpace().getVectorDimension() == 1);
            External::MMG::ScalarSolutionS res(m_mesh);
            External::MMG::Rodin_GridFunction_To_MMG5_Sol(src, res.getHandle());
            return res;
         }

      private:
         External::MMG::MeshS& m_mesh;
   };

   template <class FEC>
   class ADLCaster<
      Variational::GridFunction<FEC, Traits::Serial>, External::MMG::VectorSolutionS>
   {
      public:
         ADLCaster(External::MMG::MeshS& mesh)
            : m_mesh(mesh)
         {}

         External::MMG::VectorSolutionS cast(
               const Variational::GridFunction<FEC, Traits::Serial>& src)
         {
            assert(src.getFiniteElementSpace().getVectorDimension() == 3);
            External::MMG::VectorSolutionS res(m_mesh);
            External::MMG::Rodin_GridFunction_To_MMG5_Sol(src, res.getHandle());
            return res;
         }

      private:
         External::MMG::MeshS& m_mesh;
   };
}


#endif
