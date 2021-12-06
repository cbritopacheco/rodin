/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_RODININTEGRATION_MMG_CAST_H
#define RODIN_RODININTEGRATION_MMG_CAST_H

#include "Rodin/Cast.h"

#include "ForwardDecls.h"

#include "Rodin/Mesh.h"
#include "Rodin/Variational.h"

#include "Mesh2D.h"

namespace Rodin
{
   /**
    * @brief Specialization for converting from External::MMG::Mesh2D to Rodin::Mesh.
    *
    * @note This is a lossy cast. Data from the old object that has no direct
    * correspondence will not be present in the new object.
    *
    * @todo Which fields are not compatible?
    */
   template <>
   class Cast<External::MMG::Mesh2D, Rodin::Mesh>
   {
      public:
         /**
          * @brief Performs the cast.
          *
          * @returns New object.
          */
         Rodin::Mesh cast(const External::MMG::Mesh2D& mesh) const;
   };


   /**
    * @brief Specialization for converting from Rodin::Mesh to
    * External::MMG:Mesh2D.
    *
    * @note This is a lossy cast. Data from the old object that has no direct
    * correspondence will not be present in the new object.
    *
    * @todo Which fields are not compatible?
    */
   template <>
   class Cast<Rodin::Mesh, External::MMG::Mesh2D>
   {
      public:
         /**
          * @brief Performs the cast.
          *
          * @returns New object.
          */
         External::MMG::Mesh2D cast(const Rodin::Mesh& mesh) const;
   };

   /**
    * @brief Specialization for converting from External::MMG::ScalarSolution2D to
    * Rodin::Variational::GridFunction.
    *
    * @tparam FEC Finite element collection to which the
    * Variational::GridFunction will belong to.
    *
    * @note This is a lossy cast. Data from the old object that has no direct
    * correspondence will not be present in the new object.
    *
    * @todo Which fields are not compatible?
    */
   template <class FEC, bool HasMesh>
   class Cast<External::MMG::ScalarSolution2D<HasMesh>,
         Variational::GridFunction<FEC>>
   {
      public:
         /**
          * @brief Initializes the cast object with the finite element space
          * which is required to construct the grid function.
          * @param[in] fes Finite element space which the grid function will
          * lie in.
          */
         Cast(Variational::FiniteElementSpace<FEC>& fes);

         /**
          * @brief Performs the cast.
          *
          * @returns New object.
          */
         Variational::GridFunction<FEC>
         cast(const External::MMG::ScalarSolution2D<HasMesh>& sol) const;

      private:
         Variational::FiniteElementSpace<FEC>& m_fes;
   };

   // /**
   //  * @brief Specialization for converting from
   //  * Rodin::Variational::GridFunction to External::MMG::ScalarSolution2D.
   //  */
   // template <class FEC>
   // class Cast<Variational::GridFunction<FEC>, External::MMG::ScalarSolution2D>
   // {
   //    public:
   //       /**
   //        * @brief Performs the cast.
   //        *
   //        * @returns New object.
   //        */
   //       External::MMG::ScalarSolution2D
   //       cast(const Variational::GridFunction<FEC>& gf) const;
   // };
}

#include "Cast.hpp"

#endif
