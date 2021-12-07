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
   template <>
   Rodin::Mesh
   Cast<External::MMG::Mesh2D>::to<Rodin::Mesh>() const;

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
   template <>
   External::MMG::Mesh2D
   Cast<Rodin::Mesh>::to<External::MMG::Mesh2D>() const;

   /**
    * @brief Specialization for converting from External::MMG::ScalarSolution2D to
    * Rodin::Variational::GridFunction.
    *
    * @note This is a lossy cast. Data from the old object that has no direct
    * correspondence will not be present in the new object.
    *
    * @todo Which fields are not compatible?
    */
   template <>
   template <>
   Variational::GridFunction<>
   Cast<External::MMG::ScalarSolution2D<true>>
   ::to<Variational::GridFunction<>>() const;

   /**
    * @brief Specialization for converting from External::MMG::ScalarSolution2D to
    * Rodin::Variational::GridFunction.
    *
    * @note This is a lossy cast. Data from the old object that has no direct
    * correspondence will not be present in the new object.
    *
    * @todo Which fields are not compatible?
    */
   template <>
   template <>
   Variational::GridFunction<>
   Cast<External::MMG::ScalarSolution2D<false>>
   ::to<Variational::GridFunction<>>() const;

   /**
    * @brief Specialization for converting from
    * Rodin::Variational::GridFunction<H1> to
    * External::MMG::ScalarSolution2D<false>.
    *
    * @note This is a lossy cast. Data from the old object that has no direct
    * correspondence will not be present in the new object.
    *
    * @todo Which fields are not compatible?
    */
   template <>
   template <>
   External::MMG::ScalarSolution2D<false>
   Cast<Variational::GridFunction<Variational::H1>>
   ::to<External::MMG::ScalarSolution2D<false>>() const;

   /**
    * @brief Specialization for converting from
    * Rodin::Variational::GridFunction<> to
    * External::MMG::ScalarSolution2D<false>.
    *
    * @note This is a lossy cast. Data from the old object that has no direct
    * correspondence will not be present in the new object.
    *
    * @todo Which fields are not compatible?
    */
   template <>
   template <>
   External::MMG::ScalarSolution2D<false>
   Cast<Variational::GridFunction<>>
   ::to<External::MMG::ScalarSolution2D<false>>() const;
}


#endif
