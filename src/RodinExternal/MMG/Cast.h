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

#include "ScalarSolution2D.h"

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
   template <bool HasMesh>
   class Cast<External::MMG::ScalarSolution2D<HasMesh>>
   {
      public:
         Cast(const External::MMG::ScalarSolution2D<HasMesh>& from)
            : m_from(from)
         {}

         const auto& from() const
         {
            return m_from;
         }

         template <class To>
         To to() const;

         template <>
         Variational::GridFunction<> to<Variational::GridFunction<>>() const
         {
           auto& sol = from();
           MMG5_pSol mmgSol = sol.getHandle();
           assert(mmgSol->type == MMG5_Scalar);
           Variational::GridFunction<> res;
           double* data = new double[mmgSol->np];
           // MMG5_pSol->m is 1 indexed. We must start at m + 1 and finish at m
           // + np + 1.
           std::copy(mmgSol->m + 1, mmgSol->m + mmgSol->np + 1, data);
           res.getHandle().SetDataAndSize(data, mmgSol->np);
           res.getHandle().MakeDataOwner();
           return res;
         }

      private:
         const External::MMG::ScalarSolution2D<HasMesh>& m_from;
   };

   /**
    * @brief Specialization for converting from
    * Rodin::Variational::GridFunction<FEC> to
    * External::MMG::ScalarSolution2D<false>.
    *
    * @note This is a lossy cast. Data from the old object that has no direct
    * correspondence will not be present in the new object.
    *
    * @todo Which fields are not compatible?
    */
   template <class FEC>
   class Cast<Variational::GridFunction<FEC>>
   {
      public:
         Cast(const Variational::GridFunction<FEC>& from)
            : m_from(from)
         {}

         const auto& from() const
         {
            return m_from;
         }

         template <class To>
         To to() const;

         template <>
         External::MMG::ScalarSolution2D<false>
         to<External::MMG::ScalarSolution2D<false>>()
         const
         {
           auto& gf = from();
           auto [data, size] = gf.getData();
           if (!size)
             return External::MMG::ScalarSolution2D<false>();
           else
           {
             External::MMG::ScalarSolution2D<false> res(size);
             // MMG5_pSol->m is 1 indexed. We must start at m + 1 and finish at
             // m + size + 1.
             std::copy(data, data + size, res.getHandle()->m + 1);
             return res;
           }
         }

      private:
         const Variational::GridFunction<FEC>& m_from;
   };

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
