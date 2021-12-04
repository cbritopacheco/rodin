/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_GRIDFUNCTION_H
#define RODIN_VARIATIONAL_GRIDFUNCTION_H

#include <mfem.hpp>

#include <functional>

#include "ForwardDecls.h"

namespace Rodin::Variational
{
   class GridFunctionBase
   {
      public:
         virtual mfem::GridFunction& getHandle() = 0;
         virtual const mfem::GridFunction& getHandle() const = 0;
   };

   /**
    * @brief A grid function is a function which is defined on an unstructured
    * grid, i.e. a triangular mesh.
    *
    * @tparam FEC Finite element collection to which the function belongs.
    *
    * @note Note that the FEC template parameter is typically inferred when
    * initializing the grid function, hence it is not necessary to explicitly
    * specify it.
    */
   template <class FEC>
   class GridFunction : public GridFunctionBase
   {
      public:
         /**
          * @brief Constructs a grid function on a finite element space.
          * @param[in] fes Finite element space to which the function belongs
          * to.
          */
         GridFunction(FiniteElementSpace<FEC>& fes)
            :  m_fes(fes),
               m_gf(&fes.getFES())
         {}

         /**
          * @brief Copies the grid function.
          * @param[in] other Other grid function to copy.
          */
         GridFunction(const GridFunction& other)
            :  m_fes(other.m_fes),
               m_gf(other.m_gf)
         {}

         /**
          * @brief Returns the finite element space to which the function
          * belongs to.
          * @returns Finite element space to which the function belongs to.
          */
         FiniteElementSpace<FEC>& getFiniteElementSpace()
         {
            return m_fes;
         }

         mfem::GridFunction& getHandle() override
         {
            return m_gf;
         }

         const mfem::GridFunction& getHandle() const override
         {
            return m_gf;
         }

         void save(const std::string& filename)
         {
            m_gf.Save(filename.c_str());
         }

      private:
         std::reference_wrapper<FiniteElementSpace<FEC>> m_fes;
         mfem::GridFunction m_gf;
   };
}

#endif
