/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_OPENMPASSEMBLY_H
#define RODIN_VARIATIONAL_OPENMPASSEMBLY_H

#include <variant>
#include <ostream>

#include <mfem.hpp>

#include "AssemblyBase.h"

namespace Rodin::Variational::Assembly
{
   template <>
   class NativeAssembly<BilinearFormBase<mfem::SparseMatrix>>
      : public Assembly<BilinearFormBase<mfem::SparseMatrix>>
   {
      public:
         using Parent = Assembly<BilinearFormBase<mfem::SparseMatrix>>;
         using OperatorType = mfem::SparseMatrix;

         NativeAssembly() = default;

         NativeAssembly(const NativeAssembly& other)
            : Parent(other)
         {}

         NativeAssembly(NativeAssembly&& other)
            : Parent(std::move(other))
         {}

         OperatorType execute(const Input& input) const override;

         Assembly* copy() const noexcept override
         {
            return new NativeAssembly(*this);
         }
   };
}

#endif

