/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_GRIDFUNCTIONINDEX_H
#define RODIN_VARIATIONAL_GRIDFUNCTIONINDEX_H

#include <deque>

#include "Rodin/Core.h"

#include "H1.h"
#include "ForwardDecls.h"

namespace Rodin
{
   namespace Variational
   {
      class GridFunctionIndexBase
      {
         public:
            virtual ~GridFunctionIndexBase() = default;
            virtual std::vector<int> getIndices() const = 0;
      };

      template <>
      class GridFunctionIndex<bool> : public GridFunctionIndexBase
      {
         public:
            GridFunctionIndex(const std::deque<bool>& idx);

            std::vector<int> getIndices() const override;

            GridFunctionIndex operator||(const GridFunctionIndex& other);

         private:
            std::deque<bool> m_idx;
      };
   }

   Variational::GridFunctionIndex<bool>
   isNaN(const Variational::GridFunction<Variational::H1>& v);

   Variational::GridFunctionIndex<bool>
   isInf(const Variational::GridFunction<Variational::H1>& v);
}

#endif
