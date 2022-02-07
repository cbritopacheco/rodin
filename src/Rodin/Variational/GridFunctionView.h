/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_GRIDFUNCTIONVIEW_H
#define RODIN_VARIATIONAL_GRIDFUNCTIONVIEW_H

#include <memory>
#include <variant>
#include <utility>
#include <type_traits>

#include "ForwardDecls.h"

namespace Rodin::Variational
{
  class GridFunctionView
  {
    public:
      GridFunctionView(GridFunctionBase& gf);

      template <class T>
      std::enable_if_t<
        std::is_base_of_v<GridFunctionIndexBase, T>, GridFunctionView&>
      setIndex(T&& idx)
      {
         m_idx.reset(new T(std::forward<T>(idx)));
         return *this;
      }

      GridFunctionView& operator=(double v);

    private:
       GridFunctionBase& m_gf;
       std::unique_ptr<GridFunctionIndexBase> m_idx;
  };
}

#endif

