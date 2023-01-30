/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2023.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ARRAY_H
#define RODIN_ARRAY_H

#include <Eigen/Core>
#include <Eigen/Dense>

namespace Rodin
{
  template <class Scalar>
  class Array : public Eigen::ArrayX<Scalar>
  {
    public:
      using Parent = Eigen::ArrayX<Scalar>;

      Array(std::initializer_list<typename Parent::Scalar> l)
        : Parent(l.size())
      {
        std::copy(l.begin(), l.end(), Parent::begin());
      }

      template <class ... Args>
      Array(Args&&... args)
        : Parent(std::forward<Args>(args)...)
      {}

      template <class ... Args>
      Array& operator=(Args&&... args)
      {
         this->Parent::operator=(std::forward<Args>(args)...);
         return *this;
      }
  };
}

#endif

