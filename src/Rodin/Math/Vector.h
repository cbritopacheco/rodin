/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MATH_VECTOR_H
#define RODIN_MATH_VECTOR_H

#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>

namespace Rodin::Math
{
  class Vector : public Eigen::VectorX<double>
  {
    public:
      using Parent = Eigen::VectorX<double>;

      Vector(std::initializer_list<Parent::Scalar> l)
        : Parent(l.size())
      {
        std::copy(l.begin(), l.end(), Parent::begin());
      }

      template <class ... Args>
      Vector(Args&&... args)
        : Parent(std::forward<Args>(args)...)
      {}

      template <class ... Args>
      Vector& operator=(Args&&... args)
      {
         this->Parent::operator=(std::forward<Args>(args)...);
         return *this;
      }
  };
}

#endif
