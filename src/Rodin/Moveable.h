/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MOVEABLE_H
#define RODIN_MOVEABLE_H

/**
 * @file
 * @brief Defines the Moveable interface for polymorphic move operations.
 */

namespace Rodin
{
  /**
   * @brief Abstract base class for objects that can be moved.
   *
   * This class defines an interface for polymorphic moving of objects.
   * Classes that inherit from this interface must implement the move() 
   * method to provide transfer of ownership functionality.
   * 
   * @note Unlike Copyable which creates new instances, Moveable transfers
   *       ownership of the existing object.
   */
  class Moveable
  {
    public:
      virtual ~Moveable() = default;

      /**
       * @brief Moves this object and transfers ownership.
       * @return Pointer to this object after the move operation.
       *         The caller assumes ownership responsibility.
       * @note This method should not throw exceptions.
       */
      virtual Moveable* move() noexcept = 0;
  };
}

#endif


