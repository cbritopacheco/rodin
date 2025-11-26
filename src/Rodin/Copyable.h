/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_COPYABLE_H
#define RODIN_COPYABLE_H

/**
 * @file
 * @brief Defines the Copyable interface for polymorphic copying.
 */

namespace Rodin
{
  /**
   * @brief Abstract base class for objects that can be copied.
   *
   * This class defines an interface for polymorphic copying of objects.
   * Classes that inherit from this interface must implement the copy() 
   * method to provide deep copying functionality.
   */
  class Copyable
  {
    public:
      /// Virtual destructor for proper cleanup of derived classes
      virtual ~Copyable() = default;

      /**
       * @brief Creates a polymorphic copy of this object.
       * @return Pointer to a new instance that is a copy of this object.
       *         The caller is responsible for memory management.
       * @note This method enables polymorphic copy behavior. The copy 
       *       semantics (shallow vs deep) depend on the concrete implementation.
       *       This method should not throw exceptions.
       */
      virtual Copyable* copy() const noexcept = 0;
  };
}

#endif

