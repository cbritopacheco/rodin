#ifndef RODIN_MOVEABLE_H
#define RODIN_MOVEABLE_H

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


