#ifndef RODIN_COPYABLE_H
#define RODIN_COPYABLE_H

namespace Rodin
{
  class Copyable
  {
    public:
      virtual Copyable* copy() const noexcept = 0;
  };
}

#endif

