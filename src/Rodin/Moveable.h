#ifndef RODIN_MOVEABLE_H
#define RODIN_MOVEABLE_H

namespace Rodin
{
  class Moveable
  {
    public:
      virtual Moveable* move() noexcept = 0;
  };
}

#endif


