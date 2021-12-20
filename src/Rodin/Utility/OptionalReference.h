#ifndef RODIN_UTILITY_OPTIONALREFERENCE_H
#define RODIN_UTILITY_OPTIONALREFERENCE_H

#include <optional>
#include <functional>

namespace Rodin::Utility
{
   template <typename T>
   class OptionalReference : public std::optional<std::reference_wrapper<T>>
   {
      public:
         using Parent = std::optional<std::reference_wrapper<T>>;
         using Parent::Parent;

         T* operator->()
         {
            return &(this->Parent::operator*().get());
         }

         T& operator*()
         {
            return this->Parent::operator*().get();
         }
   };
}

#endif
