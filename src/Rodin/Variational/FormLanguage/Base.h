#ifndef RODIN_VARIATIONAL_FORMLANGUAGE_BASE_H
#define RODIN_VARIATIONAL_FORMLANGUAGE_BASE_H

#include <cassert>

namespace Rodin::Variational::FormLanguage
{
   /**
    * @brief Base class for all classes which are part of Variational::FormLanguage.
    */
   class Base
   {
      public:
         /**
          * @brief Virtual destructor.
          */
         virtual ~Base() = default;

         /**
          * @internal
          * @brief Copies the object and returns a non-owning pointer to the
          * copied object.
          * @returns Non-owning pointer to the copied object.
          */
         virtual Base* copy() const noexcept = 0;
   };

   template <class InternalValue>
   class Buildable : public Base
   {
      public:
         virtual void build() = 0;

         virtual InternalValue& get() = 0;

         virtual InternalValue* release()
         {
            assert(false); // Unimplemented method
            return nullptr;
         }
   };
}

#endif
