#ifndef RODIN_VARIATIONAL_FORMLANGUAGE_BASE_H
#define RODIN_VARIATIONAL_FORMLANGUAGE_BASE_H

namespace Rodin::Variational::FormLanguage
{
   class Base
   {
      public:
         virtual ~Base() = default;

         /**
          * @internal
          * @brief Copies the object and returns a non-owning pointer to the
          * copied object.
          * @returns Non-owning pointer to the copied object.
          */
         virtual Base* copy() const noexcept = 0;
   };
}

#endif
