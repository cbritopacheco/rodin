#ifndef RODIN_VARIATIONAL_FORMLANGUAGE_BASE_H
#define RODIN_VARIATIONAL_FORMLANGUAGE_BASE_H

namespace Rodin::Variational::FormLanguage
{
   class Base
   {
      public:
         virtual ~Base() = default;
         virtual Base* copy() const noexcept = 0;
   };
}

#endif
