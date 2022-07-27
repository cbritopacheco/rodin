#include "RangeShape.h"

namespace Rodin::Variational
{
   std::ostream& operator<<(std::ostream& os, const RangeShape& obj)
   {
      os << "{" << obj.height() << ", " << obj.width() << "}";
      return os;
   }

   std::ostream& operator<<(std::ostream& os, const RangeType& obj)
   {
      switch (obj)
      {
         case RangeType::Scalar:
         {
            os << "Scalar";
            break;
         }
         case RangeType::Vector:
         {
            os << "Vector";
            break;
         }
         case RangeType::Matrix:
         {
            os << "Matrix";
            break;
         }
      }
      return os;
   }
}
