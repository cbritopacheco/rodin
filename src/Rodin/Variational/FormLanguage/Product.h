#ifndef RODIN_VARIATIONAL_FORMLANGUAGE_PRODUCT_H
#define RODIN_VARIATIONAL_FORMLANGUAGE_PRODUCT_H

#include <memory>
#include <type_traits>

#include "ForwardDecls.h"

namespace Rodin::Variational::FormLanguage
{
   template <class Lhs, class Rhs>
   class Product
   {
      static_assert(std::is_base_of_v<Base, Lhs>,
            "Lhs must be derived from FormLanguage::Base");
      static_assert(std::is_base_of_v<Base, Rhs>,
            "Rhs must be derived from FormLanguage::Base");

      public:
         Product(const Lhs& lhs, const Rhs& rhs)
            : m_lhs(lhs.copy()), m_rhs(rhs.copy())
         {}

         Product(const Product& other)
            : m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
         {}

         const Lhs& getLHS() const
         {
            return *m_lhs;
         }

         const Rhs& getRHS() const
         {
            return *m_rhs;
         }

      private:
         std::unique_ptr<Lhs> m_lhs;
         std::unique_ptr<Rhs> m_rhs;
   };

   // ---- Rules -------------------------------------------------------------
   Product<TrialFunctionBase, TestFunctionBase>
   operator*(const TrialFunctionBase& lhs, const TestFunctionBase& rhs);

   // ---- Rules -------------------------------------------------------------
   Product<ScalarCoefficientBase, TestFunctionBase>
   operator*(const ScalarCoefficientBase& lhs, const TestFunctionBase& rhs);
}

#endif
