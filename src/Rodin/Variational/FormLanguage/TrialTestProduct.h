#ifndef RODIN_VARIATIONAL_FORMLANGUAGE_TRIALTESTPRODUCT_H
#define RODIN_VARIATIONAL_FORMLANGUAGE_TRIALTESTPRODUCT_H

#include <memory>

#include "ForwardDecls.h"

namespace Rodin::Variational::FormLanguage
{
   class TrialTestProduct
   {
      public:
         TrialTestProduct(const TrialFunctionBase& trial, const TestFunctionBase& test);

         const TrialFunctionBase& getTrialFunction() const;
         const TestFunctionBase& getTestFunction() const;

      private:
         std::unique_ptr<TrialFunctionBase> m_trial;
         std::unique_ptr<TestFunctionBase> m_test;
   };

   TrialTestProduct operator*(const TrialFunctionBase&, const TestFunctionBase&);
}

#endif
