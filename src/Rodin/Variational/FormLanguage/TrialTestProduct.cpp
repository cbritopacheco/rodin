#include "../TrialFunction.h"
#include "../TestFunction.h"

#include "TrialTestProduct.h"

namespace Rodin::Variational::FormLanguage
{
   TrialTestProduct::TrialTestProduct(
         const TrialFunctionBase& trial, const TestFunctionBase& test)
      : m_trial(trial.copy()), m_test(test.copy())
   {}

   const TrialFunctionBase& TrialTestProduct::getTrialFunction() const
   {
      return *m_trial;
   }

   const TestFunctionBase& TrialTestProduct::getTestFunction() const
   {
      return *m_test;
   }
}
