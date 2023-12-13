#include "Minus.h"

namespace Rodin::Variational
{
  FormLanguage::List<LocalBilinearFormIntegratorBase>
  operator-(
      const LocalBilinearFormIntegratorBase& lhs,
      const LocalBilinearFormIntegratorBase& rhs)
  {
    FormLanguage::List<LocalBilinearFormIntegratorBase> res;
    res.add(lhs);
    res.add(UnaryMinus(rhs));
    return res;
  }

  FormLanguage::List<LocalBilinearFormIntegratorBase>
  operator-(
      const LocalBilinearFormIntegratorBase& lhs,
      const FormLanguage::List<LocalBilinearFormIntegratorBase>& rhs)
  {
    FormLanguage::List<LocalBilinearFormIntegratorBase> res;
    res.add(lhs);
    res.add(UnaryMinus(rhs));
    return res;
  }

  FormLanguage::List<LocalBilinearFormIntegratorBase>
  operator-(
      const FormLanguage::List<LocalBilinearFormIntegratorBase>& lhs,
      const LocalBilinearFormIntegratorBase& rhs)
  {
    FormLanguage::List<LocalBilinearFormIntegratorBase> res;
    res.add(lhs);
    res.add(UnaryMinus(rhs));
    return res;
  }

  FormLanguage::List<LocalBilinearFormIntegratorBase>
  operator-(
      const FormLanguage::List<LocalBilinearFormIntegratorBase>& lhs,
      const FormLanguage::List<LocalBilinearFormIntegratorBase>& rhs)
  {
    FormLanguage::List<LocalBilinearFormIntegratorBase> res;
    res.add(lhs);
    res.add(UnaryMinus(rhs));
    return res;
  }
}
