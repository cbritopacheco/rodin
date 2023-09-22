#include "Minus.h"

namespace Rodin::Variational
{
  FormLanguage::List<BilinearFormIntegratorBase>
  operator-(
      const BilinearFormIntegratorBase& lhs,
      const BilinearFormIntegratorBase& rhs)
  {
    FormLanguage::List<BilinearFormIntegratorBase> res;
    res.add(lhs);
    res.add(UnaryMinus(rhs));
    return res;
  }

  FormLanguage::List<BilinearFormIntegratorBase>
  operator-(
      const BilinearFormIntegratorBase& lhs,
      const FormLanguage::List<BilinearFormIntegratorBase>& rhs)
  {
    FormLanguage::List<BilinearFormIntegratorBase> res;
    res.add(lhs);
    res.add(UnaryMinus(rhs));
    return res;
  }

  FormLanguage::List<BilinearFormIntegratorBase>
  operator-(
      const FormLanguage::List<BilinearFormIntegratorBase>& lhs,
      const BilinearFormIntegratorBase& rhs)
  {
    FormLanguage::List<BilinearFormIntegratorBase> res;
    res.add(lhs);
    res.add(UnaryMinus(rhs));
    return res;
  }

  FormLanguage::List<BilinearFormIntegratorBase>
  operator-(
      const FormLanguage::List<BilinearFormIntegratorBase>& lhs,
      const FormLanguage::List<BilinearFormIntegratorBase>& rhs)
  {
    FormLanguage::List<BilinearFormIntegratorBase> res;
    res.add(lhs);
    res.add(UnaryMinus(rhs));
    return res;
  }
}
