/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "RangeShape.h"

#include "BilinearFormIntegrator.h"

#include "Sum.h"

namespace Rodin::Variational
{
  FormLanguage::List<BilinearFormIntegratorBase>
  operator+(
      const BilinearFormIntegratorBase& lhs,
      const BilinearFormIntegratorBase& rhs)
  {
    FormLanguage::List<BilinearFormIntegratorBase> res;
    res.add(lhs);
    res.add(rhs);
    return res;
  }

  FormLanguage::List<BilinearFormIntegratorBase>
  operator+(
      const BilinearFormIntegratorBase& lhs,
      const FormLanguage::List<BilinearFormIntegratorBase>& rhs)
  {
    FormLanguage::List<BilinearFormIntegratorBase> res(rhs);
    res.add(lhs);
    return res;
  }

  FormLanguage::List<BilinearFormIntegratorBase>
  operator+(
      const FormLanguage::List<BilinearFormIntegratorBase>& lhs,
      const BilinearFormIntegratorBase& rhs)
  {
    FormLanguage::List<BilinearFormIntegratorBase> res(lhs);
    res.add(rhs);
    return res;
  }

  FormLanguage::List<BilinearFormIntegratorBase>
  operator+(
      const FormLanguage::List<BilinearFormIntegratorBase>& lhs,
      const FormLanguage::List<BilinearFormIntegratorBase>& rhs)
  {
    FormLanguage::List<BilinearFormIntegratorBase> res(lhs);
    for (const auto& p : rhs)
      res.add(p);
    return res;
  }
}
