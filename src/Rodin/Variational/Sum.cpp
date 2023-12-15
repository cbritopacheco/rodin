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
  FormLanguage::List<LocalBilinearFormIntegratorBase>
  operator+(
      const LocalBilinearFormIntegratorBase& lhs,
      const LocalBilinearFormIntegratorBase& rhs)
  {
    FormLanguage::List<LocalBilinearFormIntegratorBase> res;
    res.add(lhs);
    res.add(rhs);
    return res;
  }

  FormLanguage::List<LocalBilinearFormIntegratorBase>
  operator+(
      const LocalBilinearFormIntegratorBase& lhs,
      const FormLanguage::List<LocalBilinearFormIntegratorBase>& rhs)
  {
    FormLanguage::List<LocalBilinearFormIntegratorBase> res(rhs);
    res.add(lhs);
    return res;
  }

  FormLanguage::List<LocalBilinearFormIntegratorBase>
  operator+(
      const FormLanguage::List<LocalBilinearFormIntegratorBase>& lhs,
      const LocalBilinearFormIntegratorBase& rhs)
  {
    FormLanguage::List<LocalBilinearFormIntegratorBase> res(lhs);
    res.add(rhs);
    return res;
  }

  FormLanguage::List<LocalBilinearFormIntegratorBase>
  operator+(
      const FormLanguage::List<LocalBilinearFormIntegratorBase>& lhs,
      const FormLanguage::List<LocalBilinearFormIntegratorBase>& rhs)
  {
    FormLanguage::List<LocalBilinearFormIntegratorBase> res(lhs);
    for (const auto& p : rhs)
      res.add(p);
    return res;
  }
}
