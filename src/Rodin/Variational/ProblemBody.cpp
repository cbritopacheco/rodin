/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Component.h"
#include "UnaryMinus.h"

#include "ProblemBody.h"

namespace Rodin::Variational
{
  ProblemBody operator+(
      const ProblemBody& pb, const LinearFormIntegratorBase& lfi)
  {
    ProblemBody res(pb);
    res.getLFIs().add(lfi);
    return res;
  }

  ProblemBody operator-(const ProblemBody& pb, const LinearFormIntegratorBase& lfi)
  {
    ProblemBody res(pb);
    res.getLFIs().add(UnaryMinus(lfi));
    return res;
  }

  ProblemBody operator+(
    const ProblemBody& pb, const FormLanguage::List<LinearFormIntegratorBase>& lfis)
  {
    ProblemBody res(pb);
    res.getLFIs().add(lfis);
    return res;
  }

  ProblemBody operator-(
    const ProblemBody& pb, const FormLanguage::List<LinearFormIntegratorBase>& lfi)
  {
    ProblemBody res(pb);
    res.getLFIs().add(UnaryMinus(lfi));
    return res;
  }
}
