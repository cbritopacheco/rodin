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
  ProblemBody operator+(const ProblemBody& pb, const LinearFormIntegratorBase& lfi)
  {
    ProblemBody res(pb);
    res.m_lfis.add(lfi);
    return res;
  }

  ProblemBody operator+(
      const ProblemBody& pb, const FormLanguage::List<LinearFormIntegratorBase>& lfis)
  {
    ProblemBody res(pb);
    res.m_lfis.add(lfis);
    return res;
  }

  ProblemBody operator-(
      const ProblemBody& pb, const LinearFormIntegratorBase& lfi)
  {
    ProblemBody res(pb);
    res.m_lfis.add(UnaryMinus(lfi));
    return res;
  }

  ProblemBody operator-(
      const ProblemBody& pb, const FormLanguage::List<LinearFormIntegratorBase>& lfis)
  {
    ProblemBody res(pb);
    res.m_lfis.add(UnaryMinus(lfis));
    return res;
  }

  ProblemBody operator+(
      const ProblemBody& pb, const DirichletBCBase& dbc)
  {
    ProblemBody res(pb);
    res.m_essBdr.add(dbc);
    return res;
  }

  ProblemBody operator+(
      const ProblemBody& pb, const FormLanguage::List<DirichletBCBase>& dbcs)
  {
    ProblemBody res(pb);
    res.m_essBdr.add(dbcs);
    return res;
  }

  ProblemBody operator+(
      const ProblemBody& pb, const PeriodicBCBase& pbc)
  {
    ProblemBody res(pb);
    res.m_periodicBdr.add(pbc);
    return res;
  }

  ProblemBody operator+(
      const ProblemBody& pb, const FormLanguage::List<PeriodicBCBase>& pbcs)
  {
    ProblemBody res(pb);
    res.m_periodicBdr.add(pbcs);
    return res;
  }
}
