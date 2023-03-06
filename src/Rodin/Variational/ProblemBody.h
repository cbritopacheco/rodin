/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_PROBLEMBODY_H
#define RODIN_VARIATIONAL_PROBLEMBODY_H

#include <vector>
#include <memory>
#include <optional>

#include "Rodin/FormLanguage/Base.h"
#include "Rodin/FormLanguage/List.h"

#include "ForwardDecls.h"

#include "BilinearFormIntegrator.h"
#include "LinearFormIntegrator.h"
#include "DirichletBC.h"
#include "UnaryMinus.h"

namespace Rodin::Variational
{
  /**
   * @brief Represents the body of a variational problem.
   */
  class ProblemBody final : public FormLanguage::Base
  {
    public:
      using Parent = FormLanguage::Base;

      ProblemBody() = default;

      ProblemBody(const BilinearFormIntegratorBase& bfi)
      {
        m_bfis.add(bfi);
      }

      ProblemBody(const FormLanguage::List<BilinearFormIntegratorBase>& bfis)
        : m_bfis(bfis)
      {}

      ProblemBody(const ProblemBody& other)
        : Parent(other),
          m_bfis(other.m_bfis),
          m_lfis(other.m_lfis),
          m_essBdr(other.m_essBdr)
      {}

      ProblemBody(ProblemBody&& other)
        : Parent(std::move(other)),
          m_bfis(std::move(other.m_bfis)),
          m_lfis(std::move(other.m_lfis)),
          m_essBdr(std::move(other.m_essBdr))
      {}

      inline
      ProblemBody& operator=(ProblemBody&& other)
      {
        m_bfis = std::move(other.m_bfis);
        m_lfis = std::move(other.m_lfis);
        m_essBdr = std::move(other.m_essBdr);
        return *this;
      }

      const FormLanguage::List<DirichletBCBase>& getDBCs() const
      {
        return m_essBdr;
      }

      const FormLanguage::List<BilinearFormIntegratorBase>& getBFIs() const
      {
        return m_bfis;
      }

      const FormLanguage::List<LinearFormIntegratorBase>& getLFIs() const
      {
        return m_lfis;
      }

      inline ProblemBody* copy() const noexcept override
      {
        return new ProblemBody(*this);
      }

      friend
      ProblemBody operator+(
          const ProblemBody& pb, const LinearFormIntegratorBase& lfi);

      friend
      ProblemBody operator+(
          const ProblemBody& pb, const FormLanguage::List<LinearFormIntegratorBase>& lfis);

      friend
      ProblemBody operator-(
          const ProblemBody& pb, const LinearFormIntegratorBase& lfi);

      friend
      ProblemBody operator-(
          const ProblemBody& pb, const FormLanguage::List<LinearFormIntegratorBase>& lfis);

      friend
      ProblemBody operator+(
          const ProblemBody& pb, const DirichletBCBase& dbc);

      friend
      ProblemBody operator+(
          const ProblemBody& pb, const FormLanguage::List<DirichletBCBase>& dbcs);

    private:
      FormLanguage::List<BilinearFormIntegratorBase> m_bfis;
      FormLanguage::List<LinearFormIntegratorBase> m_lfis;
      FormLanguage::List<DirichletBCBase> m_essBdr;
  };

  ProblemBody operator+(
      const ProblemBody& pb, const LinearFormIntegratorBase& lfi);

  ProblemBody operator+(
      const ProblemBody& pb, const FormLanguage::List<LinearFormIntegratorBase>& lfis);

  ProblemBody operator-(
      const ProblemBody& pb, const LinearFormIntegratorBase& lfi);

  ProblemBody operator-(
      const ProblemBody& pb, const FormLanguage::List<LinearFormIntegratorBase>& lfis);

  ProblemBody operator+(
      const ProblemBody& pb, const DirichletBCBase& dbc);

  ProblemBody operator+(
      const ProblemBody& pb, const FormLanguage::List<DirichletBCBase>& dbcs);
}

#endif
