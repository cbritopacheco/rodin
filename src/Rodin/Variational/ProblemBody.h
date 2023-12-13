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

#include "UnaryMinus.h"
#include "PeriodicBC.h"
#include "DirichletBC.h"
#include "LinearFormIntegrator.h"
#include "BilinearFormIntegrator.h"
#include "Potential.h"

namespace Rodin::Variational
{
  /**
   * @brief Represents the body of a variational problem.
   */
  class ProblemBodyBase : public FormLanguage::Base
  {
    public:
      using Parent = FormLanguage::Base;

      ProblemBodyBase() = default;

      ProblemBodyBase(const ProblemBodyBase& other)
        : Parent(other),
          m_bfis(other.m_bfis),
          m_lfis(other.m_lfis),
          m_essBdr(other.m_essBdr)
      {}

      ProblemBodyBase(ProblemBodyBase&& other)
        : Parent(std::move(other)),
          m_bfis(std::move(other.m_bfis)),
          m_lfis(std::move(other.m_lfis)),
          m_essBdr(std::move(other.m_essBdr))
      {}

      inline
      ProblemBodyBase& operator=(ProblemBodyBase&& other)
      {
        m_bfis = std::move(other.m_bfis);
        m_lfis = std::move(other.m_lfis);
        m_essBdr = std::move(other.m_essBdr);
        return *this;
      }

      PeriodicBoundary& getPBCs()
      {
        return m_periodicBdr;
      }

      EssentialBoundary& getDBCs()
      {
        return m_essBdr;
      }

      FormLanguage::List<LocalBilinearFormIntegratorBase>& getBFIs()
      {
        return m_bfis;
      }

      FormLanguage::List<LinearFormIntegratorBase>& getLFIs()
      {
        return m_lfis;
      }

      const PeriodicBoundary& getPBCs() const
      {
        return m_periodicBdr;
      }

      const EssentialBoundary& getDBCs() const
      {
        return m_essBdr;
      }

      const FormLanguage::List<LocalBilinearFormIntegratorBase>& getBFIs() const
      {
        return m_bfis;
      }

      const FormLanguage::List<LinearFormIntegratorBase>& getLFIs() const
      {
        return m_lfis;
      }

      virtual ProblemBodyBase* copy() const noexcept override
      {
        return new ProblemBodyBase(*this);
      }

    private:
      FormLanguage::List<LocalBilinearFormIntegratorBase> m_bfis;
      FormLanguage::List<LinearFormIntegratorBase> m_lfis;
      EssentialBoundary m_essBdr;
      PeriodicBoundary  m_periodicBdr;
  };

  template <>
  class ProblemBody<void, void> : public ProblemBodyBase
  {
    public:
      using Parent = ProblemBodyBase;

      ProblemBody() = default;

      ProblemBody(const ProblemBody& other)
        : Parent(other)
      {}

      ProblemBody(ProblemBody&& other)
        : Parent(std::move(other))
      {}

      virtual ProblemBody* copy() const noexcept override
      {
        return new ProblemBody(*this);
      }
  };

  template <class OperatorType>
  class ProblemBody<OperatorType, void> : public ProblemBodyBase
  {
    public:
      using Parent = ProblemBodyBase;

      ProblemBody() = default;

      ProblemBody(const ProblemBody& other)
        : Parent(other),
          m_bfs(other.m_bfs)
      {}

      ProblemBody(ProblemBody&& other)
        : Parent(std::move(other)),
          m_bfs(std::move(other.m_bfs))
      {}

      FormLanguage::List<BilinearFormBase<OperatorType>>& getBFs()
      {
        return m_bfs;
      }

      const FormLanguage::List<BilinearFormBase<OperatorType>>& getBFs() const
      {
        return m_bfs;
      }

      virtual ProblemBody* copy() const noexcept override
      {
        return new ProblemBody(*this);
      }

    private:
      FormLanguage::List<BilinearFormBase<OperatorType>> m_bfs;
  };

  template <class VectorType>
  class ProblemBody<void, VectorType> : public ProblemBodyBase
  {
    public:
      using Parent = ProblemBodyBase;

      ProblemBody() = default;

      ProblemBody(const ProblemBody& other)
        : Parent(other),
          m_lfs(other.m_lfs)
      {}

      ProblemBody(ProblemBody&& other)
        : Parent(std::move(other)),
          m_lfs(std::move(other.m_lfs))
      {}

      inline
      FormLanguage::List<LinearFormBase<VectorType>>& getLFs()
      {
        return m_lfs;
      }

      inline
      const FormLanguage::List<LinearFormBase<VectorType>>& getLFs() const
      {
        return m_lfs;
      }

      virtual ProblemBody* copy() const noexcept override
      {
        return new ProblemBody(*this);
      }

    private:
      FormLanguage::List<LinearFormBase<VectorType>> m_lfs;
  };

  template <class OperatorType, class VectorType>
  class ProblemBody
    : public ProblemBodyBase
  {
    public:
      using Parent = ProblemBodyBase;

      ProblemBody(const LocalBilinearFormIntegratorBase& bfi)
      {
        getBFIs().add(bfi);
      }

      ProblemBody(const ProblemBody<OperatorType, void>& pbo)
        : Parent(pbo)
      {
        m_bfs.add(pbo.getBFs());
      }

      ProblemBody(const ProblemBody<void, VectorType>& pbv)
        : Parent(pbv)
      {
        m_lfs.add(pbv.getLFs());
      }

      ProblemBody(const ProblemBody<void, void>& parent)
        : Parent(parent)
      {}

      ProblemBody(const ProblemBody& other)
        : Parent(other),
          m_lfs(other.m_lfs),
          m_bfs(other.m_bfs)
      {}

      ProblemBody(ProblemBody&& other)
        : Parent(std::move(other)),
          m_lfs(std::move(other.m_lfs)),
          m_bfs(std::move(other.m_bfs))
      {}

      inline
      FormLanguage::List<LinearFormBase<VectorType>>& getLFs()
      {
        return m_lfs;
      }

      inline
      FormLanguage::List<BilinearFormBase<OperatorType>>& getBFs()
      {
        return m_bfs;
      }

      inline
      const FormLanguage::List<LinearFormBase<VectorType>>& getLFs() const
      {
        return m_lfs;
      }

      inline
      const FormLanguage::List<BilinearFormBase<OperatorType>>& getBFs() const
      {
        return m_bfs;
      }

      virtual ProblemBody* copy() const noexcept override
      {
        return new ProblemBody(*this);
      }

    private:
      FormLanguage::List<LinearFormBase<VectorType>> m_lfs;
      FormLanguage::List<BilinearFormBase<OperatorType>> m_bfs;
  };

  ProblemBody(const LocalBilinearFormIntegratorBase&)
    -> ProblemBody<void, void>;

  inline
  ProblemBody<void, void> operator+(
      const LocalBilinearFormIntegratorBase& bfi, const LinearFormIntegratorBase& lfi)
  {
    ProblemBody<void, void> res;
    res.getBFIs().add(bfi);
    res.getLFIs().add(lfi);
    return res;
  }

  inline
  ProblemBody<void, void> operator+(
      const LinearFormIntegratorBase& lfi, const LocalBilinearFormIntegratorBase& bfi)
  {
    ProblemBody<void, void> res;
    res.getBFIs().add(bfi);
    res.getLFIs().add(lfi);
    return res;
  }

  inline
  ProblemBody<void, void> operator-(
      const LocalBilinearFormIntegratorBase& bfi, const LinearFormIntegratorBase& lfi)
  {
    ProblemBody<void, void> res;
    res.getBFIs().add(bfi);
    res.getLFIs().add(UnaryMinus(lfi));
    return res;
  }

  inline
  ProblemBody<void, void> operator-(
      const LinearFormIntegratorBase& lfi, const LocalBilinearFormIntegratorBase& bfi)
  {
    ProblemBody<void, void> res;
    res.getBFIs().add(UnaryMinus(bfi));
    res.getLFIs().add(lfi);
    return res;
  }

  inline
  ProblemBody<void, void> operator+(
      const FormLanguage::List<LocalBilinearFormIntegratorBase>& bfis,
      const LinearFormIntegratorBase& lfi)
  {
    ProblemBody<void, void> res;
    res.getBFIs().add(bfis);
    res.getLFIs().add(lfi);
    return res;
  }

  inline
  ProblemBody<void, void> operator-(
      const FormLanguage::List<LocalBilinearFormIntegratorBase>& bfis,
      const LinearFormIntegratorBase& lfi)
  {
    ProblemBody<void, void> res;
    res.getBFIs().add(bfis);
    res.getLFIs().add(UnaryMinus(lfi));
    return res;
  }

  inline
  ProblemBody<void, void> operator+(
      const LocalBilinearFormIntegratorBase& bfi, const DirichletBCBase& dbc)
  {
    ProblemBody<void, void> res;
    res.getBFIs().add(bfi);
    res.getDBCs().add(dbc);
    return res;
  }

  inline
  ProblemBody<void, void> operator+(
      const LocalBilinearFormIntegratorBase& bfi, const FormLanguage::List<DirichletBCBase>& dbcs)
  {
    ProblemBody<void, void> res;
    res.getBFIs().add(bfi);
    res.getDBCs().add(dbcs);
    return res;
  }

  inline
  ProblemBody<void, void> operator+(
      const LocalBilinearFormIntegratorBase& bfi, const PeriodicBCBase& pbc)
  {
    ProblemBody<void, void> res;
    res.getBFIs().add(bfi);
    res.getPBCs().add(pbc);
    return res;
  }

  inline
  ProblemBody<void, void> operator+(
      const LocalBilinearFormIntegratorBase& bfi, const FormLanguage::List<PeriodicBCBase>& pbcs)
  {
    ProblemBody<void, void> res;
    res.getBFIs().add(bfi);
    res.getPBCs().add(pbcs);
    return res;
  }


  inline
  ProblemBody<void, void> operator+(
      const FormLanguage::List<LocalBilinearFormIntegratorBase>& bfis, const DirichletBCBase& dbc)
  {
    ProblemBody<void, void> res;
    res.getBFIs().add(bfis);
    res.getDBCs().add(dbc);
    return res;
  }

  inline
  ProblemBody<void, void> operator+(
      const FormLanguage::List<LocalBilinearFormIntegratorBase>& bfis, const FormLanguage::List<DirichletBCBase>& dbcs)
  {
    ProblemBody<void, void> res;
    res.getBFIs().add(bfis);
    res.getDBCs().add(dbcs);
    return res;
  }

  inline
  ProblemBody<void, void> operator+(
      const FormLanguage::List<LocalBilinearFormIntegratorBase>& bfis, const PeriodicBCBase& pbc)
  {
    ProblemBody<void, void> res;
    res.getBFIs().add(bfis);
    res.getPBCs().add(pbc);
    return res;
  }

  inline
  ProblemBody<void, void> operator+(
      const FormLanguage::List<LocalBilinearFormIntegratorBase>& bfis, const FormLanguage::List<PeriodicBCBase>& pbcs)
  {
    ProblemBody<void, void> res;
    res.getBFIs().add(bfis);
    res.getPBCs().add(pbcs);
    return res;
  }

  template <class OperatorType, class VectorType>
  ProblemBody<OperatorType, VectorType> operator+(
      const ProblemBody<OperatorType, VectorType>& pb,
      const LinearFormIntegratorBase& lfi)
  {
    ProblemBody<OperatorType, VectorType> res(pb);
    res.getLFIs().add(lfi);
    return res;
  }

  template <class OperatorType, class VectorType>
  ProblemBody<OperatorType, VectorType> operator+(
      const ProblemBody<OperatorType, VectorType>& pb,
      const FormLanguage::List<LinearFormIntegratorBase>& lfis)
  {
    ProblemBody<OperatorType, VectorType> res(pb);
    res.getLFIs().add(lfis);
    return res;
  }

  template <class OperatorType, class VectorType>
  ProblemBody<OperatorType, VectorType> operator-(
      const ProblemBody<OperatorType, VectorType>& pb,
      const LinearFormIntegratorBase& lfi)
  {
    ProblemBody<OperatorType, VectorType> res(pb);
    res.getLFIs().add(UnaryMinus(lfi));
    return res;
  }

  template <class OperatorType, class VectorType>
  ProblemBody<OperatorType, VectorType> operator-(
      const ProblemBody<OperatorType, VectorType>& pb,
      const FormLanguage::List<LinearFormIntegratorBase>& lfis)
  {
    ProblemBody<OperatorType, VectorType> res(pb);
    res.getLFIs().add(UnaryMinus(lfis));
    return res;
  }

  template <class OperatorType, class VectorType>
  ProblemBody<OperatorType, VectorType> operator+(
      const ProblemBody<OperatorType, VectorType>& pb, const DirichletBCBase& dbc)
  {
    ProblemBody<OperatorType, VectorType> res(pb);
    res.getDBCs().add(dbc);
    return res;
  }

  template <class OperatorType, class VectorType>
  ProblemBody<OperatorType, VectorType> operator+(
      const ProblemBody<OperatorType, VectorType>& pb, const FormLanguage::List<DirichletBCBase>& dbcs)
  {
    ProblemBody<OperatorType, VectorType> res(pb);
    res.getEssentialBoundary().add(dbcs);
    return res;
  }

  template <class OperatorType, class VectorType>
  ProblemBody<OperatorType, VectorType> operator+(
      const ProblemBody<OperatorType, VectorType>& pb, const PeriodicBCBase& pbc)
  {
    ProblemBody<OperatorType, VectorType> res(pb);
    res.getPBCs().add(pbc);
    return res;
  }

  template <class OperatorType, class VectorType>
  ProblemBody<OperatorType, VectorType> operator+(
      const ProblemBody<OperatorType, VectorType>& pb,
      const BilinearFormBase<OperatorType>& bf)
  {
    ProblemBody<OperatorType, VectorType> res(pb);
    res.getBFs().add(bf);
    return res;
  }

  template <class OperatorType, class VectorType>
  ProblemBody<OperatorType, VectorType> operator+(
      const ProblemBody<OperatorType, VectorType>& pb, const FormLanguage::List<PeriodicBCBase>& pbcs)
  {
    ProblemBody<OperatorType, VectorType> res(pb);
    res.getPBCs().add(pbcs);
    return res;
  }

  template <class OperatorType>
  ProblemBody<OperatorType, void> operator+(
      const LocalBilinearFormIntegratorBase& bfi, const BilinearFormBase<OperatorType>& bf)
  {
    ProblemBody<OperatorType, void> res;
    res.getBFIs().add(bfi);
    res.getBFs().add(bf);
    return res;
  }

  template <class OperatorType>
  ProblemBody<OperatorType, void> operator-(
      const BilinearFormBase<OperatorType>& bf, const LinearFormIntegratorBase& lfi)
  {
    ProblemBody<OperatorType, void> res;
    res.getBFs().add(bf);
    res.getLFIs().add(UnaryMinus(lfi));
    return res;
  }

  template <class OperatorType>
  ProblemBody<OperatorType, void> operator-(
      const FormLanguage::List<BilinearFormBase<OperatorType>>& bfs,
      const LinearFormIntegratorBase& lfi)
  {
    ProblemBody<OperatorType, void> res;
    res.getBFs().add(bfs);
    res.getLFIs().add(UnaryMinus(lfi));
    return res;
  }

  template <class OperatorType>
  ProblemBody<OperatorType, void> operator+(
      const LocalBilinearFormIntegratorBase& bfi,
      const FormLanguage::List<BilinearFormBase<OperatorType>>& bfs)
  {
    ProblemBody<OperatorType, void> res;
    res.getBFIs().add(bfi);
    res.getBFs().add(bfs);
    return res;
  }
}

#endif
