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
      using ScalarType = Real;

      using LinearFormIntegratorBaseType = LinearFormIntegratorBase<ScalarType>;

      using LocalBilinearFormIntegratorBaseType = LocalBilinearFormIntegratorBase<ScalarType>;

      using GlobalBilinearFormIntegratorBaseType = GlobalBilinearFormIntegratorBase<ScalarType>;

      using LinearFormIntegratorBaseListType = FormLanguage::List<LinearFormIntegratorBaseType>;

      using LocalBilinearFormIntegratorBaseListType = FormLanguage::List<LocalBilinearFormIntegratorBaseType>;

      using GlobalBilinearFormIntegratorBaseListType = FormLanguage::List<GlobalBilinearFormIntegratorBaseType>;

      using Parent = FormLanguage::Base;

      ProblemBodyBase() = default;

      ProblemBodyBase(const ProblemBodyBase& other)
        : Parent(other),
          m_lfis(other.m_lfis),
          m_lbfis(other.m_lbfis),
          m_gbfis(other.m_gbfis),
          m_essBdr(other.m_essBdr),
          m_periodicBdr(other.m_periodicBdr)
      {}

      ProblemBodyBase(ProblemBodyBase&& other)
        : Parent(std::move(other)),
          m_lfis(std::move(other.m_lfis)),
          m_lbfis(std::move(other.m_lbfis)),
          m_gbfis(std::move(other.m_gbfis)),
          m_essBdr(std::move(other.m_essBdr)),
          m_periodicBdr(std::move(other.m_periodicBdr))
      {}

      inline
      ProblemBodyBase& operator=(ProblemBodyBase&& other)
      {
        m_lfis = std::move(other.m_lfis);
        m_lbfis = std::move(other.m_lbfis);
        m_gbfis = std::move(other.m_gbfis);
        m_essBdr = std::move(other.m_essBdr);
        m_periodicBdr = std::move(other.m_periodicBdr);
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

      LocalBilinearFormIntegratorBaseListType& getLocalBFIs()
      {
        return m_lbfis;
      }

      GlobalBilinearFormIntegratorBaseListType& getGlobalBFIs()
      {
        return m_gbfis;
      }

      LinearFormIntegratorBaseListType& getLFIs()
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

      const LinearFormIntegratorBaseListType& getLFIs() const
      {
        return m_lfis;
      }

      const LocalBilinearFormIntegratorBaseListType& getLocalBFIs() const
      {
        return m_lbfis;
      }

      const GlobalBilinearFormIntegratorBaseListType& getGlobalBFIs() const
      {
        return m_gbfis;
      }

      virtual ProblemBodyBase* copy() const noexcept override
      {
        return new ProblemBodyBase(*this);
      }

    private:
      LinearFormIntegratorBaseListType m_lfis;
      LocalBilinearFormIntegratorBaseListType m_lbfis;
      GlobalBilinearFormIntegratorBaseListType m_gbfis;

      EssentialBoundary m_essBdr;
      PeriodicBoundary  m_periodicBdr;
  };

  template <>
  class ProblemBody<void, void>
    : public ProblemBodyBase
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

  template <class Operator>
  class ProblemBody<Operator, void> : public ProblemBodyBase
  {
    public:
      using OperatorType = Operator;

      using ScalarType = typename FormLanguage::Traits<OperatorType>::ScalarType;

      using BilinearFormBaseType = BilinearFormBase<OperatorType>;

      using BilinearFormBaseListType = FormLanguage::List<BilinearFormBaseType>;

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

      BilinearFormBaseListType& getBFs()
      {
        return m_bfs;
      }

      const BilinearFormBaseListType& getBFs() const
      {
        return m_bfs;
      }

      virtual ProblemBody* copy() const noexcept override
      {
        return new ProblemBody(*this);
      }

    private:
      BilinearFormBaseListType m_bfs;
  };

  template <class Vector>
  class ProblemBody<void, Vector> : public ProblemBodyBase
  {
    public:
      using VectorType = Vector;

      using LinearFormBaseType = LinearFormBase<VectorType>;

      using LinearFormBaseListType = FormLanguage::List<LinearFormBaseType>;

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
      LinearFormBaseListType& getLFs()
      {
        return m_lfs;
      }

      inline
      const LinearFormBaseListType& getLFs() const
      {
        return m_lfs;
      }

      virtual ProblemBody* copy() const noexcept override
      {
        return new ProblemBody(*this);
      }

    private:
      LinearFormBaseListType m_lfs;
  };

  template <class Operator, class Vector>
  class ProblemBody : public ProblemBodyBase
  {
    public:
      using VectorType = Vector;

      using OperatorType = Operator;

      using VectorScalarType = typename FormLanguage::Traits<VectorType>::ScalarType;

      using OperatorScalarType = typename FormLanguage::Traits<OperatorType>::ScalarType;

      using LinearFormBaseType = LinearFormBase<VectorType>;

      using BilinearFormBaseType = BilinearFormBase<OperatorType>;

      using LinearFormBaseListType = FormLanguage::List<LinearFormBaseType>;

      using BilinearFormBaseListType = FormLanguage::List<BilinearFormBaseType>;

      using LinearFormIntegratorBaseType = LinearFormIntegratorBase<VectorScalarType>;

      using LocalBilinearFormIntegratorBaseType = LocalBilinearFormIntegratorBase<OperatorScalarType>;

      using GlobalBilinearFormIntegratorBaseType = GlobalBilinearFormIntegratorBase<OperatorScalarType>;

      using LinearFormIntegratorBaseListType = FormLanguage::List<LinearFormIntegratorBaseType>;

      using LocalBilinearFormIntegratorBaseListType = FormLanguage::List<LocalBilinearFormIntegratorBaseType>;

      using GlobalBilinearFormIntegratorBaseListType = FormLanguage::List<GlobalBilinearFormIntegratorBaseType>;

      using Parent = ProblemBodyBase;

      ProblemBody(const LocalBilinearFormIntegratorBaseType& bfi)
      {
        this->getLocalBFIs().add(bfi);
      }

      ProblemBody(const GlobalBilinearFormIntegratorBaseType& bfi)
      {
        this->getGlobalBFIs().add(bfi);
      }

      ProblemBody(const LocalBilinearFormIntegratorBaseListType& bfis)
      {
        this->getLocalBFIs().add(bfis);
      }

      ProblemBody(const GlobalBilinearFormIntegratorBaseListType& bfis)
      {
        this->getGlobalBFIs().add(bfis);
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
      LinearFormBaseListType& getLFs()
      {
        return m_lfs;
      }

      inline
      BilinearFormBaseListType& getBFs()
      {
        return m_bfs;
      }

      inline
      const LinearFormBaseListType& getLFs() const
      {
        return m_lfs;
      }

      inline
      const BilinearFormBaseListType& getBFs() const
      {
        return m_bfs;
      }

      virtual ProblemBody* copy() const noexcept override
      {
        return new ProblemBody(*this);
      }

    private:
      LinearFormBaseListType m_lfs;
      BilinearFormBaseListType m_bfs;
  };

  ProblemBody(const LocalBilinearFormIntegratorBase<Real>&)
    -> ProblemBody<void, void>;

  template <class LHSNumber, class RHSNumber>
  inline
  ProblemBody<void, void> operator+(
      const LocalBilinearFormIntegratorBase<LHSNumber>& bfi, const LinearFormIntegratorBase<RHSNumber>& lfi)
  {
    ProblemBody<void, void> res;
    res.getLocalBFIs().add(bfi);
    res.getLFIs().add(lfi);
    return res;
  }

  template <class LHSNumber, class RHSNumber>
  inline
  ProblemBody<void, void> operator+(
      const LinearFormIntegratorBase<LHSNumber>& lfi, const LocalBilinearFormIntegratorBase<RHSNumber>& bfi)
  {
    ProblemBody<void, void> res;
    res.getLocalBFIs().add(bfi);
    res.getLFIs().add(lfi);
    return res;
  }

  template <class LHSNumber, class RHSNumber>
  inline
  ProblemBody<void, void> operator-(
      const LocalBilinearFormIntegratorBase<LHSNumber>& bfi, const LinearFormIntegratorBase<RHSNumber>& lfi)
  {
    ProblemBody<void, void> res;
    res.getLocalBFIs().add(bfi);
    res.getLFIs().add(UnaryMinus(lfi));
    return res;
  }

  template <class LHSNumber, class RHSNumber>
  inline
  ProblemBody<void, void> operator-(
      const GlobalBilinearFormIntegratorBase<LHSNumber>& bfi, const LinearFormIntegratorBase<RHSNumber>& lfi)
  {
    ProblemBody<void, void> res;
    res.getGlobalBFIs().add(bfi);
    res.getLFIs().add(UnaryMinus(lfi));
    return res;
  }

  template <class LHSNumber, class RHSNumber>
  inline
  ProblemBody<void, void> operator-(
      const GlobalBilinearFormIntegratorBase<LHSNumber>& bfi,
      const FormLanguage::List<LinearFormIntegratorBase<RHSNumber>>& lfis)
  {
    ProblemBody<void, void> res;
    res.getGlobalBFIs().add(bfi);
    res.getLFIs().add(UnaryMinus(lfis));
    return res;
  }

  template <class LHSNumber, class RHSNumber>
  inline
  ProblemBody<void, void> operator-(
      const LinearFormIntegratorBase<LHSNumber>& lfi, const LocalBilinearFormIntegratorBase<RHSNumber>& bfi)
  {
    ProblemBody<void, void> res;
    res.getLocalBFIs().add(UnaryMinus(bfi));
    res.getLFIs().add(lfi);
    return res;
  }

  template <class LHSNumber, class RHSNumber>
  inline
  ProblemBody<void, void> operator+(
      const FormLanguage::List<LocalBilinearFormIntegratorBase<LHSNumber>>& bfis,
      const LinearFormIntegratorBase<RHSNumber>& lfi)
  {
    ProblemBody<void, void> res;
    res.getLocalBFIs().add(bfis);
    res.getLFIs().add(lfi);
    return res;
  }

  template <class LHSNumber, class RHSNumber>
  inline
  ProblemBody<void, void> operator+(
      const FormLanguage::List<LocalBilinearFormIntegratorBase<LHSNumber>>& lbfis,
      const GlobalBilinearFormIntegratorBase<RHSNumber>& gbfi)
  {
    ProblemBody<void, void> res;
    res.getLocalBFIs().add(lbfis);
    res.getGlobalBFIs().add(gbfi);
    return res;
  }

  template <class LHSNumber, class RHSNumber>
  inline
  ProblemBody<void, void> operator+(
      const FormLanguage::List<LocalBilinearFormIntegratorBase<LHSNumber>>& lbfis,
      const FormLanguage::List<GlobalBilinearFormIntegratorBase<RHSNumber>>& gbfis)
  {
    ProblemBody<void, void> res;
    res.getLocalBFIs().add(lbfis);
    res.getGlobalBFIs().add(gbfis);
    return res;
  }

  template <class LHSNumber, class RHSNumber>
  inline
  ProblemBody<void, void> operator+(
      const LocalBilinearFormIntegratorBase<LHSNumber>& lbfi,
      const FormLanguage::List<GlobalBilinearFormIntegratorBase<RHSNumber>>& gbfis)
  {
    ProblemBody<void, void> res;
    res.getLocalBFIs().add(lbfi);
    res.getGlobalBFIs().add(gbfis);
    return res;
  }

  template <class LHSNumber, class RHSNumber>
  inline
  ProblemBody<void, void> operator+(
      const LocalBilinearFormIntegratorBase<LHSNumber>& lbfi,
      const GlobalBilinearFormIntegratorBase<RHSNumber>& gbfi)
  {
    ProblemBody<void, void> res;
    res.getLocalBFIs().add(lbfi);
    res.getGlobalBFIs().add(gbfi);
    return res;
  }

  template <class LHSNumber, class RHSNumber>
  inline
  ProblemBody<void, void> operator-(
      const FormLanguage::List<LocalBilinearFormIntegratorBase<LHSNumber>>& bfis,
      const LinearFormIntegratorBase<RHSNumber>& lfi)
  {
    ProblemBody<void, void> res;
    res.getLocalBFIs().add(bfis);
    res.getLFIs().add(UnaryMinus(lfi));
    return res;
  }

  template <class LHSNumber>
  inline
  ProblemBody<void, void> operator+(
      const LocalBilinearFormIntegratorBase<LHSNumber>& bfi, const DirichletBCBase& dbc)
  {
    ProblemBody<void, void> res;
    res.getLocalBFIs().add(bfi);
    res.getDBCs().add(dbc);
    return res;
  }

  template <class LHSNumber>
  inline
  ProblemBody<void, void> operator+(
      const LocalBilinearFormIntegratorBase<LHSNumber>& bfi, const FormLanguage::List<DirichletBCBase>& dbcs)
  {
    ProblemBody<void, void> res;
    res.getLocalBFIs().add(bfi);
    res.getDBCs().add(dbcs);
    return res;
  }

  template <class LHSNumber, class RHSNumber>
  inline
  ProblemBody<void, void> operator+(
      const LocalBilinearFormIntegratorBase<LHSNumber>& bfi, const PeriodicBCBase& pbc)
  {
    ProblemBody<void, void> res;
    res.getLocalBFIs().add(bfi);
    res.getPBCs().add(pbc);
    return res;
  }

  template <class LHSNumber>
  inline
  ProblemBody<void, void> operator+(
      const LocalBilinearFormIntegratorBase<LHSNumber>& bfi, const FormLanguage::List<PeriodicBCBase>& pbcs)
  {
    ProblemBody<void, void> res;
    res.getLocalBFIs().add(bfi);
    res.getPBCs().add(pbcs);
    return res;
  }

  template <class LHSNumber>
  inline
  ProblemBody<void, void> operator+(
      const FormLanguage::List<LocalBilinearFormIntegratorBase<LHSNumber>>& bfis, const DirichletBCBase& dbc)
  {
    ProblemBody<void, void> res;
    res.getLocalBFIs().add(bfis);
    res.getDBCs().add(dbc);
    return res;
  }

  template <class LHSNumber>
  inline
  ProblemBody<void, void> operator+(
      const FormLanguage::List<LocalBilinearFormIntegratorBase<LHSNumber>>& bfis, const FormLanguage::List<DirichletBCBase>& dbcs)
  {
    ProblemBody<void, void> res;
    res.getLocalBFIs().add(bfis);
    res.getDBCs().add(dbcs);
    return res;
  }

  template <class LHSNumber>
  inline
  ProblemBody<void, void> operator+(
      const FormLanguage::List<LocalBilinearFormIntegratorBase<LHSNumber>>& bfis, const PeriodicBCBase& pbc)
  {
    ProblemBody<void, void> res;
    res.getLocalBFIs().add(bfis);
    res.getPBCs().add(pbc);
    return res;
  }

  template <class LHSNumber>
  inline
  ProblemBody<void, void> operator+(
      const FormLanguage::List<LocalBilinearFormIntegratorBase<LHSNumber>>& bfis, const FormLanguage::List<PeriodicBCBase>& pbcs)
  {
    ProblemBody<void, void> res;
    res.getLocalBFIs().add(bfis);
    res.getPBCs().add(pbcs);
    return res;
  }

  template <class RHSNumber>
  inline
  ProblemBody<void, void> operator+(
      const ProblemBody<void, void>& pb,
      const LinearFormIntegratorBase<RHSNumber>& lfi)
  {
    ProblemBody<void, void> res(pb);
    res.getLFIs().add(lfi);
    return res;
  }

  template <class RHSNumber>
  inline
  ProblemBody<void, void> operator+(
      const ProblemBody<void, void>& pb,
      const LocalBilinearFormIntegratorBase<RHSNumber>& lbfi)
  {
    ProblemBody<void, void> res(pb);
    res.getLocalBFIs().add(lbfi);
    return res;
  }

  template <class RHSNumber>
  inline
  ProblemBody<void, void> operator-(
      const ProblemBody<void, void>& pb,
      const LocalBilinearFormIntegratorBase<RHSNumber>& lbfi)
  {
    ProblemBody<void, void> res(pb);
    res.getLocalBFIs().add(UnaryMinus(lbfi));
    return res;
  }

  template <class RHSNumber>
  inline
  ProblemBody<void, void> operator+(
      const ProblemBody<void, void>& pb,
      const GlobalBilinearFormIntegratorBase<RHSNumber>& gbfi)
  {
    ProblemBody<void, void> res(pb);
    res.getGlobalBFIs().add(gbfi);
    return res;
  }

  template <class OperatorType, class VectorType, class RHSNumber>
  ProblemBody<OperatorType, VectorType> operator+(
      const ProblemBody<OperatorType, VectorType>& pb,
      const LinearFormIntegratorBase<RHSNumber>& lfi)
  {
    ProblemBody<OperatorType, VectorType> res(pb);
    res.getLFIs().add(lfi);
    return res;
  }

  template <class OperatorType, class VectorType, class RHSNumber>
  ProblemBody<OperatorType, VectorType> operator+(
      const ProblemBody<OperatorType, VectorType>& pb,
      const FormLanguage::List<LinearFormIntegratorBase<RHSNumber>>& lfis)
  {
    ProblemBody<OperatorType, VectorType> res(pb);
    res.getLFIs().add(lfis);
    return res;
  }

  template <class OperatorType, class VectorType, class RHSNumber>
  ProblemBody<OperatorType, VectorType> operator-(
      const ProblemBody<OperatorType, VectorType>& pb,
      const LinearFormIntegratorBase<RHSNumber>& lfi)
  {
    ProblemBody<OperatorType, VectorType> res(pb);
    res.getLFIs().add(UnaryMinus(lfi));
    return res;
  }

  template <class OperatorType, class VectorType, class RHSNumber>
  ProblemBody<OperatorType, VectorType> operator-(
      const ProblemBody<OperatorType, VectorType>& pb,
      const FormLanguage::List<LinearFormIntegratorBase<RHSNumber>>& lfis)
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

  template <class LHSNumber, class OperatorType>
  ProblemBody<OperatorType, void> operator+(
      const LocalBilinearFormIntegratorBase<LHSNumber>& bfi, const BilinearFormBase<OperatorType>& bf)
  {
    ProblemBody<OperatorType, void> res;
    res.getLocalBFIs().add(bfi);
    res.getBFs().add(bf);
    return res;
  }

  template <class OperatorType, class Number>
  ProblemBody<OperatorType, void> operator-(
      const BilinearFormBase<OperatorType>& bf, const LinearFormIntegratorBase<Number>& lfi)
  {
    ProblemBody<OperatorType, void> res;
    res.getBFs().add(bf);
    res.getLFIs().add(UnaryMinus(lfi));
    return res;
  }

  template <class OperatorType, class Number>
  ProblemBody<OperatorType, void> operator-(
      const FormLanguage::List<BilinearFormBase<OperatorType>>& bfs,
      const LinearFormIntegratorBase<Number>& lfi)
  {
    ProblemBody<OperatorType, void> res;
    res.getBFs().add(bfs);
    res.getLFIs().add(UnaryMinus(lfi));
    return res;
  }

  template <class LHSNumber, class OperatorType>
  ProblemBody<OperatorType, void> operator+(
      const LocalBilinearFormIntegratorBase<LHSNumber>& bfi,
      const FormLanguage::List<BilinearFormBase<OperatorType>>& bfs)
  {
    ProblemBody<OperatorType, void> res;
    res.getLocalBFIs().add(bfi);
    res.getBFs().add(bfs);
    return res;
  }
}

#endif
