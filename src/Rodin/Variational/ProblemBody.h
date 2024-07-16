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
  template <class Scalar>
  class ProblemBodyBase : public FormLanguage::Base
  {
    public:
      using ScalarType = Scalar;

      using LinearFormIntegratorBaseType = LinearFormIntegratorBase<ScalarType>;

      using LocalBilinearFormIntegratorBaseType = LocalBilinearFormIntegratorBase<ScalarType>;

      using GlobalBilinearFormIntegratorBaseType = GlobalBilinearFormIntegratorBase<ScalarType>;

      using LinearFormIntegratorBaseListType = FormLanguage::List<LinearFormIntegratorBaseType>;

      using LocalBilinearFormIntegratorBaseListType = FormLanguage::List<LocalBilinearFormIntegratorBaseType>;

      using GlobalBilinearFormIntegratorBaseListType = FormLanguage::List<GlobalBilinearFormIntegratorBaseType>;

      using EssentialBoundaryType = EssentialBoundary<ScalarType>;

      using PeriodicBoundaryType = PeriodicBoundary<ScalarType>;

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

      ProblemBodyBase& operator=(ProblemBodyBase&& other)
      {
        m_lfis = std::move(other.m_lfis);
        m_lbfis = std::move(other.m_lbfis);
        m_gbfis = std::move(other.m_gbfis);
        m_essBdr = std::move(other.m_essBdr);
        m_periodicBdr = std::move(other.m_periodicBdr);
        return *this;
      }

      PeriodicBoundaryType& getPBCs()
      {
        return m_periodicBdr;
      }

      EssentialBoundaryType& getDBCs()
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

      const PeriodicBoundaryType& getPBCs() const
      {
        return m_periodicBdr;
      }

      const EssentialBoundaryType& getDBCs() const
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

      EssentialBoundaryType m_essBdr;
      PeriodicBoundaryType  m_periodicBdr;
  };

  template <class Scalar>
  class ProblemBody<void, void, Scalar>
    : public ProblemBodyBase<Scalar>
  {
    public:
      using Parent = ProblemBodyBase<Scalar>;

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

  template <class Operator, class Scalar>
  class ProblemBody<Operator, void, Scalar>
    : public ProblemBodyBase<Scalar>
  {
    public:
      using OperatorType = Operator;

      using ScalarType = typename FormLanguage::Traits<OperatorType>::ScalarType;

      using BilinearFormBaseType = BilinearFormBase<OperatorType>;

      using BilinearFormBaseListType = FormLanguage::List<BilinearFormBaseType>;

      using Parent = ProblemBodyBase<Scalar>;

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

  template <class Vector, class Scalar>
  class ProblemBody<void, Vector, Scalar> : public ProblemBodyBase<Scalar>
  {
    public:
      using VectorType = Vector;

      using LinearFormBaseType = LinearFormBase<VectorType>;

      using LinearFormBaseListType = FormLanguage::List<LinearFormBaseType>;

      using Parent = ProblemBodyBase<Scalar>;

      ProblemBody() = default;

      ProblemBody(const ProblemBody& other)
        : Parent(other),
          m_lfs(other.m_lfs)
      {}

      ProblemBody(ProblemBody&& other)
        : Parent(std::move(other)),
          m_lfs(std::move(other.m_lfs))
      {}

      LinearFormBaseListType& getLFs()
      {
        return m_lfs;
      }

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

  template <class Operator, class Vector, class Scalar>
  class ProblemBody : public ProblemBodyBase<Scalar>
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

      using Parent = ProblemBodyBase<Scalar>;

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

      ProblemBody(const ProblemBody<OperatorType, void, Scalar>& pbo)
        : Parent(pbo)
      {
        m_bfs.add(pbo.getBFs());
      }

      ProblemBody(const ProblemBody<void, VectorType, Scalar>& pbv)
        : Parent(pbv)
      {
        m_lfs.add(pbv.getLFs());
      }

      ProblemBody(const ProblemBody<void, void, Scalar>& parent)
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

      LinearFormBaseListType& getLFs()
      {
        return m_lfs;
      }

      BilinearFormBaseListType& getBFs()
      {
        return m_bfs;
      }

      const LinearFormBaseListType& getLFs() const
      {
        return m_lfs;
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
      LinearFormBaseListType m_lfs;
      BilinearFormBaseListType m_bfs;
  };

  template <class Scalar>
  ProblemBody(const LocalBilinearFormIntegratorBase<Scalar>&)
    -> ProblemBody<void, void, Scalar>;

  template <class LHSScalar, class RHSScalar>
  auto
  operator+(const LocalBilinearFormIntegratorBase<LHSScalar>& bfi, const LinearFormIntegratorBase<RHSScalar>& lfi)
  {
    using ScalarType = typename FormLanguage::Sum<LHSScalar, RHSScalar>::Type;
    ProblemBody<void, void, ScalarType> res;
    res.getLocalBFIs().add(bfi);
    res.getLFIs().add(lfi);
    return res;
  }

  template <class LHSScalar, class RHSScalar>
  auto
  operator+(const LinearFormIntegratorBase<LHSScalar>& lfi, const LocalBilinearFormIntegratorBase<RHSScalar>& bfi)
  {
    using ScalarType = typename FormLanguage::Sum<LHSScalar, RHSScalar>::Type;
    ProblemBody<void, void, ScalarType> res;
    res.getLocalBFIs().add(bfi);
    res.getLFIs().add(lfi);
    return res;
  }

  template <class LHSScalar, class RHSScalar>
  auto
  operator-(const LocalBilinearFormIntegratorBase<LHSScalar>& bfi, const LinearFormIntegratorBase<RHSScalar>& lfi)
  {
    using ScalarType = typename FormLanguage::Minus<LHSScalar, RHSScalar>::Type;
    ProblemBody<void, void, ScalarType> res;
    res.getLocalBFIs().add(bfi);
    res.getLFIs().add(UnaryMinus(lfi));
    return res;
  }

  template <class LHSScalar, class RHSScalar>
  auto
  operator-(const GlobalBilinearFormIntegratorBase<LHSScalar>& bfi, const LinearFormIntegratorBase<RHSScalar>& lfi)
  {
    using ScalarType = typename FormLanguage::Minus<LHSScalar, RHSScalar>::Type;
    ProblemBody<void, void, ScalarType> res;
    res.getGlobalBFIs().add(bfi);
    res.getLFIs().add(UnaryMinus(lfi));
    return res;
  }

  template <class LHSScalar, class RHSScalar>
  auto
  operator-(
      const GlobalBilinearFormIntegratorBase<LHSScalar>& bfi,
      const FormLanguage::List<LinearFormIntegratorBase<RHSScalar>>& lfis)
  {
    using ScalarType = typename FormLanguage::Minus<LHSScalar, RHSScalar>::Type;
    ProblemBody<void, void, ScalarType> res;
    res.getGlobalBFIs().add(bfi);
    res.getLFIs().add(UnaryMinus(lfis));
    return res;
  }

  template <class LHSScalar, class RHSScalar>
  auto operator-(
    const LinearFormIntegratorBase<LHSScalar>& lfi, const LocalBilinearFormIntegratorBase<RHSScalar>& bfi)
  {
    using ScalarType = typename FormLanguage::Minus<LHSScalar, RHSScalar>::Type;
    ProblemBody<void, void, ScalarType> res;
    res.getLocalBFIs().add(UnaryMinus(bfi));
    res.getLFIs().add(lfi);
    return res;
  }

  template <class LHSScalar, class RHSScalar>
  auto operator+(
    const FormLanguage::List<LocalBilinearFormIntegratorBase<LHSScalar>>& bfis,
    const LinearFormIntegratorBase<RHSScalar>& lfi)
  {
    using ScalarType = typename FormLanguage::Sum<LHSScalar, RHSScalar>::Type;
    ProblemBody<void, void, ScalarType> res;
    res.getLocalBFIs().add(bfis);
    res.getLFIs().add(lfi);
    return res;
  }

  template <class LHSScalar, class RHSScalar>
  auto operator+(
      const FormLanguage::List<LocalBilinearFormIntegratorBase<LHSScalar>>& lbfis,
      const GlobalBilinearFormIntegratorBase<RHSScalar>& gbfi)
  {
    using ScalarType = typename FormLanguage::Sum<LHSScalar, RHSScalar>::Type;
    ProblemBody<void, void, ScalarType> res;
    res.getLocalBFIs().add(lbfis);
    res.getGlobalBFIs().add(gbfi);
    return res;
  }

  template <class LHSScalar, class RHSScalar>
  auto operator+(
      const FormLanguage::List<LocalBilinearFormIntegratorBase<LHSScalar>>& lbfis,
      const FormLanguage::List<GlobalBilinearFormIntegratorBase<RHSScalar>>& gbfis)
  {
    using ScalarType = typename FormLanguage::Sum<LHSScalar, RHSScalar>::Type;
    ProblemBody<void, void, ScalarType> res;
    res.getLocalBFIs().add(lbfis);
    res.getGlobalBFIs().add(gbfis);
    return res;
  }

  template <class LHSScalar, class RHSScalar>
  auto operator+(
      const LocalBilinearFormIntegratorBase<LHSScalar>& lbfi,
      const FormLanguage::List<GlobalBilinearFormIntegratorBase<RHSScalar>>& gbfis)
  {
    using ScalarType = typename FormLanguage::Sum<LHSScalar, RHSScalar>::Type;
    ProblemBody<void, void, ScalarType> res;
    res.getLocalBFIs().add(lbfi);
    res.getGlobalBFIs().add(gbfis);
    return res;
  }

  template <class LHSScalar, class RHSScalar>
  auto operator+(
      const LocalBilinearFormIntegratorBase<LHSScalar>& lbfi,
      const GlobalBilinearFormIntegratorBase<RHSScalar>& gbfi)
  {
    using ScalarType = typename FormLanguage::Sum<LHSScalar, RHSScalar>::Type;
    ProblemBody<void, void, ScalarType> res;
    res.getLocalBFIs().add(lbfi);
    res.getGlobalBFIs().add(gbfi);
    return res;
  }

  template <class LHSScalar, class RHSScalar>
  auto operator-(
      const FormLanguage::List<LocalBilinearFormIntegratorBase<LHSScalar>>& bfis,
      const LinearFormIntegratorBase<RHSScalar>& lfi)
  {
    using ScalarType = typename FormLanguage::Minus<LHSScalar, RHSScalar>::Type;
    ProblemBody<void, void, ScalarType> res;
    res.getLocalBFIs().add(bfis);
    res.getLFIs().add(UnaryMinus(lfi));
    return res;
  }

  template <class LHSScalar, class RHSScalar>
  auto
  operator+(const LocalBilinearFormIntegratorBase<LHSScalar>& bfi, const DirichletBCBase<RHSScalar>& dbc)
  {
    using ScalarType = RHSScalar;
    ProblemBody<void, void, ScalarType> res;
    res.getLocalBFIs().add(bfi);
    res.getDBCs().add(dbc);
    return res;
  }

  template <class LHSScalar, class RHSScalar>
  auto
  operator+(const LocalBilinearFormIntegratorBase<LHSScalar>& bfi, const FormLanguage::List<DirichletBCBase<RHSScalar>>& dbcs)
  {
    using ScalarType = RHSScalar;
    ProblemBody<void, void, ScalarType> res;
    res.getLocalBFIs().add(bfi);
    res.getDBCs().add(dbcs);
    return res;
  }

  template <class LHSScalar, class RHSScalar>
  auto operator+(
      const LocalBilinearFormIntegratorBase<LHSScalar>& bfi, const PeriodicBCBase<RHSScalar>& pbc)
  {
    using ScalarType = RHSScalar;
    ProblemBody<void, void, ScalarType> res;
    res.getLocalBFIs().add(bfi);
    res.getPBCs().add(pbc);
    return res;
  }

  template <class LHSScalar, class RHSScalar>
  auto
  operator+(const LocalBilinearFormIntegratorBase<LHSScalar>& bfi, const FormLanguage::List<PeriodicBCBase<RHSScalar>>& pbcs)
  {
    using ScalarType = RHSScalar;
    ProblemBody<void, void, ScalarType> res;
    res.getLocalBFIs().add(bfi);
    res.getPBCs().add(pbcs);
    return res;
  }

  template <class LHSScalar, class RHSScalar>
  auto operator+(
    const FormLanguage::List<LocalBilinearFormIntegratorBase<LHSScalar>>& bfis, const DirichletBCBase<RHSScalar>& dbc)
  {
    using ScalarType = RHSScalar;
    ProblemBody<void, void, ScalarType> res;
    res.getLocalBFIs().add(bfis);
    res.getDBCs().add(dbc);
    return res;
  }

  template <class LHSScalar, class RHSScalar>
  auto
  operator+(
    const FormLanguage::List<LocalBilinearFormIntegratorBase<LHSScalar>>& bfis,
    const FormLanguage::List<DirichletBCBase<RHSScalar>>& dbcs)
  {
    using ScalarType = RHSScalar;
    ProblemBody<void, void, ScalarType> res;
    res.getLocalBFIs().add(bfis);
    res.getDBCs().add(dbcs);
    return res;
  }

  template <class LHSScalar, class RHSScalar>
  auto operator+(
      const FormLanguage::List<LocalBilinearFormIntegratorBase<LHSScalar>>& bfis,
      const FormLanguage::List<PeriodicBCBase<RHSScalar>>& pbcs)
  {
    using ScalarType = RHSScalar;
    ProblemBody<void, void, ScalarType> res;
    res.getLocalBFIs().add(bfis);
    res.getPBCs().add(pbcs);
    return res;
  }

  template <class Operator, class Vector, class LHSScalar, class RHSScalar>
  auto operator+(
      const ProblemBody<Operator, Vector, LHSScalar>& pb,
      const LinearFormIntegratorBase<RHSScalar>& lfi)
  {
    using ScalarType = typename FormLanguage::Sum<LHSScalar, RHSScalar>::Type;
    ProblemBody<Operator, Vector, ScalarType> res(pb);
    res.getLFIs().add(lfi);
    return res;
  }

  template <class Operator, class Vector, class LHSScalar, class RHSScalar>
  auto operator+(
      const ProblemBody<Operator, Vector, LHSScalar>& pb,
      const LocalBilinearFormIntegratorBase<RHSScalar>& lbfi)
  {
    using ScalarType = typename FormLanguage::Sum<LHSScalar, RHSScalar>::Type;
    ProblemBody<Operator, Vector, ScalarType> res(pb);
    res.getLocalBFIs().add(lbfi);
    return res;
  }

  template <class Operator, class Vector, class LHSScalar, class RHSScalar>
  auto operator-(
      const ProblemBody<Operator, Vector, LHSScalar>& pb,
      const LocalBilinearFormIntegratorBase<RHSScalar>& lbfi)
  {
    using ScalarType = typename FormLanguage::Minus<LHSScalar, RHSScalar>::Type;
    ProblemBody<Operator, Vector, ScalarType> res(pb);
    res.getLocalBFIs().add(UnaryMinus(lbfi));
    return res;
  }

  template <class Operator, class Vector, class LHSScalar, class RHSScalar>
  auto
  operator+(
      const ProblemBody<Operator, Vector, LHSScalar>& pb,
      const GlobalBilinearFormIntegratorBase<RHSScalar>& gbfi)
  {
    using ScalarType = typename FormLanguage::Sum<LHSScalar, RHSScalar>::Type;
    ProblemBody<Operator, Vector, ScalarType> res(pb);
    res.getGlobalBFIs().add(gbfi);
    return res;
  }

  template <class OperatorType, class VectorType, class LHSScalar, class RHSScalar>
  auto
  operator+(
      const ProblemBody<OperatorType, VectorType, LHSScalar>& pb,
      const FormLanguage::List<LinearFormIntegratorBase<RHSScalar>>& lfis)
  {
    using ScalarType = typename FormLanguage::Sum<LHSScalar, RHSScalar>::Type;
    ProblemBody<OperatorType, VectorType, ScalarType> res(pb);
    res.getLFIs().add(lfis);
    return res;
  }

  template <class OperatorType, class VectorType, class LHSScalar, class RHSScalar>
  auto
  operator-(
      const ProblemBody<OperatorType, VectorType, LHSScalar>& pb,
      const LinearFormIntegratorBase<RHSScalar>& lfi)
  {
    using ScalarType = typename FormLanguage::Minus<LHSScalar, RHSScalar>::Type;
    ProblemBody<OperatorType, VectorType, ScalarType> res(pb);
    res.getLFIs().add(UnaryMinus(lfi));
    return res;
  }

  template <class OperatorType, class VectorType, class LHSScalar, class RHSScalar>
  auto
  operator-(
      const ProblemBody<OperatorType, VectorType, LHSScalar>& pb,
      const FormLanguage::List<LinearFormIntegratorBase<RHSScalar>>& lfis)
  {
    using ScalarType = typename FormLanguage::Minus<LHSScalar, RHSScalar>::Type;
    ProblemBody<OperatorType, VectorType, ScalarType> res(pb);
    res.getLFIs().add(UnaryMinus(lfis));
    return res;
  }

  template <class OperatorType, class VectorType, class LHSScalar, class RHSScalar>
  auto
  operator+(
      const ProblemBody<OperatorType, VectorType, LHSScalar>& pb, const DirichletBCBase<RHSScalar>& dbc)
  {
    using ScalarType = RHSScalar;
    ProblemBody<OperatorType, VectorType, ScalarType> res(pb);
    res.getDBCs().add(dbc);
    return res;
  }

  template <class OperatorType, class VectorType, class LHSScalar, class RHSScalar>
  auto
  operator+(
      const ProblemBody<OperatorType, VectorType, LHSScalar>& pb, const FormLanguage::List<DirichletBCBase<RHSScalar>>& dbcs)
  {
    using ScalarType = RHSScalar;
    ProblemBody<OperatorType, VectorType, ScalarType> res(pb);
    res.getEssentialBoundary().add(dbcs);
    return res;
  }

  template <class OperatorType, class VectorType, class LHSScalar, class RHSScalar>
  auto
  operator+(
      const ProblemBody<OperatorType, VectorType, LHSScalar>& pb,
      const PeriodicBCBase<RHSScalar>& pbc)
  {
    using ScalarType = RHSScalar;
    ProblemBody<OperatorType, VectorType, ScalarType> res(pb);
    res.getPBCs().add(pbc);
    return res;
  }

  template <class OperatorType, class VectorType, class LHSScalar>
  auto
  operator+(
      const ProblemBody<OperatorType, VectorType, LHSScalar>& pb,
      const BilinearFormBase<OperatorType>& bf)
  {
    using RHSScalar = typename FormLanguage::Traits<BilinearFormBase<OperatorType>>::ScalarType;
    using ScalarType = typename FormLanguage::Sum<LHSScalar, RHSScalar>::Type;
    ProblemBody<OperatorType, VectorType, ScalarType> res(pb);
    res.getBFs().add(bf);
    return res;
  }

  template <class OperatorType, class VectorType, class LHSScalar, class RHSScalar>
  auto
  operator+(
      const ProblemBody<OperatorType, VectorType, LHSScalar>& pb,
      const FormLanguage::List<PeriodicBCBase<RHSScalar>>& pbcs)
  {
    using ScalarType = RHSScalar;
    ProblemBody<OperatorType, VectorType, ScalarType> res(pb);
    res.getPBCs().add(pbcs);
    return res;
  }

  template <class OperatorType, class LHSScalar>
  auto
  operator+(
      const LocalBilinearFormIntegratorBase<LHSScalar>& bfi, const BilinearFormBase<OperatorType>& bf)
  {
    using RHSScalar = typename FormLanguage::Traits<BilinearFormBase<OperatorType>>::ScalarType;
    using ScalarType = typename FormLanguage::Sum<LHSScalar, RHSScalar>::Type;
    ProblemBody<OperatorType, void, ScalarType> res;
    res.getLocalBFIs().add(bfi);
    res.getBFs().add(bf);
    return res;
  }

  template <class OperatorType, class RHSScalar>
  auto
  operator-(
      const BilinearFormBase<OperatorType>& bf, const LinearFormIntegratorBase<RHSScalar>& lfi)
  {
    using LHSScalar = typename FormLanguage::Traits<BilinearFormBase<OperatorType>>::ScalarType;
    using ScalarType = typename FormLanguage::Minus<LHSScalar, RHSScalar>::Type;
    ProblemBody<OperatorType, void, ScalarType> res;
    res.getBFs().add(bf);
    res.getLFIs().add(UnaryMinus(lfi));
    return res;
  }

  template <class OperatorType, class RHSScalar>
  auto
  operator-(
      const FormLanguage::List<BilinearFormBase<OperatorType>>& bfs,
      const LinearFormIntegratorBase<RHSScalar>& lfi)
  {
    using LHSScalar = typename FormLanguage::Traits<BilinearFormBase<OperatorType>>::ScalarType;
    using ScalarType = typename FormLanguage::Minus<LHSScalar, RHSScalar>::Type;
    ProblemBody<OperatorType, void, ScalarType> res;
    res.getBFs().add(bfs);
    res.getLFIs().add(UnaryMinus(lfi));
    return res;
  }

  template <class LHSScalar, class OperatorType>
  auto
  operator+(
      const LocalBilinearFormIntegratorBase<LHSScalar>& bfi,
      const FormLanguage::List<BilinearFormBase<OperatorType>>& bfs)
  {
    using RHSScalar = typename FormLanguage::Traits<BilinearFormBase<OperatorType>>::ScalarType;
    using ScalarType = typename FormLanguage::Sum<LHSScalar, RHSScalar>::Type;
    ProblemBody<OperatorType, void, ScalarType> res;
    res.getLocalBFIs().add(bfi);
    res.getBFs().add(bfs);
    return res;
  }
}

#endif
