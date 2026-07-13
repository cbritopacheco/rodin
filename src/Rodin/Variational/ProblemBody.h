/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ProblemBody.h
 * @brief Problem body class managing integrators and boundary conditions.
 *
 * This file defines the ProblemBody class, which manages the collection of
 * integrators (bilinear and linear forms) and boundary conditions that
 * comprise a variational problem. The ProblemBody serves as a container
 * and coordinator for all components of a finite element formulation.
 *
 * ## Role in Problem Assembly
 * ProblemBody acts as an intermediate layer that:
 * 1. Collects all bilinear form integrators (matrix terms)
 * 2. Collects all linear form integrators (load vector terms)
 * 3. Manages essential (Dirichlet) boundary conditions
 * 4. Manages periodic boundary conditions
 * 5. Coordinates the assembly process
 *
 * ## Design Pattern
 * The separation of ProblemBody from Problem follows the Bridge pattern,
 * allowing:
 * - Shared state between different problem types
 * - Flexible composition of problem components
 * - Independent evolution of assembly and solution strategies
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
   * @ingroup RodinVariational
   * @brief Base class representing the body of a variational problem.
   *
   * ProblemBodyBase manages all integrators and boundary conditions for a
   * variational problem, providing a unified interface for problem assembly.
   * It serves as the container for all mathematical components of the
   * discrete system.
   *
   * @tparam Scalar Scalar type for problem coefficients
   */
  template <class Scalar>
  class ProblemBodyBase : public FormLanguage::Base
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = Scalar;

      /// @brief Linear form integrator base type.
      using LinearFormIntegratorBaseType = LinearFormIntegratorBase<ScalarType>;

      /// @brief Local bilinear form integrator base type.
      using LocalBilinearFormIntegratorBaseType = LocalBilinearFormIntegratorBase<ScalarType>;

      /// @brief Global bilinear form integrator base type.
      using GlobalBilinearFormIntegratorBaseType = GlobalBilinearFormIntegratorBase<ScalarType>;

      /// @brief Linear form integrator list type.
      using LinearFormIntegratorBaseListType = FormLanguage::List<LinearFormIntegratorBaseType>;

      /// @brief Local bilinear form integrator list type.
      using LocalBilinearFormIntegratorBaseListType = FormLanguage::List<LocalBilinearFormIntegratorBaseType>;

      /// @brief Global bilinear form integrator list type.
      using GlobalBilinearFormIntegratorBaseListType = FormLanguage::List<GlobalBilinearFormIntegratorBaseType>;

      /// @brief Essential boundary condition collection type.
      using EssentialBoundaryType = EssentialBoundary<ScalarType>;

      /// @brief Periodic boundary condition collection type.
      using PeriodicBoundaryType = PeriodicBoundary<ScalarType>;

      /// @brief Parent class type.
      using Parent = FormLanguage::Base;

      /// @brief Default constructor.
      ProblemBodyBase() = default;

      /// @brief Copy constructor.
      ProblemBodyBase(const ProblemBodyBase& other)
        : Parent(other),
          m_lfis(other.m_lfis),
          m_lbfis(other.m_lbfis),
          m_gbfis(other.m_gbfis),
          m_essBdr(other.m_essBdr),
          m_periodicBdr(other.m_periodicBdr)
      {}

      /// @brief Copy assignment operator.
      ProblemBodyBase& operator=(const ProblemBodyBase& other)
      {
        if (this != &other)
        {
          m_lfis = other.m_lfis;
          m_lbfis = other.m_lbfis;
          m_gbfis = other.m_gbfis;
          m_essBdr = other.m_essBdr;
          m_periodicBdr = other.m_periodicBdr;
        }
        return *this;
      }

      /// @brief Move constructor.
      ProblemBodyBase(ProblemBodyBase&& other)
        : Parent(std::move(other)),
          m_lfis(std::move(other.m_lfis)),
          m_lbfis(std::move(other.m_lbfis)),
          m_gbfis(std::move(other.m_gbfis)),
          m_essBdr(std::move(other.m_essBdr)),
          m_periodicBdr(std::move(other.m_periodicBdr))
      {}

      /// @brief Move assignment operator.
      ProblemBodyBase& operator=(ProblemBodyBase&& other)
      {
        m_lfis = std::move(other.m_lfis);
        m_lbfis = std::move(other.m_lbfis);
        m_gbfis = std::move(other.m_gbfis);
        m_essBdr = std::move(other.m_essBdr);
        m_periodicBdr = std::move(other.m_periodicBdr);
        return *this;
      }

      /// @brief Returns periodic boundary conditions.
      PeriodicBoundaryType& getPBCs()
      {
        return m_periodicBdr;
      }

      /// @brief Returns essential boundary conditions.
      EssentialBoundaryType& getDBCs()
      {
        return m_essBdr;
      }

      /// @brief Returns local bilinear form integrators.
      LocalBilinearFormIntegratorBaseListType& getLocalBFIs()
      {
        return m_lbfis;
      }

      /// @brief Returns global bilinear form integrators.
      GlobalBilinearFormIntegratorBaseListType& getGlobalBFIs()
      {
        return m_gbfis;
      }

      /// @brief Returns linear form integrators.
      LinearFormIntegratorBaseListType& getLFIs()
      {
        return m_lfis;
      }

      /// @brief Returns periodic boundary conditions.
      const PeriodicBoundaryType& getPBCs() const
      {
        return m_periodicBdr;
      }

      /// @brief Returns essential boundary conditions.
      const EssentialBoundaryType& getDBCs() const
      {
        return m_essBdr;
      }

      /// @brief Returns linear form integrators.
      const LinearFormIntegratorBaseListType& getLFIs() const
      {
        return m_lfis;
      }

      /// @brief Returns local bilinear form integrators.
      const LocalBilinearFormIntegratorBaseListType& getLocalBFIs() const
      {
        return m_lbfis;
      }

      /// @brief Returns global bilinear form integrators.
      const GlobalBilinearFormIntegratorBaseListType& getGlobalBFIs() const
      {
        return m_gbfis;
      }

      /// @brief Polymorphically copies this problem body base.
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

  /// @brief Problem body containing only integrators and boundary conditions.
  template <class Scalar>
  class ProblemBody<void, void, Scalar>
    : public ProblemBodyBase<Scalar>
  {
    public:
      /// @brief Parent class type.
      using Parent = ProblemBodyBase<Scalar>;

      /// @brief Default constructor.
      ProblemBody() = default;

      /// @brief Copy constructor.
      ProblemBody(const ProblemBody& other)
        : Parent(other)
      {}

      /// @brief Move constructor.
      ProblemBody(ProblemBody&& other)
        : Parent(std::move(other))
      {}

      /// @brief Copy assignment operator.
      ProblemBody& operator=(const ProblemBody& other)
      {
        if (this != &other)
        {
          Parent::operator=(other);
        }
        return *this;
      }

      /// @brief Move assignment operator.
      ProblemBody& operator=(ProblemBody&& other)
      {
        if (this != &other)
        {
          Parent::operator=(std::move(other));
        }
        return *this;
      }

      /// @brief Polymorphically copies this problem body.
      virtual ProblemBody* copy() const noexcept override
      {
        return new ProblemBody(*this);
      }
  };

  /// @brief Problem body containing operator terms only.
  template <class Operator, class Scalar>
  class ProblemBody<Operator, void, Scalar>
    : public ProblemBodyBase<Scalar>
  {
    public:
      /// @brief Assembled operator type.
      using OperatorType = Operator;

      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<OperatorType>::ScalarType;

      /// @brief Bilinear form base type.
      using BilinearFormBaseType = BilinearFormBase<OperatorType>;

      /// @brief Bilinear form list type.
      using BilinearFormBaseListType = FormLanguage::List<BilinearFormBaseType>;

      /// @brief Parent class type.
      using Parent = ProblemBodyBase<Scalar>;

      /// @brief Default constructor.
      ProblemBody() = default;

      /// @brief Copy constructor.
      ProblemBody(const ProblemBody& other)
        : Parent(other),
          m_bfs(other.m_bfs)
      {}

      /// @brief Move constructor.
      ProblemBody(ProblemBody&& other)
        : Parent(std::move(other)),
          m_bfs(std::move(other.m_bfs))
      {}

      /// @brief Copy assignment operator.
      ProblemBody& operator=(const ProblemBody& other)
      {
        if (this != &other)
        {
          Parent::operator=(other);
          m_bfs = other.m_bfs;
        }
        return *this;
      }

      /// @brief Move assignment operator.
      ProblemBody& operator=(ProblemBody&& other)
      {
        if (this != &other)
        {
          Parent::operator=(std::move(other));
          m_bfs = std::move(other.m_bfs);
        }
        return *this;
      }

      /// @brief Returns bilinear forms.
      BilinearFormBaseListType& getBFs()
      {
        return m_bfs;
      }

      /// @brief Returns bilinear forms.
      const BilinearFormBaseListType& getBFs() const
      {
        return m_bfs;
      }

      /// @brief Polymorphically copies this problem body.
      virtual ProblemBody* copy() const noexcept override
      {
        return new ProblemBody(*this);
      }

    private:
      BilinearFormBaseListType m_bfs;
  };

  /// @brief Problem body containing vector terms only.
  template <class Vector, class Scalar>
  class ProblemBody<void, Vector, Scalar> : public ProblemBodyBase<Scalar>
  {
    public:
      /// @brief Vector type of the linear system.
      using VectorType = Vector;

      /// @brief Linear form base type.
      using LinearFormBaseType = LinearFormBase<VectorType>;

      /// @brief Linear form list type.
      using LinearFormBaseListType = FormLanguage::List<LinearFormBaseType>;

      /// @brief Parent class type.
      using Parent = ProblemBodyBase<Scalar>;

      /// @brief Default constructor.
      ProblemBody() = default;

      /// @brief Copy constructor.
      ProblemBody(const ProblemBody& other)
        : Parent(other),
          m_lfs(other.m_lfs)
      {}

      /// @brief Move constructor.
      ProblemBody(ProblemBody&& other)
        : Parent(std::move(other)),
          m_lfs(std::move(other.m_lfs))
      {}

      /// @brief Copy assignment operator.
      ProblemBody& operator=(const ProblemBody& other)
      {
        if (this != &other)
        {
          Parent::operator=(other);
          m_lfs = other.m_lfs;
        }
        return *this;
      }

      /// @brief Move assignment operator.
      ProblemBody& operator=(ProblemBody&& other)
      {
        if (this != &other)
        {
          Parent::operator=(std::move(other));
          m_lfs = std::move(other.m_lfs);
        }
        return *this;
      }

      /// @brief Returns linear forms.
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
      /// @brief Vector type of the linear system.
      using VectorType = Vector;

      /// @brief Assembled operator type.
      using OperatorType = Operator;

      using VectorScalarType =
        typename FormLanguage::Traits<
          std::remove_reference_t<VectorType>>::ScalarType;

      using OperatorScalarType =
        typename FormLanguage::Traits<
          std::remove_reference_t<OperatorType>>::ScalarType;

      using LinearFormBaseType = LinearFormBase<VectorType>;

      using BilinearFormBaseType = BilinearFormBase<OperatorType>;

      using LinearFormBaseListType = FormLanguage::List<LinearFormBaseType>;

      using BilinearFormBaseListType = FormLanguage::List<BilinearFormBaseType>;

      /// @brief Linear form integrator base type.
      using LinearFormIntegratorBaseType = LinearFormIntegratorBase<VectorScalarType>;

      using LocalBilinearFormIntegratorBaseType = LocalBilinearFormIntegratorBase<OperatorScalarType>;

      using GlobalBilinearFormIntegratorBaseType = GlobalBilinearFormIntegratorBase<OperatorScalarType>;

      using LinearFormIntegratorBaseListType = FormLanguage::List<LinearFormIntegratorBaseType>;

      using LocalBilinearFormIntegratorBaseListType = FormLanguage::List<LocalBilinearFormIntegratorBaseType>;

      using GlobalBilinearFormIntegratorBaseListType = FormLanguage::List<GlobalBilinearFormIntegratorBaseType>;

      /// @brief Parent class type.
      using Parent = ProblemBodyBase<Scalar>;

      ProblemBody() = default;

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

      ProblemBody(const BilinearFormBaseType& bf)
      {
        m_bfs.add(bf);
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

      ProblemBody& operator=(const ProblemBody& other)
      {
        if (this != &other)
        {
          Parent::operator=(other);
          m_lfs = other.m_lfs;
          m_bfs = other.m_bfs;
        }
        return *this;
      }

      ProblemBody& operator=(ProblemBody&& other)
      {
        if (this != &other)
        {
          Parent::operator=(std::move(other));
          m_lfs = std::move(other.m_lfs);
          m_bfs = std::move(other.m_bfs);
        }
        return *this;
      }

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
    /// @brief Scalar value type.
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
    /// @brief Scalar value type.
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
    /// @brief Scalar value type.
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
    /// @brief Scalar value type.
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
    /// @brief Scalar value type.
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
    /// @brief Scalar value type.
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
    /// @brief Scalar value type.
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
    /// @brief Scalar value type.
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
    /// @brief Scalar value type.
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
    /// @brief Scalar value type.
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
    /// @brief Scalar value type.
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
    /// @brief Scalar value type.
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
    /// @brief Scalar value type.
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
    /// @brief Scalar value type.
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
    /// @brief Scalar value type.
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
    /// @brief Scalar value type.
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
    /// @brief Scalar value type.
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
    /// @brief Scalar value type.
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
    /// @brief Scalar value type.
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
    /// @brief Scalar value type.
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
    /// @brief Scalar value type.
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
    /// @brief Scalar value type.
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
    /// @brief Scalar value type.
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
    /// @brief Scalar value type.
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
    /// @brief Scalar value type.
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
    /// @brief Scalar value type.
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
    /// @brief Scalar value type.
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
    /// @brief Scalar value type.
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
    /// @brief Scalar value type.
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
    using RHSScalar = typename FormLanguage::Traits<std::remove_reference_t<OperatorType>>::ScalarType;
    /// @brief Scalar value type.
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
    /// @brief Scalar value type.
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
    using RHSScalar = typename FormLanguage::Traits<std::remove_reference_t<OperatorType>>::ScalarType;
    /// @brief Scalar value type.
    using ScalarType = typename FormLanguage::Sum<LHSScalar, RHSScalar>::Type;
    ProblemBody<OperatorType, void, ScalarType> res;
    res.getLocalBFIs().add(bfi);
    res.getBFs().add(bf);
    return res;
  }

  /**
   * @brief Combines a preassembled BilinearForm with a DirichletBC into a
   * ProblemBody.
   *
   * Enables expression chains such as:
   * @code
   *   problem = preassembledBF + DirichletBC(u, g);
   * @endcode
   * where the BilinearForm has already been assembled and the Dirichlet
   * boundary condition is added to constrain the problem.
   */
  template <class OperatorType, class RHSScalar>
  auto
  operator+(
      const BilinearFormBase<OperatorType>& bf, const DirichletBCBase<RHSScalar>& dbc)
  {
    /// @brief Scalar value type.
    using ScalarType = RHSScalar;
    ProblemBody<OperatorType, void, ScalarType> res;
    res.getBFs().add(bf);
    res.getDBCs().add(dbc);
    return res;
  }

  template <class OperatorType, class RHSScalar>
  auto
  operator-(
      const BilinearFormBase<OperatorType>& bf, const LinearFormIntegratorBase<RHSScalar>& lfi)
  {
    using LHSScalar = typename FormLanguage::Traits<std::remove_reference_t<OperatorType>>::ScalarType;
    /// @brief Scalar value type.
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
    using LHSScalar = typename FormLanguage::Traits<std::remove_reference_t<OperatorType>>::ScalarType;
    /// @brief Scalar value type.
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
    using RHSScalar = typename FormLanguage::Traits<std::remove_reference_t<OperatorType>>::ScalarType;
    /// @brief Scalar value type.
    using ScalarType = typename FormLanguage::Sum<LHSScalar, RHSScalar>::Type;
    ProblemBody<OperatorType, void, ScalarType> res;
    res.getLocalBFIs().add(bfi);
    res.getBFs().add(bfs);
    return res;
  }

  /**
   * @brief Combines a list of local bilinear form integrators with a
   * preassembled BilinearForm into a ProblemBody.
   *
   * Enables expression chains such as:
   * @code
   *   problem = Integral(u, v) - Integral(p, v) + preassembledBF - ...;
   * @endcode
   * where the first two Integral terms produce a List<LocalBFI> (via Sum)
   * and the preassembled BilinearForm is added afterwards.
   */
  template <class LHSScalar, class OperatorType>
  auto
  operator+(
      const FormLanguage::List<LocalBilinearFormIntegratorBase<LHSScalar>>& lbfis,
      const BilinearFormBase<OperatorType>& bf)
  {
    using RHSScalar = typename FormLanguage::Traits<
      std::remove_reference_t<OperatorType>>::ScalarType;
    /// @brief Scalar value type.
    using ScalarType = typename FormLanguage::Sum<LHSScalar, RHSScalar>::Type;
    ProblemBody<OperatorType, void, ScalarType> res;
    res.getLocalBFIs().add(lbfis);
    res.getBFs().add(bf);
    return res;
  }

  /**
   * @brief Subtracts a preassembled LinearForm from a single local bilinear
   * form integrator to create a ProblemBody.
   *
   * Enables expression chains such as:
   * @code
   *   problem = Integral(Grad(u), Grad(v)) - assembledLoadForm;
   * @endcode
   * where a single Integral term (LocalBFI) is combined with a preassembled
   * LinearForm.
   */
  template <class LHSScalar, class VectorType>
  auto
  operator-(
      const LocalBilinearFormIntegratorBase<LHSScalar>& bfi,
      const LinearFormBase<VectorType>& lf)
  {
    using RHSScalar = typename FormLanguage::Traits<
      std::remove_reference_t<VectorType>>::ScalarType;
    /// @brief Scalar value type.
    using ScalarType = typename FormLanguage::Minus<LHSScalar, RHSScalar>::Type;
    ProblemBody<void, VectorType, ScalarType> res;
    res.getLocalBFIs().add(bfi);
    res.getLFs().add(lf);
    return res;
  }

  /**
   * @brief Subtracts a preassembled LinearForm from a list of local bilinear
   * form integrators to create a ProblemBody.
   *
   * Enables expression chains such as:
   * @code
   *   problem = Integral(u, v) - Integral(p, v) + Integral(p, q) - loadP0;
   * @endcode
   * where the Integral terms produce a List<LocalBFI> and the preassembled
   * LinearForm is subtracted afterwards.
   */
  template <class LHSScalar, class VectorType>
  auto
  operator-(
      const FormLanguage::List<LocalBilinearFormIntegratorBase<LHSScalar>>& lbfis,
      const LinearFormBase<VectorType>& lf)
  {
    using RHSScalar = typename FormLanguage::Traits<
      std::remove_reference_t<VectorType>>::ScalarType;
    /// @brief Scalar value type.
    using ScalarType = typename FormLanguage::Minus<LHSScalar, RHSScalar>::Type;
    ProblemBody<void, VectorType, ScalarType> res;
    res.getLocalBFIs().add(lbfis);
    res.getLFs().add(lf);
    return res;
  }

  /**
   * @brief Subtracts a preassembled LinearForm from a ProblemBody that has
   * an OperatorType but no VectorType yet.
   *
   * Enables expression chains such as:
   * @code
   *   problem = Integral(u, v) - Integral(p, v) + preassembledBF - loadP0;
   * @endcode
   * where the chain first produces ProblemBody<Op, void, S> (after adding
   * the preassembled BF) and then the LinearForm introduces the VectorType.
   */
  template <class OperatorType, class LHSScalar, class VectorType>
  auto
  operator-(
      const ProblemBody<OperatorType, void, LHSScalar>& pb,
      const LinearFormBase<VectorType>& lf)
  {
    using RHSScalar = typename FormLanguage::Traits<
      std::remove_reference_t<VectorType>>::ScalarType;
    /// @brief Scalar value type.
    using ScalarType = typename FormLanguage::Minus<LHSScalar, RHSScalar>::Type;
    ProblemBody<OperatorType, VectorType, ScalarType> res(pb);
    res.getLFs().add(lf);
    return res;
  }

  template <class LHSScalar, class RHSScalar>
  auto operator+(
      const FormLanguage::List<LinearFormIntegratorBase<LHSScalar>>& lfis,
      const LocalBilinearFormIntegratorBase<RHSScalar>& bfi)
  {
    /// @brief Scalar value type.
    using ScalarType = typename FormLanguage::Sum<LHSScalar, RHSScalar>::Type;

    ProblemBody<void, void, ScalarType> res;
    res.getLFIs().add(lfis);
    res.getLocalBFIs().add(bfi);
    return res;
  }

  template <class LHSScalar, class RHSScalar>
  auto operator-(
      const FormLanguage::List<LinearFormIntegratorBase<LHSScalar>>& lfis,
      const LocalBilinearFormIntegratorBase<RHSScalar>& bfi)
  {
    /// @brief Scalar value type.
    using ScalarType = typename FormLanguage::Minus<LHSScalar, RHSScalar>::Type;

    ProblemBody<void, void, ScalarType> res;
    res.getLFIs().add(lfis);
    res.getLocalBFIs().add(UnaryMinus(bfi));
    return res;
  }

  template <class Operator, class Vector, class LHSScalar, class RHSScalar>
  auto operator+(
      const ProblemBody<Operator, Vector, LHSScalar>& pb,
      const FormLanguage::List<LocalBilinearFormIntegratorBase<RHSScalar>>& bfis)
  {
    /// @brief Scalar value type.
    using ScalarType = typename FormLanguage::Sum<LHSScalar, RHSScalar>::Type;

    ProblemBody<Operator, Vector, ScalarType> res(pb);
    res.getLocalBFIs().add(bfis);
    return res;
  }

  template <class Operator, class Vector, class LHSScalar, class RHSScalar>
  auto operator-(
      const ProblemBody<Operator, Vector, LHSScalar>& pb,
      const FormLanguage::List<LocalBilinearFormIntegratorBase<RHSScalar>>& bfis)
  {
    /// @brief Scalar value type.
    using ScalarType = typename FormLanguage::Minus<LHSScalar, RHSScalar>::Type;

    ProblemBody<Operator, Vector, ScalarType> res(pb);
    res.getLocalBFIs().add(UnaryMinus(bfis));
    return res;
  }
}

#endif
