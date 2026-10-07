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

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      ProblemBodyBase(const ProblemBodyBase& other)
        : Parent(other),
          m_lfis(other.m_lfis),
          m_lbfis(other.m_lbfis),
          m_gbfis(other.m_gbfis),
          m_essBdr(other.m_essBdr),
          m_periodicBdr(other.m_periodicBdr)
      {}

      /**
       * @brief Copy assignment operator.
       * @param other Object to copy from.
       * @returns Reference to this object after the operation.
       */
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

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      ProblemBodyBase(ProblemBodyBase&& other)
        : Parent(std::move(other)),
          m_lfis(std::move(other.m_lfis)),
          m_lbfis(std::move(other.m_lbfis)),
          m_gbfis(std::move(other.m_gbfis)),
          m_essBdr(std::move(other.m_essBdr)),
          m_periodicBdr(std::move(other.m_periodicBdr))
      {}

      /**
       * @brief Move assignment operator.
       * @param other Object to move from.
       * @returns Reference to this object after the operation.
       */
      ProblemBodyBase& operator=(ProblemBodyBase&& other)
      {
        m_lfis = std::move(other.m_lfis);
        m_lbfis = std::move(other.m_lbfis);
        m_gbfis = std::move(other.m_gbfis);
        m_essBdr = std::move(other.m_essBdr);
        m_periodicBdr = std::move(other.m_periodicBdr);
        return *this;
      }

      /**
       * @brief Returns periodic boundary conditions.
       * @returns Periodic boundary conditions.
       */
      PeriodicBoundaryType& getPBCs()
      {
        return m_periodicBdr;
      }

      /**
       * @brief Returns essential boundary conditions.
       * @returns Essential boundary conditions.
       */
      EssentialBoundaryType& getDBCs()
      {
        return m_essBdr;
      }

      /**
       * @brief Returns local bilinear form integrators.
       * @returns Local bilinear form integrators.
       */
      LocalBilinearFormIntegratorBaseListType& getLocalBFIs()
      {
        return m_lbfis;
      }

      /**
       * @brief Returns global bilinear form integrators.
       * @returns Global bilinear form integrators.
       */
      GlobalBilinearFormIntegratorBaseListType& getGlobalBFIs()
      {
        return m_gbfis;
      }

      /**
       * @brief Returns linear form integrators.
       * @returns Linear form integrators.
       */
      LinearFormIntegratorBaseListType& getLFIs()
      {
        return m_lfis;
      }

      /**
       * @brief Returns periodic boundary conditions.
       * @returns Periodic boundary conditions.
       */
      const PeriodicBoundaryType& getPBCs() const
      {
        return m_periodicBdr;
      }

      /**
       * @brief Returns essential boundary conditions.
       * @returns Essential boundary conditions.
       */
      const EssentialBoundaryType& getDBCs() const
      {
        return m_essBdr;
      }

      /**
       * @brief Returns linear form integrators.
       * @returns Linear form integrators.
       */
      const LinearFormIntegratorBaseListType& getLFIs() const
      {
        return m_lfis;
      }

      /**
       * @brief Returns local bilinear form integrators.
       * @returns Local bilinear form integrators.
       */
      const LocalBilinearFormIntegratorBaseListType& getLocalBFIs() const
      {
        return m_lbfis;
      }

      /**
       * @brief Returns global bilinear form integrators.
       * @returns Global bilinear form integrators.
       */
      const GlobalBilinearFormIntegratorBaseListType& getGlobalBFIs() const
      {
        return m_gbfis;
      }

      /**
       * @brief Polymorphically copies this problem body base.
       * @returns Pointer to a newly allocated copy; the caller owns the returned object.
       */
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

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      ProblemBody(const ProblemBody& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      ProblemBody(ProblemBody&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Copy assignment operator.
       * @param other Object to copy from.
       * @returns Reference to this object after the operation.
       */
      ProblemBody& operator=(const ProblemBody& other)
      {
        if (this != &other)
        {
          Parent::operator=(other);
        }
        return *this;
      }

      /**
       * @brief Move assignment operator.
       * @param other Object to move from.
       * @returns Reference to this object after the operation.
       */
      ProblemBody& operator=(ProblemBody&& other)
      {
        if (this != &other)
        {
          Parent::operator=(std::move(other));
        }
        return *this;
      }

      /**
       * @brief Polymorphically copies this problem body.
       * @returns Pointer to a newly allocated copy; the caller owns the returned object.
       */
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

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      ProblemBody(const ProblemBody& other)
        : Parent(other),
          m_bfs(other.m_bfs)
      {}

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      ProblemBody(ProblemBody&& other)
        : Parent(std::move(other)),
          m_bfs(std::move(other.m_bfs))
      {}

      /**
       * @brief Copy assignment operator.
       * @param other Object to copy from.
       * @returns Reference to this object after the operation.
       */
      ProblemBody& operator=(const ProblemBody& other)
      {
        if (this != &other)
        {
          Parent::operator=(other);
          m_bfs = other.m_bfs;
        }
        return *this;
      }

      /**
       * @brief Move assignment operator.
       * @param other Object to move from.
       * @returns Reference to this object after the operation.
       */
      ProblemBody& operator=(ProblemBody&& other)
      {
        if (this != &other)
        {
          Parent::operator=(std::move(other));
          m_bfs = std::move(other.m_bfs);
        }
        return *this;
      }

      /**
       * @brief Returns bilinear forms.
       * @returns Bilinear forms.
       */
      BilinearFormBaseListType& getBFs()
      {
        return m_bfs;
      }

      /**
       * @brief Returns bilinear forms.
       * @returns Bilinear forms.
       */
      const BilinearFormBaseListType& getBFs() const
      {
        return m_bfs;
      }

      /**
       * @brief Polymorphically copies this problem body.
       * @returns Pointer to a newly allocated copy; the caller owns the returned object.
       */
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

      /**
       * @brief Retains inline integrators when introducing a preassembled vector.
       * @param other Object to copy from.
       */
      ProblemBody(const ProblemBody<void, void, Scalar>& other)
        : Parent(other)
      {}

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      ProblemBody(const ProblemBody& other)
        : Parent(other),
          m_lfs(other.m_lfs)
      {}

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      ProblemBody(ProblemBody&& other)
        : Parent(std::move(other)),
          m_lfs(std::move(other.m_lfs))
      {}

      /**
       * @brief Copy assignment operator.
       * @param other Object to copy from.
       * @returns Reference to this object after the operation.
       */
      ProblemBody& operator=(const ProblemBody& other)
      {
        if (this != &other)
        {
          Parent::operator=(other);
          m_lfs = other.m_lfs;
        }
        return *this;
      }

      /**
       * @brief Move assignment operator.
       * @param other Object to move from.
       * @returns Reference to this object after the operation.
       */
      ProblemBody& operator=(ProblemBody&& other)
      {
        if (this != &other)
        {
          Parent::operator=(std::move(other));
          m_lfs = std::move(other.m_lfs);
        }
        return *this;
      }

      /**
       * @brief Returns linear forms.
       * @returns Linear forms.
       */
      LinearFormBaseListType& getLFs()
      {
        return m_lfs;
      }

      /**
       * @brief Gets the linear forms of the problem.
       * @returns The linear forms of the problem.
       */
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

  /**
   * @brief Accumulated bilinear forms, linear forms and boundary conditions of a
   * Problem.
   */
  template <class Operator, class Vector, class Scalar>
  class ProblemBody : public ProblemBodyBase<Scalar>
  {
    public:
      /// @brief Vector type of the linear system.
      using VectorType = Vector;

      /// @brief Assembled operator type.
      using OperatorType = Operator;

      /// @brief Scalar value type of the vector.
      using VectorScalarType =
        typename FormLanguage::Traits<
          std::remove_reference_t<VectorType>>::ScalarType;

      /// @brief Scalar value type of the operator.
      using OperatorScalarType =
        typename FormLanguage::Traits<
          std::remove_reference_t<OperatorType>>::ScalarType;

      /// @brief Linear form base type.
      using LinearFormBaseType = LinearFormBase<VectorType>;

      /// @brief Bilinear form base type.
      using BilinearFormBaseType = BilinearFormBase<OperatorType>;

      /// @brief List type of linear forms.
      using LinearFormBaseListType = FormLanguage::List<LinearFormBaseType>;

      /// @brief List type of bilinear forms.
      using BilinearFormBaseListType = FormLanguage::List<BilinearFormBaseType>;

      /// @brief Linear form integrator base type.
      using LinearFormIntegratorBaseType = LinearFormIntegratorBase<VectorScalarType>;

      /// @brief Local bilinear form integrator base type.
      using LocalBilinearFormIntegratorBaseType = LocalBilinearFormIntegratorBase<OperatorScalarType>;

      /// @brief Global bilinear form integrator base type.
      using GlobalBilinearFormIntegratorBaseType = GlobalBilinearFormIntegratorBase<OperatorScalarType>;

      /// @brief List type of linear form integrators.
      using LinearFormIntegratorBaseListType = FormLanguage::List<LinearFormIntegratorBaseType>;

      /// @brief List type of local bilinear form integrators.
      using LocalBilinearFormIntegratorBaseListType = FormLanguage::List<LocalBilinearFormIntegratorBaseType>;

      /// @brief List type of global bilinear form integrators.
      using GlobalBilinearFormIntegratorBaseListType = FormLanguage::List<GlobalBilinearFormIntegratorBaseType>;

      /// @brief Parent class type.
      using Parent = ProblemBodyBase<Scalar>;

      ProblemBody() = default;

      /**
       * @brief Constructs the ProblemBody from the given arguments.
       * @param bfi Bilinear form integrator.
       */
      ProblemBody(const LocalBilinearFormIntegratorBaseType& bfi)
      {
        this->getLocalBFIs().add(bfi);
      }

      /**
       * @brief Constructs the ProblemBody from the given arguments.
       * @param bfi Bilinear form integrator.
       */
      ProblemBody(const GlobalBilinearFormIntegratorBaseType& bfi)
      {
        this->getGlobalBFIs().add(bfi);
      }

      /**
       * @brief Constructs the ProblemBody from the given arguments.
       * @param bfis Bilinear form integrators.
       */
      ProblemBody(const LocalBilinearFormIntegratorBaseListType& bfis)
      {
        this->getLocalBFIs().add(bfis);
      }

      /**
       * @brief Constructs the ProblemBody from the given arguments.
       * @param bfis Bilinear form integrators.
       */
      ProblemBody(const GlobalBilinearFormIntegratorBaseListType& bfis)
      {
        this->getGlobalBFIs().add(bfis);
      }

      /**
       * @brief Constructs the ProblemBody from the given arguments.
       * @param pbo Problem body supplying the operator terms.
       */
      ProblemBody(const ProblemBody<OperatorType, void, Scalar>& pbo)
        : Parent(pbo)
      {
        m_bfs.add(pbo.getBFs());
      }

      /**
       * @brief Constructs the ProblemBody from the given arguments.
       * @param bf Bilinear form.
       */
      ProblemBody(const BilinearFormBaseType& bf)
      {
        m_bfs.add(bf);
      }

      /**
       * @brief Constructs the ProblemBody from the given arguments.
       * @param pbv Problem body supplying the vector terms.
       */
      ProblemBody(const ProblemBody<void, VectorType, Scalar>& pbv)
        : Parent(pbv)
      {
        m_lfs.add(pbv.getLFs());
      }

      /**
       * @brief Constructs the ProblemBody from the given arguments.
       * @param parent Problem body supplying the boundary conditions.
       */
      ProblemBody(const ProblemBody<void, void, Scalar>& parent)
        : Parent(parent)
      {}

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      ProblemBody(const ProblemBody& other)
        : Parent(other),
          m_lfs(other.m_lfs),
          m_bfs(other.m_bfs)
      {}

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      ProblemBody(ProblemBody&& other)
        : Parent(std::move(other)),
          m_lfs(std::move(other.m_lfs)),
          m_bfs(std::move(other.m_bfs))
      {}

      /**
       * @brief Copy assignment.
       * @param other Object to copy from.
       * @returns Reference to this object after the operation.
       */
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

      /**
       * @brief Move assignment.
       * @param other Object to move from.
       * @returns Reference to this object after the operation.
       */
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

      /**
       * @brief Gets the linear forms of the problem.
       * @returns The linear forms of the problem.
       */
      LinearFormBaseListType& getLFs()
      {
        return m_lfs;
      }

      /**
       * @brief Gets the bilinear forms of the problem.
       * @returns The bilinear forms of the problem.
       */
      BilinearFormBaseListType& getBFs()
      {
        return m_bfs;
      }

      /**
       * @brief Gets the linear forms of the problem.
       * @returns The linear forms of the problem.
       */
      const LinearFormBaseListType& getLFs() const
      {
        return m_lfs;
      }

      /**
       * @brief Gets the bilinear forms of the problem.
       * @returns The bilinear forms of the problem.
       */
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

  /**
   * @brief Deduction guide for @c ProblemBody.
   * @param pbo Problem body supplying the operator terms.
   */
  template <class Scalar>
  ProblemBody(const LocalBilinearFormIntegratorBase<Scalar>& pbo)
    -> ProblemBody<void, void, Scalar>;

  /**
   * @brief Adds variational terms to construct a problem body.
   * @param bfi Bilinear form integrator.
   * @param lfi Linear form integrator.
   * @returns Problem body owning the supplied terms with their algebraic signs.
   */
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

  /**
   * @brief Adds variational terms to construct a problem body.
   * @param lfi Linear form integrator.
   * @param bfi Bilinear form integrator.
   * @returns Problem body owning the supplied terms with their algebraic signs.
   */
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

  /**
   * @brief Subtracts variational terms to construct a problem body.
   * @param bfi Bilinear form integrator.
   * @param lfi Linear form integrator.
   * @returns Problem body owning the supplied terms with their algebraic signs.
   */
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

  /**
   * @brief Subtracts variational terms to construct a problem body.
   * @param bfi Bilinear form integrator.
   * @param lfi Linear form integrator.
   * @returns Problem body owning the supplied terms with their algebraic signs.
   */
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

  /**
   * @brief Subtracts variational terms to construct a problem body.
   * @param bfi Bilinear form integrator.
   * @param lfis Linear form integrators.
   * @returns Problem body owning the supplied terms with their algebraic signs.
   */
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

  /**
   * @brief Subtracts variational terms to construct a problem body.
   * @param lfi Linear form integrator.
   * @param bfi Bilinear form integrator.
   * @returns Problem body owning the supplied terms with their algebraic signs.
   */
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

  /**
   * @brief Adds variational terms to construct a problem body.
   * @param bfis Bilinear form integrators.
   * @param lfi Linear form integrator.
   * @returns Problem body owning the supplied terms with their algebraic signs.
   */
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

  /**
   * @brief Adds variational terms to construct a problem body.
   * @param lbfis Local bilinear form integrators.
   * @param gbfi Global bilinear form integrator.
   * @returns Problem body owning the supplied terms with their algebraic signs.
   */
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

  /**
   * @brief Adds variational terms to construct a problem body.
   * @param lbfis Local bilinear form integrators.
   * @param gbfis Global bilinear form integrators.
   * @returns Problem body owning the supplied terms with their algebraic signs.
   */
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

  /**
   * @brief Adds variational terms to construct a problem body.
   * @param lbfi Local bilinear form integrator.
   * @param gbfis Global bilinear form integrators.
   * @returns Problem body owning the supplied terms with their algebraic signs.
   */
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

  /**
   * @brief Adds variational terms to construct a problem body.
   * @param lbfi Local bilinear form integrator.
   * @param gbfi Global bilinear form integrator.
   * @returns Problem body owning the supplied terms with their algebraic signs.
   */
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

  /**
   * @brief Subtracts variational terms to construct a problem body.
   * @param bfis Bilinear form integrators.
   * @param lfi Linear form integrator.
   * @returns Problem body owning the supplied terms with their algebraic signs.
   */
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

  /**
   * @brief Adds variational terms to construct a problem body.
   * @param bfi Bilinear form integrator.
   * @param dbc Dirichlet boundary condition.
   * @returns Problem body owning the supplied terms with their algebraic signs.
   */
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

  /**
   * @brief Adds variational terms to construct a problem body.
   * @param bfi Bilinear form integrator.
   * @param dbcs Dirichlet boundary conditions.
   * @returns Problem body owning the supplied terms with their algebraic signs.
   */
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

  /**
   * @brief Adds variational terms to construct a problem body.
   * @param bfi Bilinear form integrator.
   * @param pbc Periodic boundary condition.
   * @returns Problem body owning the supplied terms with their algebraic signs.
   */
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

  /**
   * @brief Adds variational terms to construct a problem body.
   * @param bfi Bilinear form integrator.
   * @param pbcs Periodic boundary conditions.
   * @returns Problem body owning the supplied terms with their algebraic signs.
   */
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

  /**
   * @brief Adds variational terms to construct a problem body.
   * @param bfis Bilinear form integrators.
   * @param dbc Dirichlet boundary condition.
   * @returns Problem body owning the supplied terms with their algebraic signs.
   */
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

  /**
   * @brief Adds variational terms to construct a problem body.
   * @param bfis Bilinear form integrators.
   * @param dbcs Dirichlet boundary conditions.
   * @returns Problem body owning the supplied terms with their algebraic signs.
   */
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

  /**
   * @brief Adds variational terms to construct a problem body.
   * @param bfis Bilinear form integrators.
   * @param pbcs Periodic boundary conditions.
   * @returns Problem body owning the supplied terms with their algebraic signs.
   */
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

  /**
   * @brief Adds variational terms to construct a problem body.
   * @param pb Problem body to extend.
   * @param lfi Linear form integrator.
   * @returns Problem body owning the supplied terms with their algebraic signs.
   */
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

  /**
   * @brief Adds variational terms to construct a problem body.
   * @param pb Problem body to extend.
   * @param lbfi Local bilinear form integrator.
   * @returns Problem body owning the supplied terms with their algebraic signs.
   */
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

  /**
   * @brief Subtracts variational terms to construct a problem body.
   * @param pb Problem body to extend.
   * @param lbfi Local bilinear form integrator.
   * @returns Problem body owning the supplied terms with their algebraic signs.
   */
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

  /**
   * @brief Adds variational terms to construct a problem body.
   * @param pb Problem body to extend.
   * @param gbfi Global bilinear form integrator.
   * @returns Problem body owning the supplied terms with their algebraic signs.
   */
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

  /**
   * @brief Adds variational terms to construct a problem body.
   * @param pb Problem body to extend.
   * @param lfis Linear form integrators.
   * @returns Problem body owning the supplied terms with their algebraic signs.
   */
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

  /**
   * @brief Subtracts variational terms to construct a problem body.
   * @param pb Problem body to extend.
   * @param lfi Linear form integrator.
   * @returns Problem body owning the supplied terms with their algebraic signs.
   */
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

  /**
   * @brief Subtracts variational terms to construct a problem body.
   * @param pb Problem body to extend.
   * @param lfis Linear form integrators.
   * @returns Problem body owning the supplied terms with their algebraic signs.
   */
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

  /**
   * @brief Adds variational terms to construct a problem body.
   * @param pb Problem body to extend.
   * @param dbc Dirichlet boundary condition.
   * @returns Problem body owning the supplied terms with their algebraic signs.
   */
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

  /**
   * @brief Adds variational terms to construct a problem body.
   * @param pb Problem body to extend.
   * @param dbcs Dirichlet boundary conditions.
   * @returns Problem body owning the supplied terms with their algebraic signs.
   */
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

  /**
   * @brief Adds variational terms to construct a problem body.
   * @param pb Problem body to extend.
   * @param pbc Periodic boundary condition.
   * @returns Problem body owning the supplied terms with their algebraic signs.
   */
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

  /**
   * @brief Adds variational terms to construct a problem body.
   * @param pb Problem body to extend.
   * @param bf Preassembled bilinear form.
   * @returns Problem body owning the supplied terms with their algebraic signs.
   */
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

  /**
   * @brief Adds variational terms to construct a problem body.
   * @param pb Problem body to extend.
   * @param pbcs Periodic boundary conditions.
   * @returns Problem body owning the supplied terms with their algebraic signs.
   */
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

  /**
   * @brief Adds variational terms to construct a problem body.
   * @param bfi Bilinear form integrator.
   * @param bf Preassembled bilinear form.
   * @returns Problem body owning the supplied terms with their algebraic signs.
   */
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
    * @param bf Bilinear form.
    * @param dbc Dirichlet boundary condition.
    * @returns Sum of the operands.
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

  /**
   * @brief Subtracts a linear integrator from a preassembled bilinear form.
   * @param bf Preassembled bilinear form.
   * @param lfi Linear integrator whose negated expression is appended.
   * @returns Problem body owning copies of the supplied terms.
   */
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

  /**
   * @brief Subtracts a linear integrator from preassembled bilinear forms.
   * @param bfs Preassembled bilinear forms.
   * @param lfi Linear integrator whose negated expression is appended.
   * @returns Problem body owning copies of the supplied terms.
   */
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

  /**
   * @brief Combines a local bilinear integrator with preassembled bilinear forms.
   * @param bfi Local bilinear integrator.
   * @param bfs Preassembled bilinear forms.
   * @returns Problem body owning copies of the supplied terms.
   */
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
    * @param lbfis Local bilinear form integrators.
    * @param bf Bilinear form.
    * @returns Sum of the operands.
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
    * @param bfi Bilinear form integrator.
    * @param lf Linear form.
    * @returns Difference of the operands, or the negated operand for the unary overload.
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
    * @param lbfis Local bilinear form integrators.
    * @param lf Linear form.
    * @returns Difference of the operands, or the negated operand for the unary overload.
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
    * @param pb Variational problem to operate on.
    * @param lf Linear form.
    * @returns Difference of the operands, or the negated operand for the unary overload.
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

  /**
   * @brief Subtracts a preassembled LinearForm from a ProblemBody that
   * already carries a VectorType.
   *
   * The @f$ \mathrm{void} @f$-vector overload above covers the case where the
   * LinearForm is what introduces the vector; this one covers a body that
   * already holds linear-form content, so that several preassembled forms can
   * be reused in one chain:
   * @code
   *   problem = preassembledBF + Integral(u, v) - Integral(f, v) - preassembledLF;
   * @endcode
    * @param pb Variational problem to operate on.
    * @param lf Linear form.
    * @returns Difference of the operands, or the negated operand for the unary overload.
   */
  template <class OperatorType, class VectorType, class LHSScalar>
  auto operator-(const ProblemBody<OperatorType, VectorType, LHSScalar>& pb,
    const LinearFormBase<VectorType>& lf)
  {
    using RHSScalar =
      typename FormLanguage::Traits<std::remove_reference_t<VectorType>>::ScalarType;
    /// @brief Scalar value type.
    using ScalarType = typename FormLanguage::Minus<LHSScalar, RHSScalar>::Type;
    ProblemBody<OperatorType, VectorType, ScalarType> res(pb);
    res.getLFs().add(lf);
    return res;
  }

  /**
   * @brief Subtracts an owned snapshot of a preassembled bilinear form.
   * @param pb Problem body to extend.
   * @param bf Preassembled bilinear form whose negated snapshot is appended.
   * @returns Extended problem body owning the negated bilinear form.
   */
  template <class Operator, class Vector, class Scalar>
    requires requires(Operator& op) { op *= -1; }
  auto operator-(
    const ProblemBody<Operator, Vector, Scalar>& pb, const BilinearFormBase<Operator>& bf)
  {
    ProblemBody<Operator, Vector, Scalar> res(pb);
    std::unique_ptr<BilinearFormBase<Operator>> negative(bf.copy());
    negative->getOperator() *= -1;
    res.getBFs().add(*negative);
    return res;
  }

  /**
   * @brief Adds a preassembled linear form to the residual, hence negates its load.
   * @param pb Problem body to extend.
   * @param lf Preassembled linear form whose negated load snapshot is appended.
   * @returns Extended problem body owning the negated load.
   */
  template <class Operator, class Vector, class Scalar>
    requires requires(Vector& vec) { vec *= -1; }
  auto operator+(
    const ProblemBody<Operator, Vector, Scalar>& pb, const LinearFormBase<Vector>& lf)
  {
    ProblemBody<Operator, Vector, Scalar> res(pb);
    std::unique_ptr<LinearFormBase<Vector>> negative(lf.copy());
    negative->getVector() *= -1;
    res.getLFs().add(*negative);
    return res;
  }

  /**
   * @brief Combines linear integrators with a local bilinear integrator.
   * @param lfis Linear integrators.
   * @param bfi Local bilinear integrator.
   * @returns Problem body owning copies of the supplied terms.
   */
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

  /**
   * @brief Subtracts variational terms to construct a problem body.
   * @param lfis Linear form integrators.
   * @param bfi Bilinear form integrator.
   * @returns Problem body owning the supplied terms with their algebraic signs.
   */
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

  /**
   * @brief Adds variational terms to construct a problem body.
   * @param pb Problem body to extend.
   * @param bfis Bilinear form integrators.
   * @returns Problem body owning the supplied terms with their algebraic signs.
   */
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

  /**
   * @brief Subtracts variational terms to construct a problem body.
   * @param pb Problem body to extend.
   * @param bfis Bilinear form integrators.
   * @returns Problem body owning the supplied terms with their algebraic signs.
   */
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
