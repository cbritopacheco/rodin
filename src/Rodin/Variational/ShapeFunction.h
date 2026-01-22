/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2024.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ShapeFunction.h
 * @brief Base class for finite element shape functions and basis functions.
 *
 * This file defines the ShapeFunctionBase and ShapeFunction classes, which 
 * provide the foundation for representing finite element basis functions in
 * both trial and test spaces. These classes form the core of the finite element
 * discretization in Rodin's variational framework.
 */
#ifndef RODIN_VARIATIONAL_SHAPEFUNCTION_H
#define RODIN_VARIATIONAL_SHAPEFUNCTION_H

#include "ForwardDecls.h"

#include "Rodin/Geometry/Point.h"
#include "Rodin/FormLanguage/Base.h"
#include "Rodin/FormLanguage/Traits.h"
#include "Rodin/Variational/Traits.h"

namespace Rodin::FormLanguage
{
  template <class Derived, class FES, Variational::ShapeFunctionSpaceType Space>
  struct Traits<Variational::ShapeFunctionBase<Derived, FES, Space>>
  {
    using DerivedType = Derived;

    using FESType = FES;
    static constexpr const Variational::ShapeFunctionSpaceType SpaceType = Space;

    using ResultType =
      typename ResultOf<Variational::ShapeFunctionBase<Derived, FES, SpaceType>>::Type;

    using RangeType =
      typename RangeOf<Variational::ShapeFunctionBase<Derived, FES, SpaceType>>::Type;

    using ScalarType = typename FormLanguage::Traits<RangeType>::ScalarType;
  };

  template <class Derived, class FES, Variational::ShapeFunctionSpaceType Space>
  struct Traits<Variational::ShapeFunction<Derived, FES, Space>>
  {
    using DerivedType = Derived;

    using FESType = FES;
    static constexpr const Variational::ShapeFunctionSpaceType SpaceType = Space;

    using ResultType =
      typename ResultOf<
        Variational::ShapeFunctionBase<
          Variational::ShapeFunction<Derived, FES, SpaceType>, FES, SpaceType>>::Type;

    using RangeType =
      typename RangeOf<
        Variational::ShapeFunctionBase<
          Variational::ShapeFunction<Derived, FES, SpaceType>, FES, SpaceType>>::Type;

    using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;
  };
}

namespace Rodin::Variational
{
  /**
   * @defgroup ShapeFunctionSpecializations ShapeFunction Template Specializations
   * @brief Template specializations of the ShapeFunction class.
   * @see ShapeFunction
   */

  /**
   * @brief Type trait to check if a type is a TrialFunction.
   * @tparam T Type to check
   */
  template <class T>
  struct IsTrialFunction
  {
    /// @brief Value is true only if T is a TrialFunction specialization
    static constexpr Boolean Value = false;
  };

  /// @brief Specialization for TrialFunction types
  template <class Solution, class FES>
  struct IsTrialFunction<TrialFunction<Solution, FES>>
  {
    static constexpr Boolean Value = true;
  };

  /**
   * @brief Type trait to check if a type is a TestFunction.
   * @tparam T Type to check
   */
  template <class T>
  struct IsTestFunction
  {
    /// @brief Value is true only if T is a TestFunction specialization
    static constexpr Boolean Value = false;
  };

  /// @brief Specialization for TestFunction types
  template <class FES>
  struct IsTestFunction<TestFunction<FES>>
  {
    static constexpr Boolean Value = true;
  };

  /**
   * @ingroup RodinVariational
   * @brief Base class for shape functions in finite element spaces.
   *
   * ShapeFunctionBase provides the foundation for representing finite element
   * basis functions (shape functions) in both trial and test spaces. Shape 
   * functions are the building blocks of finite element approximations.
   *
   * ## Mathematical Foundation
   * In the finite element method, functions are approximated using a linear
   * combination of basis functions (shape functions):
   * @f[
   *   u_h(x) = \sum_{i=1}^N u_i \phi_i(x)
   * @f]
   * where:
   * - @f$ \phi_i(x) @f$ are the shape functions
   * - @f$ u_i @f$ are the degrees of freedom (coefficients)
   * - @f$ N @f$ is the number of basis functions
   *
   * Shape functions have specific properties:
   * - **Local support**: Each @f$ \phi_i @f$ is non-zero only on elements adjacent to node @f$ i @f$
   * - **Partition of unity**: @f$ \sum_{i=1}^N \phi_i(x) = 1 @f$ for all @f$ x @f$
   * - **Nodal property**: For Lagrange elements, @f$ \phi_i(x_j) = \delta_{ij} @f$
   *
   * ## Usage
   * ShapeFunctionBase is a CRTP base class. Derived classes include:
   * - TrialFunction: Represents @f$ u @f$ in the trial space
   * - TestFunction: Represents @f$ v @f$ in the test space
   * - Various operators: Grad, Div, Curl, etc.
   *
   * @tparam Derived Derived class (CRTP pattern)
   * @tparam FES Finite element space type
   * @tparam SpaceType Either TrialSpace or TestSpace
   *
   * @see TrialFunction, TestFunction
   */
  template <
    class Derived,
    class FES = typename FormLanguage::Traits<Derived>::FESType,
    ShapeFunctionSpaceType SpaceType = FormLanguage::Traits<Derived>::SpaceType>
  class ShapeFunctionBase : public FormLanguage::Base
  {
    public:
      /// @brief Finite element space type
      using FESType = FES;

      /// @brief Space type (trial or test)
      static constexpr ShapeFunctionSpaceType Space = SpaceType;

      /// @brief Scalar type from the finite element space
      using ScalarType =
        typename FormLanguage::Traits<FESType>::ScalarType;

      /// @brief Parent class type
      using Parent =
        FormLanguage::Base;

      /**
       * @brief Constructs a shape function in the given finite element space.
       * @param[in] fes Finite element space
       */
      constexpr
      ShapeFunctionBase(const FES& fes)
        : m_fes(fes)
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Shape function to copy
       */
      constexpr
      ShapeFunctionBase(const ShapeFunctionBase& other)
        : Parent(other),
          m_fes(other.m_fes)
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Shape function to move
       */
      constexpr
      ShapeFunctionBase(ShapeFunctionBase&& other)
        : Parent(std::move(other)),
          m_fes(std::move(other.m_fes))
      {}

      /**
       * @brief Indicates whether the shape function is part of a trial or test space.
       * @returns The space type (TrialSpace or TestSpace)
       */
      constexpr
      ShapeFunctionSpaceType getSpaceType() const
      {
        return Space;
      }

      /**
       * @brief Gets the x-component of a vector-valued shape function.
       * @returns Component expression representing the x-component
       */
      auto x() const
      {
        return Component(*this, 0);
      }

      /**
       * @brief Gets the y-component of a vector-valued shape function.
       * @returns Component expression representing the y-component
       */
      auto y() const
      {
        return Component(*this, 1);
      }

      /**
       * @brief Gets the z-component of a vector-valued shape function.
       * @returns Component expression representing the z-component
       */
      auto z() const
      {
        return Component(*this, 2);
      }

      /**
       * @brief Gets the transpose of a matrix-valued shape function.
       * @returns Transpose expression
       */
      constexpr
      auto T() const
      {
        return Transpose(*this);
      }

      /**
       * @brief Gets the operand in the shape function expression.
       * @note CRTP function to be overriden in the Derived class.
       */
      constexpr
      const auto& getLeaf() const
      {
        return static_cast<const Derived&>(*this).getLeaf();
      }

      /**
       * @brief Gets the number of degrees of freedom for the given polytope.
       * @param[in] polytope Polytope (element) for which to count DOFs
       * @returns Number of degrees of freedom
       * @note CRTP function to be overriden in the Derived class.
       */
      constexpr
      size_t getDOFs(const Geometry::Polytope& polytope) const
      {
        return static_cast<const Derived&>(*this).getDOFs(polytope);
      }

      /**
       * @brief Gets the current evaluation point.
       * @returns Reference to the point
       * @note CRTP function to be overriden in the Derived class.
       */
      constexpr
      const Geometry::Point& getPoint() const
      {
        return static_cast<const Derived&>(*this).getPoint();
      }

      /**
       * @brief Sets the evaluation point for the shape function.
       * @param[in] p Point at which to evaluate
       * @returns Reference to the derived object
       * @note CRTP function to be overriden in the Derived class.
       */
      constexpr
      Derived& setPoint(const Geometry::Point& p)
      {
        return static_cast<Derived&>(*this).setPoint(p);
      }

      /**
       * @brief Gets the basis function at the specified local index.
       * @param[in] local Local index of the basis function
       * @returns Expression representing the basis function
       * @note CRTP function to be overriden in the Derived class.
       */
      constexpr
      decltype(auto) getBasis(size_t local) const
      {
        return static_cast<const Derived&>(*this).getBasis(local);
      }

      /**
       * @brief Call operator to get the basis function at the specified index.
       * @param[in] local Local index of the basis function
       * @returns Expression representing the basis function
       *
       * This is a synonym for getBasis(size_t).
       */
      constexpr
      decltype(auto) operator()(size_t local) const
      {
        return static_cast<const Derived&>(*this).getBasis(local);
      }

      /**
       * @brief Gets the finite element space to which the shape function belongs.
       * @returns Const reference to the finite element space
       */
      constexpr
      const FES& getFiniteElementSpace() const
      {
        return m_fes.get();
      }

      /**
       * @brief Gets the derived object (mutable).
       * @returns Reference to derived object
       */
      Derived& getDerived() noexcept
      {
        return static_cast<Derived&>(*this);
      }

      /**
       * @brief Gets the derived object (const).
       * @returns Const reference to derived object
       */
      const Derived& getDerived() const noexcept
      {
        return static_cast<const Derived&>(*this);
      }

      /**
       * @brief Creates a polymorphic copy of the shape function.
       * @returns Pointer to newly allocated copy
       */
      virtual ShapeFunctionBase* copy() const noexcept override
      {
        return static_cast<const Derived&>(*this).copy();
      }

    private:
      std::reference_wrapper<const FES> m_fes;
  };

  /**
   * @ingroup ShapeFunctionSpecializations
   * @brief Intermediate shape function class for CRTP derivation.
   *
   * This class provides an intermediate layer in the ShapeFunction hierarchy,
   * allowing for further specialization while maintaining the common interface.
   * It is used as a base for specific shape function implementations like
   * TestFunction and TrialFunction.
   *
   * @tparam Derived Final derived class
   * @tparam FES Finite element space type
   * @tparam Space Space type (trial or test)
   */
  template <class Derived, class FES, ShapeFunctionSpaceType Space>
  class ShapeFunction
    : public ShapeFunctionBase<ShapeFunction<Derived, FES, Space>, FES, Space>
  {
    public:
      /// @brief Finite element space type
      using FESType = FES;
      
      /// @brief Space type (trial or test)
      static constexpr ShapeFunctionSpaceType SpaceType = Space;

      /// @brief Scalar type from the finite element space
      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      /// @brief Range type from the finite element space
      using RangeType = typename FormLanguage::Traits<FESType>::RangeType;

      /// @brief Parent class type
      using Parent =
        ShapeFunctionBase<
          ShapeFunction<Derived, FESType, SpaceType>, FESType, SpaceType>;

      /// @brief Default constructor is deleted
      ShapeFunction() = delete;

      constexpr
      ShapeFunction(const FESType& fes)
        : Parent(fes),
          m_p(nullptr)
      {}

      constexpr
      ShapeFunction(const ShapeFunction& other)
        : Parent(other),
          m_basis(other.m_basis),
          m_p(nullptr)
      {}

      constexpr
      ShapeFunction(ShapeFunction&& other)
        : Parent(std::move(other)),
          m_basis(std::move(other.m_basis)),
          m_p(std::exchange(other.m_p, nullptr))
      {}

      constexpr
      size_t getDOFs(const Geometry::Polytope& polytope) const
      {
        const size_t d = polytope.getDimension();
        const size_t i = polytope.getIndex();
        return this->getFiniteElementSpace().getFiniteElement(d, i).getCount();
      }

      constexpr
      const Geometry::Point& getPoint() const
      {
        assert(m_p);
        return *m_p;
      }

      ShapeFunction& setPoint(const Geometry::Point& p)
      {
        if (m_p == &p)
          return *this;
        m_p = &p;
        const auto& polytope = p.getPolytope();
        const size_t d = polytope.getDimension();
        const Index i = polytope.getIndex();
        const auto& fes = this->getFiniteElementSpace();
        const auto& fe = fes.getFiniteElement(d, i);
        const size_t count = fe.getCount();
        m_basis.resize(count);
        for (size_t local = 0; local < count; local++)
          m_basis[local] = fes.getPushforward({ d, i }, fe.getBasis(local))(p);
        return *this;
      }

      constexpr
      const RangeType& getBasis(size_t local) const
      {
        assert(local < m_basis.size());
        return m_basis[local];
      }

      constexpr
      const auto& getLeaf() const
      {
        return static_cast<const Derived&>(*this).getLeaf();
      }

      virtual ShapeFunction* copy() const noexcept override
      {
        return static_cast<const Derived&>(*this).copy();
      }

    private:
      std::vector<RangeType> m_basis;

      const Geometry::Point* m_p;
  };
}

#endif
