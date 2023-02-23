/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FORMLANGUAGE_SUM_H
#define RODIN_VARIATIONAL_FORMLANGUAGE_SUM_H

#include "Rodin/FormLanguage/Base.h"
#include "Rodin/FormLanguage/List.h"

#include "ForwardDecls.h"
#include "Function.h"
#include "ShapeFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup SumSpecializations Sum Template Specializations
   * @brief Template specializations of the Sum class.
   * @see Sum
   */

  /**
   * @ingroup SumSpecializations
   */
  template <class LHSDerived, class RHSDerived>
  class Sum<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>
    : public FunctionBase<Sum<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>
  {
    public:
      using LHS = FunctionBase<LHSDerived>;
      using RHS = FunctionBase<RHSDerived>;
      using Parent = FunctionBase<Sum<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>;

      Sum(const LHS& lhs, const RHS& rhs)
        : m_lhs(lhs), m_rhs(rhs)
      {
        assert(lhs.getRangeShape() == rhs.getRangeShape());
      }

      Sum(const Sum& other)
        : Parent(other),
          m_lhs(other.m_lhs), m_rhs(other.m_rhs)
      {}

      Sum(Sum&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
      {}

      inline
      constexpr
      RangeShape getRangeShape() const
      {
        assert(m_lhs.getRangeShape() == m_rhs.getRangeShape());
        return m_lhs.getRangeShape();
      }

      inline
      constexpr
      Sum& traceOf(Geometry::Attribute attrs)
      {
        m_lhs.traceOf(attrs);
        m_rhs.traceOf(attrs);
        return *this;
      }

      inline
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return m_lhs.getValue(p) + m_rhs.getValue(p);
      }

    private:
      LHS m_lhs;
      RHS m_rhs;
  };

  template <class LHSDerived, class RHSDerived>
  Sum(const FunctionBase<LHSDerived>&, const FunctionBase<RHSDerived>&)
    -> Sum<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>;

  template <class LHSDerived, class RHSDerived>
  inline
  constexpr
  auto
  operator+(const FunctionBase<LHSDerived>& lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return Sum(lhs, rhs);
  }

  template <class LHSDerived, class Number, typename = std::enable_if_t<std::is_arithmetic_v<Number>>>
  inline
  constexpr
  auto
  operator+(const FunctionBase<LHSDerived>& lhs, Number rhs)
  {
    return Sum(lhs, ScalarFunction(rhs));
  }

  template <class Number, class RHSDerived, typename = std::enable_if_t<std::is_arithmetic_v<Number>>>
  inline
  constexpr
  auto
  operator+(Number lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return Sum(ScalarFunction(lhs), rhs);
  }

  /**
   * @ingroup SumSpecializations
   */
  template <class LHSDerived, class RHSDerived, class FES, ShapeFunctionSpaceType Space>
  class Sum<ShapeFunctionBase<LHSDerived, FES, Space>, ShapeFunctionBase<RHSDerived, FES, Space>> final
    : public ShapeFunctionBase<Sum<ShapeFunctionBase<LHSDerived, FES, Space>, ShapeFunctionBase<RHSDerived, FES, Space>>, FES, Space>
  {
    public:
      using LHS = ShapeFunctionBase<LHSDerived, FES, Space>;
      using RHS = ShapeFunctionBase<RHSDerived, FES, Space>;
      using Parent = ShapeFunctionBase<Sum<LHS, RHS>, FES, Space>;

      constexpr
      Sum(const LHS& lhs, const RHS& rhs)
        : m_lhs(lhs), m_rhs(rhs)
      {
        assert(lhs.getRangeShape() == rhs.getRangeShape());
        assert(lhs.getLeaf().getUUID() == rhs.getLeaf().getUUID());
      }

      constexpr
      Sum(const Sum& other)
        : Parent(other),
          m_lhs(other.m_lhs), m_rhs(other.m_rhs)
      {}

      constexpr
      Sum(Sum&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
      {}

      inline
      constexpr
      const LHS& getLHS() const
      {
        return m_lhs;
      }

      inline
      constexpr
      const RHS& getRHS() const
      {
        return m_rhs;
      }

      inline
      constexpr
      const auto& getLeaf() const
      {
        return getRHS().getLeaf();
      }

      inline
      constexpr
      RangeShape getRangeShape() const
      {
        assert(m_lhs.getRangeShape() == m_rhs.getRangeShape());
        return m_lhs.getRangeShape();
      }

      inline
      constexpr
      size_t getDOFs(const Geometry::Simplex& element) const
      {
        assert(getLHS().getDOFs(element) == getRHS().getDOFs(element));
        return getLHS().getDOFs(element);
      }

      inline
      constexpr
      auto getTensorBasis(const Geometry::Point& p) const
      {
        return getLHS().getTensorBasis( p) + getRHS().getTensorBasis(p);
      }

      inline
      constexpr
      auto& getFiniteElementSpace()
      {
        return getLHS().getFiniteElementSpace();
      }

      inline
      constexpr
      const auto& getFiniteElementSpace() const
      {
        return getLHS().getFiniteElementSpace();
      }

    private:
      LHS m_lhs;
      RHS m_rhs;
  };

  template <class LHSDerived, class RHSDerived, class FES, ShapeFunctionSpaceType Space>
  Sum(const ShapeFunctionBase<LHSDerived, FES, Space>&, const ShapeFunctionBase<RHSDerived, FES, Space>&)
    -> Sum<ShapeFunctionBase<LHSDerived, FES, Space>, ShapeFunctionBase<RHSDerived, FES, Space>>;

  template <class LHSDerived, class RHSDerived, class FES, ShapeFunctionSpaceType Space>
  inline
  constexpr
  auto
  operator+(const ShapeFunctionBase<LHSDerived, FES, Space>& lhs,
            const ShapeFunctionBase<RHSDerived, FES, Space>& rhs)
  {
    return Sum(lhs, rhs);
  }

  FormLanguage::List<BilinearFormIntegratorBase>
  operator+(
      const BilinearFormIntegratorBase& lhs,
      const BilinearFormIntegratorBase& rhs);

  FormLanguage::List<BilinearFormIntegratorBase>
  operator+(
      const BilinearFormIntegratorBase& lhs,
      const FormLanguage::List<BilinearFormIntegratorBase>& rhs);

  FormLanguage::List<BilinearFormIntegratorBase>
  operator+(
      const FormLanguage::List<BilinearFormIntegratorBase>& lhs,
      const BilinearFormIntegratorBase& rhs);

  FormLanguage::List<BilinearFormIntegratorBase>
  operator+(
      const FormLanguage::List<BilinearFormIntegratorBase>& lhs,
      const FormLanguage::List<BilinearFormIntegratorBase>& rhs);

  FormLanguage::List<LinearFormIntegratorBase>
  operator+(
      const LinearFormIntegratorBase& lhs,
      const LinearFormIntegratorBase& rhs);

  FormLanguage::List<LinearFormIntegratorBase>
  operator+(
      const LinearFormIntegratorBase& lhs,
      const FormLanguage::List<LinearFormIntegratorBase>& rhs);

  FormLanguage::List<LinearFormIntegratorBase>
  operator+(
      const FormLanguage::List<LinearFormIntegratorBase>& lhs,
      const LinearFormIntegratorBase& rhs);

  FormLanguage::List<LinearFormIntegratorBase>
  operator+(
      const FormLanguage::List<LinearFormIntegratorBase>& lhs,
      const FormLanguage::List<LinearFormIntegratorBase>& rhs);

}

#endif
