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

namespace Rodin::FormLanguage
{
  template <class LHSDerived, class RHSDerived, class FESType, Variational::ShapeFunctionSpaceType SpaceType>
  struct Traits<
    Variational::Sum<
      Variational::ShapeFunctionBase<LHSDerived, FESType, SpaceType>,
      Variational::ShapeFunctionBase<RHSDerived, FESType, SpaceType>>>
  {
    using FES = FESType;
    static constexpr Variational::ShapeFunctionSpaceType Space = SpaceType;
  };
}

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
  class Sum<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>> final
    : public FunctionBase<Sum<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>
  {
    public:
      using LHS = FunctionBase<LHSDerived>;
      using RHS = FunctionBase<RHSDerived>;
      using Parent = FunctionBase<Sum<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>;
      using LHSRange = typename FormLanguage::Traits<LHS>::RangeType;
      using RHSRange = typename FormLanguage::Traits<RHS>::RangeType;
      static_assert(std::is_same_v<LHSRange, RHSRange>);

      constexpr
      Sum(const LHS& lhs, const RHS& rhs)
        : m_lhs(lhs.copy()), m_rhs(rhs.copy())
      {
        assert(lhs.getRangeShape() == rhs.getRangeShape());
      }

      constexpr
      Sum(const Sum& other)
        : Parent(other),
          m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
      {}

      constexpr
      Sum(Sum&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
      {}

      inline
      constexpr
      RangeShape getRangeShape() const
      {
        assert(getLHS().getRangeShape() == getLHS().getRangeShape());
        return getLHS().getRangeShape();
      }

      inline
      constexpr
      const auto& getLHS() const
      {
        assert(m_lhs);
        return *m_lhs;
      }

      inline
      constexpr
      const auto& getRHS() const
      {
        assert(m_rhs);
        return *m_rhs;
      }

      inline
      constexpr
      Sum& traceOf(Geometry::Attribute& attr)
      {
        Parent::traceOf(attr);
        getLHS().traceOf(attr);
        getRHS().traceOf(attr);
        return *this;
      }

      inline
      constexpr
      Sum& traceOf(const std::set<Geometry::Attribute>& attrs)
      {
        Parent::traceOf(attrs);
        getLHS().traceOf(attrs);
        getRHS().traceOf(attrs);
        return *this;
      }

      inline
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return this->object(getLHS().getValue(p)) + this->object(getRHS().getValue(p));
      }

      inline
      constexpr
      void getValueByReference(Math::Vector& res, const Geometry::Point& p) const
      {
        static_assert(FormLanguage::IsVectorRange<LHSRange>::Value);
        getLHS().getValue(res, p);
        res += getRHS().getValue(p);
      }

      inline
      constexpr
      void getValueByReference(Math::Matrix& res, const Geometry::Point& p) const
      {
        static_assert(FormLanguage::IsMatrixRange<LHSRange>::Value);
        getLHS().getValue(res, p);
        res += getRHS().getValue(p);
      }

      inline Sum* copy() const noexcept override
      {
        return new Sum(*this);
      }

    private:
      std::unique_ptr<LHS> m_lhs;
      std::unique_ptr<RHS> m_rhs;
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
  template <class LHSDerived, class RHSDerived, class FESType, ShapeFunctionSpaceType SpaceType>
  class Sum<ShapeFunctionBase<LHSDerived, FESType, SpaceType>, ShapeFunctionBase<RHSDerived, FESType, SpaceType>> final
    : public ShapeFunctionBase<Sum<ShapeFunctionBase<LHSDerived, FESType, SpaceType>, ShapeFunctionBase<RHSDerived, FESType, SpaceType>>>
  {
    public:
      using FES = FESType;
      static constexpr ShapeFunctionSpaceType Space = SpaceType;

      using LHS = ShapeFunctionBase<LHSDerived, FES, Space>;
      using RHS = ShapeFunctionBase<RHSDerived, FES, Space>;
      using Parent = ShapeFunctionBase<Sum<LHS, RHS>, FES, Space>;
      using LHSRange = typename FormLanguage::Traits<LHS>::RangeType;
      using RHSRange = typename FormLanguage::Traits<RHS>::RangeType;
      static_assert(std::is_same_v<LHSRange, RHSRange>);

      constexpr
      Sum(const LHS& lhs, const RHS& rhs)
        : Parent(lhs.getFiniteElementSpace()),
          m_lhs(lhs.copy()), m_rhs(rhs.copy())
      {
        assert(lhs.getRangeShape() == rhs.getRangeShape());
        assert(lhs.getLeaf().getUUID() == rhs.getLeaf().getUUID());
        assert(lhs.getFiniteElementSpace() == rhs.getFiniteElementSpace());
      }

      constexpr
      Sum(const Sum& other)
        : Parent(other),
          m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
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
        assert(m_lhs);
        return *m_lhs;
      }

      inline
      constexpr
      const RHS& getRHS() const
      {
        assert(m_rhs);
        return *m_rhs;
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
        assert(getLHS().getRangeShape() == getRHS().getRangeShape());
        return getLHS().getRangeShape();
      }

      inline
      constexpr
      size_t getDOFs(const Geometry::Polytope& element) const
      {
        assert(getLHS().getDOFs(element) == getRHS().getDOFs(element));
        return getLHS().getDOFs(element);
      }

      inline
      constexpr
      auto getTensorBasis(const Geometry::Point& p) const
      {
        const auto& lhs = this->object(getLHS().getTensorBasis(p));
        const auto& rhs = this->object(getRHS().getTensorBasis(p));
        return lhs + rhs;
      }

      inline
      constexpr
      const auto& getFiniteElementSpace() const
      {
        return getLHS().getFiniteElementSpace();
      }

      inline Sum* copy() const noexcept override
      {
        return new Sum(*this);
      }

    private:
      std::unique_ptr<LHS> m_lhs;
      std::unique_ptr<RHS> m_rhs;
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
