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

#include "LinearFormIntegrator.h"
#include "BilinearFormIntegrator.h"

namespace Rodin::FormLanguage
{
  template <class LHSDerived, class RHSDerived>
  struct Traits<
    Variational::Sum<Variational::FunctionBase<LHSDerived>, Variational::FunctionBase<RHSDerived>>>
  {
    using LHSType = Variational::FunctionBase<LHSDerived>;
    using RHSType = Variational::FunctionBase<RHSDerived>;
  };

  template <class LHSDerived, class RHSDerived, class FES, Variational::ShapeFunctionSpaceType Space>
  struct Traits<
    Variational::Sum<
      Variational::ShapeFunctionBase<LHSDerived, FES, Space>,
      Variational::ShapeFunctionBase<RHSDerived, FES, Space>>>
  {
    using FESType = FES;
    using LHSType = Variational::FunctionBase<LHSDerived>;
    using RHSType = Variational::FunctionBase<RHSDerived>;
    static constexpr Variational::ShapeFunctionSpaceType SpaceType = Space;
  };

  template <class LHSNumber, class RHSNumber>
  struct Traits<
    Variational::Sum<
      Variational::LinearFormIntegratorBase<LHSNumber>,
      Variational::LinearFormIntegratorBase<RHSNumber>>>
  {
    using LHSScalarType = LHSNumber;

    using RHSScalarType = RHSNumber;

    using LHSType = Variational::LinearFormIntegratorBase<LHSScalarType>;

    using RHSType = Variational::LinearFormIntegratorBase<RHSScalarType>;

    using ScalarType = decltype(std::declval<LHSScalarType>() + std::declval<RHSScalarType>());
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
      using LHSType = FunctionBase<LHSDerived>;
      using RHSType = FunctionBase<RHSDerived>;

      using LHSRangeType = typename FormLanguage::Traits<LHSType>::RangeType;
      using RHSRangeType = typename FormLanguage::Traits<RHSType>::RangeType;

      using Parent = FunctionBase<Sum<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>;
      static_assert(std::is_same_v<LHSRangeType, RHSRangeType>);

      constexpr
      Sum(const LHSType& lhs, const RHSType& rhs)
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
      Sum& traceOf(Geometry::Attribute attr)
      {
        Parent::traceOf(attr);
        getLHS().traceOf(attr);
        getRHS().traceOf(attr);
        return *this;
      }

      inline
      constexpr
      Sum& traceOf(const FlatSet<Geometry::Attribute>& attrs)
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
      void getValue(Math::Vector<Real>& res, const Geometry::Point& p) const
      {
        static_assert(FormLanguage::IsVectorRange<LHSRangeType>::Value);
        getLHS().getValue(res, p);
        res += getRHS().getValue(p);
      }

      inline
      constexpr
      void getValue(Math::Matrix<Real>& res, const Geometry::Point& p) const
      {
        static_assert(FormLanguage::IsMatrixRange<LHSRangeType>::Value);
        getLHS().getValue(res, p);
        res += getRHS().getValue(p);
      }

      inline Sum* copy() const noexcept override
      {
        return new Sum(*this);
      }

    private:
      std::unique_ptr<LHSType> m_lhs;
      std::unique_ptr<RHSType> m_rhs;
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
    return Sum(lhs, RealFunction(rhs));
  }

  template <class Number, class RHSDerived, typename = std::enable_if_t<std::is_arithmetic_v<Number>>>
  inline
  constexpr
  auto
  operator+(Number lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return Sum(RealFunction(lhs), rhs);
  }

  /**
   * @ingroup SumSpecializations
   */
  template <class LHSDerived, class RHSDerived, class FES, ShapeFunctionSpaceType Space>
  class Sum<ShapeFunctionBase<LHSDerived, FES, Space>, ShapeFunctionBase<RHSDerived, FES, Space>> final
    : public ShapeFunctionBase<Sum<ShapeFunctionBase<LHSDerived, FES, Space>, ShapeFunctionBase<RHSDerived, FES, Space>>>
  {
    public:
      using FESType = FES;
      static constexpr ShapeFunctionSpaceType SpaceType = Space;

      using LHSType = ShapeFunctionBase<LHSDerived, FES, Space>;

      using RHSType = ShapeFunctionBase<RHSDerived, FES, Space>;

      using LHSRangeType = typename FormLanguage::Traits<LHSType>::RangeType;

      using RHSRangeType = typename FormLanguage::Traits<RHSType>::RangeType;

      using Parent = ShapeFunctionBase<Sum<LHSType, RHSType>, FES, Space>;
      static_assert(std::is_same_v<LHSRangeType, RHSRangeType>);

      constexpr
      Sum(const LHSType& lhs, const RHSType& rhs)
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
      const LHSType& getLHS() const
      {
        assert(m_lhs);
        return *m_lhs;
      }

      inline
      constexpr
      const RHSType& getRHS() const
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
      Sum& setPoint(const Geometry::Point& p)
      {
        m_lhs->setPoint(p);
        m_rhs->setPoint(p);
        return *this;
      }

      inline
      constexpr
      auto getBasis(size_t local) const
      {
        return this->object(getLHS().getBasis(local)) + this->object(getRHS().getBasis(local));
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
      std::unique_ptr<LHSType> m_lhs;
      std::unique_ptr<RHSType> m_rhs;
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

  template <class LHSNumber, class RHSNumber>
  class Sum<LinearFormIntegratorBase<LHSNumber>, LinearFormIntegratorBase<RHSNumber>>
    : public FormLanguage::List<
              LinearFormIntegratorBase<decltype(std::declval<LHSNumber>() + std::declval<RHSNumber>())>>
  {
    public:
      using LHSScalarType = LHSNumber;

      using RHSScalarType = RHSNumber;

      using LHSType = LinearFormIntegratorBase<LHSScalarType>;

      using RHSType = LinearFormIntegratorBase<RHSScalarType>;

      using ScalarType = decltype(std::declval<LHSScalarType>() + std::declval<RHSScalarType>());

      using Parent = FormLanguage::List<LinearFormIntegratorBase<ScalarType>>;

      Sum(const LHSType& lhs, const RHSType& rhs)
      {
        this->add(lhs);
        this->add(rhs);
      }

      Sum(const Sum& other)
        : Parent(other)
      {}

      Sum(Sum&& other)
        : Parent(std::move(other))
      {}
  };

  template <class LHSNumber, class RHSNumber>
  Sum(const LinearFormIntegratorBase<LHSNumber>&, const LinearFormIntegratorBase<RHSNumber>&)
    -> Sum<LinearFormIntegratorBase<LHSNumber>, LinearFormIntegratorBase<RHSNumber>>;

  template <class LHSNumber, class RHSNumber>
  inline
  constexpr
  auto
  operator+(
      const LinearFormIntegratorBase<LHSNumber>& lhs,
      const LinearFormIntegratorBase<RHSNumber>& rhs)
  {
    return Sum(lhs, rhs);
  }

  template <class LHSNumber, class RHSNumber>
  class Sum<LinearFormIntegratorBase<LHSNumber>, FormLanguage::List<LinearFormIntegratorBase<RHSNumber>>>
    : public FormLanguage::List<
              LinearFormIntegratorBase<decltype(std::declval<LHSNumber>() + std::declval<RHSNumber>())>>
  {
    public:
      using LHSScalarType = LHSNumber;

      using RHSScalarType = RHSNumber;

      using LHSType = LinearFormIntegratorBase<LHSScalarType>;

      using RHSType = FormLanguage::List<LinearFormIntegratorBase<RHSNumber>>;

      using ScalarType = decltype(std::declval<LHSScalarType>() + std::declval<RHSScalarType>());

      using Parent = FormLanguage::List<LinearFormIntegratorBase<ScalarType>>;

      Sum(const LHSType& lhs, const RHSType& rhs)
      {
        this->add(lhs);
        this->add(rhs);
      }

      Sum(const Sum& other)
        : Parent(other)
      {}

      Sum(Sum&& other)
        : Parent(std::move(other))
      {}
  };

  template <class LHSNumber, class RHSNumber>
  Sum(const LinearFormIntegratorBase<LHSNumber>&, const FormLanguage::List<LinearFormIntegratorBase<RHSNumber>&>)
    -> Sum<LinearFormIntegratorBase<LHSNumber>, FormLanguage::List<LinearFormIntegratorBase<RHSNumber>>>;

  template <class LHSNumber, class RHSNumber>
  inline
  constexpr
  auto
  operator+(
      const LinearFormIntegratorBase<LHSNumber>& lhs,
      const FormLanguage::List<LinearFormIntegratorBase<RHSNumber>>& rhs)
  {
    return Sum(lhs, rhs);
  }

  template <class LHSNumber, class RHSNumber>
  class Sum<
    FormLanguage::List<LinearFormIntegratorBase<LHSNumber>>, LinearFormIntegratorBase<RHSNumber>>
    : public FormLanguage::List<
              LinearFormIntegratorBase<decltype(std::declval<LHSNumber>() + std::declval<RHSNumber>())>>
  {
    public:
      using LHSScalarType = LHSNumber;

      using RHSScalarType = RHSNumber;

      using LHSType = FormLanguage::List<LinearFormIntegratorBase<LHSNumber>>;

      using RHSType = LinearFormIntegratorBase<RHSScalarType>;

      using ScalarType = decltype(std::declval<LHSScalarType>() + std::declval<RHSScalarType>());

      using Parent = FormLanguage::List<LinearFormIntegratorBase<ScalarType>>;

      Sum(const LHSType& lhs, const RHSType& rhs)
      {
        this->add(lhs);
        this->add(rhs);
      }

      Sum(const Sum& other)
        : Parent(other)
      {}

      Sum(Sum&& other)
        : Parent(std::move(other))
      {}
  };

  template <class LHSNumber, class RHSNumber>
  Sum(const FormLanguage::List<LinearFormIntegratorBase<RHSNumber>>& lhs,
      const LinearFormIntegratorBase<LHSNumber>& rhs)
    -> Sum<FormLanguage::List<LinearFormIntegratorBase<LHSNumber>>, LinearFormIntegratorBase<RHSNumber>>;

  template <class LHSNumber, class RHSNumber>
  inline
  constexpr
  auto
  operator+(
      const FormLanguage::List<LinearFormIntegratorBase<RHSNumber>>& lhs,
      const LinearFormIntegratorBase<LHSNumber>& rhs)
  {
    return Sum(lhs, rhs);
  }

  template <class LHSNumber, class RHSNumber>
  class Sum<
    FormLanguage::List<LinearFormIntegratorBase<LHSNumber>>,
    FormLanguage::List<LinearFormIntegratorBase<RHSNumber>>>
      : public FormLanguage::List<
          LinearFormIntegratorBase<decltype(std::declval<LHSNumber>() + std::declval<RHSNumber>())>>
  {
    public:
      using LHSScalarType = LHSNumber;

      using RHSScalarType = RHSNumber;

      using LHSType = FormLanguage::List<LinearFormIntegratorBase<RHSNumber>>;

      using RHSType = FormLanguage::List<LinearFormIntegratorBase<RHSNumber>>;

      using ScalarType = decltype(std::declval<LHSScalarType>() + std::declval<RHSScalarType>());

      using Parent = FormLanguage::List<LinearFormIntegratorBase<ScalarType>>;

      Sum(const LHSType& lhs, const RHSType& rhs)
      {
        this->add(lhs);
        this->add(rhs);
      }

      Sum(const Sum& other)
        : Parent(other)
      {}

      Sum(Sum&& other)
        : Parent(std::move(other))
      {}
  };

  template <class LHSNumber, class RHSNumber>
  Sum(const FormLanguage::List<LinearFormIntegratorBase<LHSNumber>>&,
      const FormLanguage::List<LinearFormIntegratorBase<RHSNumber>>&)
    -> Sum<
        FormLanguage::List<LinearFormIntegratorBase<LHSNumber>>,
        FormLanguage::List<LinearFormIntegratorBase<RHSNumber>>>;

  template <class LHSNumber, class RHSNumber>
  inline
  constexpr
  auto
  operator+(
      const FormLanguage::List<LinearFormIntegratorBase<LHSNumber>>& lhs,
      const FormLanguage::List<LinearFormIntegratorBase<RHSNumber>>& rhs)
  {
    return Sum(lhs, rhs);
  }

  template <class LHSNumber, class RHSNumber>
  class Sum<LocalBilinearFormIntegratorBase<LHSNumber>, LocalBilinearFormIntegratorBase<RHSNumber>>
    : public FormLanguage::List<
              LocalBilinearFormIntegratorBase<decltype(std::declval<LHSNumber>() + std::declval<RHSNumber>())>>
  {
    public:
      using LHSScalarType = LHSNumber;

      using RHSScalarType = RHSNumber;

      using LHSType = LocalBilinearFormIntegratorBase<LHSScalarType>;

      using RHSType = LocalBilinearFormIntegratorBase<RHSScalarType>;

      using ScalarType = decltype(std::declval<LHSScalarType>() + std::declval<RHSScalarType>());

      using Parent = FormLanguage::List<LocalBilinearFormIntegratorBase<ScalarType>>;

      Sum(const LHSType& lhs, const RHSType& rhs)
      {
        this->add(lhs);
        this->add(rhs);
      }

      Sum(const Sum& other)
        : Parent(other)
      {}

      Sum(Sum&& other)
        : Parent(std::move(other))
      {}
  };

  template <class LHSNumber, class RHSNumber>
  Sum(const LocalBilinearFormIntegratorBase<LHSNumber>& lhs,
      const LocalBilinearFormIntegratorBase<RHSNumber>& rhs)
    -> Sum<LocalBilinearFormIntegratorBase<LHSNumber>, LocalBilinearFormIntegratorBase<RHSNumber>>;

  template <class LHSNumber, class RHSNumber>
  inline
  constexpr
  auto
  operator+(
    const LocalBilinearFormIntegratorBase<LHSNumber>& lhs,
    const LocalBilinearFormIntegratorBase<RHSNumber>& rhs)
  {
    return Sum(lhs, rhs);
  }

  template <class LHSNumber, class RHSNumber>
  class Sum<LocalBilinearFormIntegratorBase<LHSNumber>, FormLanguage::List<LocalBilinearFormIntegratorBase<RHSNumber>>>
    : public FormLanguage::List<
              LocalBilinearFormIntegratorBase<decltype(std::declval<LHSNumber>() + std::declval<RHSNumber>())>>
  {
    public:
      using LHSScalarType = LHSNumber;

      using RHSScalarType = RHSNumber;

      using LHSType = LocalBilinearFormIntegratorBase<LHSScalarType>;

      using RHSType = FormLanguage::List<LocalBilinearFormIntegratorBase<RHSNumber>>;

      using ScalarType = decltype(std::declval<LHSScalarType>() + std::declval<RHSScalarType>());

      using Parent = FormLanguage::List<LocalBilinearFormIntegratorBase<ScalarType>>;

      Sum(const LHSType& lhs, const RHSType& rhs)
      {
        this->add(lhs);
        this->add(rhs);
      }

      Sum(const Sum& other)
        : Parent(other)
      {}

      Sum(Sum&& other)
        : Parent(std::move(other))
      {}
  };

  template <class LHSNumber, class RHSNumber>
  Sum(const LocalBilinearFormIntegratorBase<LHSNumber>&,
      const FormLanguage::List<LocalBilinearFormIntegratorBase<RHSNumber>>&)
    -> Sum<
        LocalBilinearFormIntegratorBase<LHSNumber>,
        FormLanguage::List<LocalBilinearFormIntegratorBase<RHSNumber>>>;

  template <class LHSNumber, class RHSNumber>
  inline
  constexpr
  auto
  operator+(
      const LocalBilinearFormIntegratorBase<LHSNumber>& lhs,
      const FormLanguage::List<LocalBilinearFormIntegratorBase<RHSNumber>>& rhs)
  {
    return Sum(lhs, rhs);
  }

  template <class LHSNumber, class RHSNumber>
  class Sum<
    FormLanguage::List<LocalBilinearFormIntegratorBase<LHSNumber>>, LocalBilinearFormIntegratorBase<RHSNumber>>
    : public FormLanguage::List<
              LocalBilinearFormIntegratorBase<decltype(std::declval<LHSNumber>() + std::declval<RHSNumber>())>>
  {
    public:
      using LHSScalarType = LHSNumber;

      using RHSScalarType = RHSNumber;

      using LHSType = FormLanguage::List<LocalBilinearFormIntegratorBase<LHSNumber>>;

      using RHSType = LocalBilinearFormIntegratorBase<RHSScalarType>;

      using ScalarType = decltype(std::declval<LHSScalarType>() + std::declval<RHSScalarType>());

      using Parent = FormLanguage::List<LocalBilinearFormIntegratorBase<ScalarType>>;

      Sum(const LHSType& lhs, const RHSType& rhs)
      {
        this->add(lhs);
        this->add(rhs);
      }

      Sum(const Sum& other)
        : Parent(other)
      {}

      Sum(Sum&& other)
        : Parent(std::move(other))
      {}
  };

  template <class LHSNumber, class RHSNumber>
  Sum(const FormLanguage::List<LocalBilinearFormIntegratorBase<LHSNumber>>&,
      const LocalBilinearFormIntegratorBase<RHSNumber>&)
    -> Sum<
        FormLanguage::List<LocalBilinearFormIntegratorBase<LHSNumber>>,
        LocalBilinearFormIntegratorBase<RHSNumber>>;

  template <class LHSNumber, class RHSNumber>
  inline
  constexpr
  auto
  operator+(
      const FormLanguage::List<LocalBilinearFormIntegratorBase<LHSNumber>>& lhs,
      const LocalBilinearFormIntegratorBase<RHSNumber>& rhs)
  {
    return Sum(lhs, rhs);
  }

  template <class LHSNumber, class RHSNumber>
  class Sum<
    FormLanguage::List<LocalBilinearFormIntegratorBase<LHSNumber>>,
    FormLanguage::List<LocalBilinearFormIntegratorBase<RHSNumber>>>
      : public FormLanguage::List<
          LocalBilinearFormIntegratorBase<decltype(std::declval<LHSNumber>() + std::declval<RHSNumber>())>>
  {
    public:
      using LHSScalarType = LHSNumber;

      using RHSScalarType = RHSNumber;

      using LHSType = FormLanguage::List<LocalBilinearFormIntegratorBase<RHSNumber>>;

      using RHSType = FormLanguage::List<LocalBilinearFormIntegratorBase<RHSNumber>>;

      using ScalarType = decltype(std::declval<LHSScalarType>() + std::declval<RHSScalarType>());

      using Parent = FormLanguage::List<LocalBilinearFormIntegratorBase<ScalarType>>;

      Sum(const LHSType& lhs, const RHSType& rhs)
      {
        this->add(lhs);
        this->add(rhs);
      }

      Sum(const Sum& other)
        : Parent(other)
      {}

      Sum(Sum&& other)
        : Parent(std::move(other))
      {}
  };

  template <class LHSNumber, class RHSNumber>
  Sum(const FormLanguage::List<LocalBilinearFormIntegratorBase<LHSNumber>>&,
      const FormLanguage::List<LocalBilinearFormIntegratorBase<RHSNumber>>&)
    -> Sum<
        FormLanguage::List<LocalBilinearFormIntegratorBase<LHSNumber>>,
        FormLanguage::List<LocalBilinearFormIntegratorBase<RHSNumber>>>;

  template <class LHSNumber, class RHSNumber>
  inline
  constexpr
  auto
  operator+(
      const FormLanguage::List<LocalBilinearFormIntegratorBase<LHSNumber>>& lhs,
      const FormLanguage::List<LocalBilinearFormIntegratorBase<RHSNumber>>& rhs)
  {
    return Sum(lhs, rhs);
  }

  template <class Operator>
  class Sum<BilinearFormBase<Operator>, BilinearFormBase<Operator>>
    : public FormLanguage::List<BilinearFormBase<Operator>>
  {
    public:
      using LHSType = BilinearFormBase<Operator>;

      using RHSType = BilinearFormBase<Operator>;

      using Parent = FormLanguage::List<BilinearFormBase<Operator>>;

      Sum(const LHSType& lhs, const RHSType& rhs)
      {
        this->add(lhs);
        this->add(rhs);
      }

      Sum(const Sum& other)
        : Parent(other)
      {}

      Sum(Sum&& other)
        : Parent(std::move(other))
      {}
  };

  template <class Operator>
  Sum(const BilinearFormBase<Operator>& lhs, const BilinearFormBase<Operator>& rhs)
    -> Sum<BilinearFormBase<Operator>, BilinearFormBase<Operator>>;

  template <class Operator>
  inline
  constexpr
  auto
  operator+(const BilinearFormBase<Operator>& lhs, const BilinearFormBase<Operator>& rhs)
  {
    return Sum(lhs, rhs);
  }
}

#endif
