/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_MULT_H
#define RODIN_VARIATIONAL_MULT_H

#include <memory>
#include <type_traits>

#include "Rodin/Alert.h"
#include "Rodin/FormLanguage/Base.h"


#include "ForwardDecls.h"
#include "Function.h"
#include "RealFunction.h"
#include "MatrixFunction.h"
#include "ShapeFunction.h"

#include "LinearFormIntegrator.h"
#include "BilinearFormIntegrator.h"

namespace Rodin::FormLanguage
{
  template <class LHSDerived, class RHSDerived, class FES, Variational::ShapeFunctionSpaceType Space>
  struct Traits<
    Variational::Mult<
      Variational::FunctionBase<LHSDerived>,
      Variational::ShapeFunctionBase<RHSDerived, FES, Space>>>
  {
    using FESType = FES;
    static constexpr Variational::ShapeFunctionSpaceType SpaceType = Space;

    using LHSType =
      Variational::FunctionBase<LHSDerived>;

    using RHSType =
      Variational::ShapeFunctionBase<RHSDerived, FESType, SpaceType>;

    using LHSRangeType =
      typename FormLanguage::Traits<LHSType>::RangeType;

    using RHSRangeType =
      typename FormLanguage::Traits<RHSType>::RangeType;

    using RangeType =
      std::conditional_t<
        // If
        std::is_same_v<LHSRangeType, Real>,
        // Then <----------------------------------------------- LHS is Real
        RHSRangeType,
        // -------------------------------------------------------------------
        // Else
        std::conditional_t<
          // If
          std::is_same_v<LHSRangeType, Math::Vector<Real>>,
          // Then <--------------------------------------------- LHS is Vector
          std::conditional_t<
            // If
            std::is_same_v<RHSRangeType, Real>,
            // Then
            Math::Vector<Real>,
            // Else
            std::conditional_t<
              // If
              std::is_same_v<RHSRangeType, Math::Vector<Real>>,
              // Then
              void,
              // Else
              std::conditional_t<
                // If
                std::is_same_v<RHSRangeType, Math::Matrix<Real>>,
                // Then
                Math::Matrix<Real>,
                // Else
                void
              >
            >
          >,
          // -----------------------------------------------------------------
          // Else
          std::conditional_t<
            // If
            std::is_same_v<LHSRangeType, Math::Matrix<Real>>,
            // Then <------------------------------------------- LHS is Matrix
            std::conditional_t<
              std::is_same_v<RHSRangeType, Real>,
              Math::Matrix<Real>,
              std::conditional_t<
                std::is_same_v<RHSRangeType, Math::Vector<Real>>,
                Math::Vector<Real>,
                std::conditional_t<
                  std::is_same_v<RHSRangeType, Math::Matrix<Real>>,
                    Math::Matrix<Real>,
                    void
                  >
                >
              >,
            // ---------------------------------------------------------------
            // Else
            void
          >
        >
      >;
  };

  template <class LHSDerived, class RHSDerived, class FES, Variational::ShapeFunctionSpaceType Space>
  struct Traits<
    Variational::Mult<
      Variational::ShapeFunctionBase<LHSDerived, FES, Space>,
      Variational::FunctionBase<RHSDerived>>>
  {
    using FESType = FES;
    static constexpr Variational::ShapeFunctionSpaceType SpaceType = Space;
  };
}

namespace Rodin::Variational
{
  /**
   * @defgroup MultSpecializations Mult Template Specializations
   * @brief Template specializations of the Mult class.
   * @see Mult
   */

  /**
   * @ingroup MultSpecializations
   * @brief Multiplication of two FunctionBase instances.
   */
  template <class LHSDerived, class RHSDerived>
  class Mult<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>> final
    : public FunctionBase<Mult<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>
  {
    public:
      using LHSType = FunctionBase<LHSDerived>;

      using RHSType = FunctionBase<RHSDerived>;

      using LHSRangeType = typename FormLanguage::Traits<LHSType>::RangeType;

      using RHSRangeType = typename FormLanguage::Traits<RHSType>::RangeType;

      using Parent = FunctionBase<Mult<LHSType, RHSType>>;

      Mult(const LHSType& lhs, const RHSType& rhs)
        : m_lhs(lhs.copy()), m_rhs(rhs.copy())
      {
        assert(lhs.getRangeType() == RangeType::Real
            || rhs.getRangeType() == RangeType::Real
            || lhs.getRangeShape().width() == rhs.getRangeShape().height());
      }

      Mult(const Mult& other)
        : Parent(other),
          m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
      {}

      Mult(Mult&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
      {}

      inline
      constexpr
      RangeShape getRangeShape() const
      {
        if constexpr (std::is_same_v<LHSRangeType, Real>)
        {
          return getRHS().getRangeShape();
        }
        else if constexpr (std::is_same_v<RHSRangeType, Real>)
        {
          return getLHS().getRangeShape();
        }
        else if constexpr (std::is_same_v<LHSRangeType, Math::Vector<Real>> || std::is_same_v<LHSRangeType, Math::Matrix<Real>>)
        {
          return getLHS().getRangeShape().product(getRHS().getRangeShape());
        }
        else
        {
          assert(false);
          return { 0, 0 };
        }
      }

      inline
      constexpr
      Mult& traceOf(Geometry::Attribute attr)
      {
        getLHS().traceOf(attr);
        getRHS().traceOf(attr);
        return *this;
      }

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
      auto getValue(const Geometry::Point& p) const
      {
        return this->object(getLHS().getValue(p)) * this->object(getRHS().getValue(p));
      }

      inline
      constexpr
      void getValue(Math::Vector<Real>& out, const Geometry::Point& p) const
      {
        if constexpr (std::is_same_v<LHSRangeType, Real> && std::is_same_v<RHSRangeType, Math::Vector<Real>>)
        {
          getRHS().getDerived().getValue(out, p);
          out *= getLHS().getValue(p);
        }
        else if constexpr (std::is_same_v<LHSRangeType, Math::Vector<Real>> && std::is_same_v<RHSRangeType, Real>)
        {
          getLHS().getDerived().getValue(out, p);
          out *= getRHS().getValue(p);
        }
        else
        {
          out = getValue(p);
        }
      }

      inline
      constexpr
      void getValue(Math::Matrix<Real>& out, const Geometry::Point& p) const
      {
        if constexpr (std::is_same_v<LHSRangeType, Real> && std::is_same_v<RHSRangeType, Math::Matrix<Real>>)
        {
          getRHS().getValue(out, p);
          out *= getLHS().getValue(p);
        }
        else if constexpr (std::is_same_v<LHSRangeType, Math::Matrix<Real>> && std::is_same_v<RHSRangeType, Real>)
        {
          getLHS().getValue(out, p);
          out *= getRHS().getValue(p);
        }
        else
        {
          out = getValue(p);
        }
      }

      inline Mult* copy() const noexcept override
      {
        return new Mult(*this);
      }

    private:
      std::unique_ptr<LHSType> m_lhs;
      std::unique_ptr<RHSType> m_rhs;
  };

  template <class LHSDerived, class RHSDerived>
  Mult(const FunctionBase<LHSDerived>&, const FunctionBase<RHSDerived>&)
    -> Mult<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>;

  template <class LHSDerived, class RHSDerived>
  inline
  constexpr
  auto
  operator*(const FunctionBase<LHSDerived>& lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return Mult(lhs, rhs);
  }

  template <class Number, class RHSDerived, typename = std::enable_if_t<std::is_arithmetic_v<Number>>>
  inline
  constexpr
  auto
  operator*(Number lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return Mult(RealFunction(lhs), rhs);
  }

  template <class Number, class LHSDerived, typename = std::enable_if_t<std::is_arithmetic_v<Number>>>
  inline
  constexpr
  auto
  operator*(const FunctionBase<LHSDerived>& lhs, Number rhs)
  {
    return Mult(lhs, RealFunction(rhs));
  }

  template <class LHSDerived>
  inline
  auto
  operator*(const FunctionBase<LHSDerived>& lhs, std::reference_wrapper<const Math::Matrix<Real>> rhs)
  {
    return Mult(lhs, MatrixFunction(rhs));
  }

  template <class LHSDerived>
  inline
  auto
  operator*(std::reference_wrapper<const Math::Matrix<Real>> lhs, const FunctionBase<LHSDerived>& rhs)
  {
    return Mult(MatrixFunction(lhs), rhs);
  }

  /**
   * @ingroup MultSpecializations
   * @brief Left Multiplication of a ShapeFunctionBase by a FunctionBase
   *
   * Represents the following expression:
   * @f[
   *   f A(u)
   * @f]
   * where @f$ f @f$ is the function (scalar, vector or matrix valued), and
   * $A(u)$ is the shape function (scalar, vector, or matrix valued).
   */
  template <class LHSDerived, class RHSDerived, class FES, ShapeFunctionSpaceType Space>
  class Mult<FunctionBase<LHSDerived>, ShapeFunctionBase<RHSDerived, FES, Space>> final
    : public ShapeFunctionBase<Mult<FunctionBase<LHSDerived>, ShapeFunctionBase<RHSDerived, FES, Space>>>
  {
    public:
      using FESType = FES;
      static constexpr ShapeFunctionSpaceType SpaceType = Space;

      using LHSType = FunctionBase<LHSDerived>;

      using RHSType = ShapeFunctionBase<RHSDerived, FES, SpaceType>;

      using LHSRangeType = typename FormLanguage::Traits<LHSType>::RangeType;

      using RHSRangeType = typename FormLanguage::Traits<RHSType>::RangeType;

      using Parent = ShapeFunctionBase<Mult<LHSType, RHSType>, FES, SpaceType>;

      constexpr
      Mult(const LHSType& lhs, const RHSType& rhs)
        : Parent(rhs.getFiniteElementSpace()),
          m_lhs(lhs.copy()), m_rhs(rhs.copy())
      {}

      constexpr
      Mult(const Mult& other)
        : Parent(other),
          m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
      {}

      constexpr
      Mult(Mult&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
      {}

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
        if constexpr (std::is_same_v<LHSRangeType, Real>)
        {
          return getRHS().getRangeShape();
        }
        else if constexpr (std::is_same_v<RHSRangeType, Real>)
        {
          return getLHS().getRangeShape();
        }
        else if constexpr (
            std::is_same_v<LHSRangeType, Math::Vector<Real>> ||
            std::is_same_v<LHSRangeType, Math::Matrix<Real>>)
        {
          return getLHS().getRangeShape().product(getRHS().getRangeShape());
        }
        else
        {
          assert(false);
          return { 0, 0 };
        }
      }

      inline
      constexpr
      size_t getDOFs(const Geometry::Polytope& element) const
      {
        return getRHS().getDOFs(element);
      }

      inline
      constexpr
      const auto& getFiniteElementSpace() const
      {
        return getRHS().getFiniteElementSpace();
      }

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
      const Geometry::Point& getPoint() const
      {
        return m_rhs->getPoint();
      }

      Mult& setPoint(const Geometry::Point& p)
      {
        m_rhs->setPoint(p);
        return *this;
      }

      inline
      constexpr
      auto getBasis(size_t local) const
      {
        const auto& p = getPoint();
        return this->object(getLHS().getValue(p)) * this->object(getRHS().getBasis(local));
      }

      inline Mult* copy() const noexcept override
      {
        return new Mult(*this);
      }

    private:
      std::unique_ptr<LHSType> m_lhs;
      std::unique_ptr<RHSType> m_rhs;
  };

  template <class LHSDerived, class RHSDerived, class FES, ShapeFunctionSpaceType Space>
  Mult(const FunctionBase<LHSDerived>&, const ShapeFunctionBase<RHSDerived, FES, Space>&)
    -> Mult<FunctionBase<LHSDerived>, ShapeFunctionBase<RHSDerived, FES, Space>>;

  template <class LHSDerived, class RHSDerived, class FES, ShapeFunctionSpaceType Space>
  inline
  constexpr
  auto
  operator*(const FunctionBase<LHSDerived>& lhs, const ShapeFunctionBase<RHSDerived, FES, Space>& rhs)
  {
    return Mult(lhs, rhs);
  }

  template <class RHSDerived, class FES, ShapeFunctionSpaceType Space>
  inline
  constexpr
  auto
  operator*(Real lhs, const ShapeFunctionBase<RHSDerived, FES, Space>& rhs)
  {
    return Mult(RealFunction(lhs), rhs);
  }

  /**
   * @ingroup MultSpecializations
   * @brief Right multiplication of a ShapeFunctionBase by a FunctionBase
   *
   * Represents the following expression:
   * @f[
   *    A(u) f
   * @f]
   * where @f$ f @f$ is the function (scalar, vector or matrix valued), and
   * $A(u)$ is the shape function (scalar, vector, or matrix valued).
   */
  template <class LHSDerived, class RHSDerived, class FES, ShapeFunctionSpaceType Space>
  class Mult<ShapeFunctionBase<LHSDerived, FES, Space>, FunctionBase<RHSDerived>> final
    : public ShapeFunctionBase<Mult<ShapeFunctionBase<LHSDerived, FES, Space>, FunctionBase<RHSDerived>>>
  {
    public:
      using FESType = FES;
      static constexpr ShapeFunctionSpaceType SpaceType = Space;

      using LHSType = ShapeFunctionBase<LHSDerived, FES, SpaceType>;

      using RHSType = FunctionBase<RHSDerived>;

      using LHSRangeType = typename FormLanguage::Traits<LHSType>::RangeType;

      using RHSRangeType = typename FormLanguage::Traits<RHSType>::RangeType;

      using LHSScalarType = typename FormLanguage::Traits<LHSRangeType>::ScalarType;

      using RHSScalarType = typename FormLanguage::Traits<RHSRangeType>::ScalarType;

      using ScalarType = decltype(std::declval<LHSScalarType>() * std::declval<RHSScalarType>());

      using Parent = ShapeFunctionBase<Mult<LHSType, RHSType>, FES, SpaceType>;

      constexpr
      Mult(const LHSType& lhs, const RHSType& rhs)
        : Parent(lhs.getFiniteElementSpace()),
          m_lhs(lhs.copy()), m_rhs(rhs.copy())
      {}

      constexpr
      Mult(const Mult& other)
        : Parent(other),
          m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
      {}

      constexpr
      Mult(Mult&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
      {}

      inline
      constexpr
      const auto& getLeaf() const
      {
        return getLHS().getLeaf();
      }

      inline
      constexpr
      RangeShape getRangeShape() const
      {
        if constexpr (std::is_same_v<LHSRangeType, ScalarType>)
        {
          return getRHS().getRangeShape();
        }
        else if constexpr (std::is_same_v<RHSRangeType, ScalarType>)
        {
          return getLHS().getRangeShape();
        }
        else if constexpr (
            std::is_same_v<LHSRangeType, Math::Vector<ScalarType>> ||
            std::is_same_v<LHSRangeType, Math::Matrix<ScalarType>>)
        {
          return getLHS().getRangeShape().product(getRHS().getRangeShape());
        }
        else
        {
          assert(false);
          return { 0, 0 };
        }
      }

      inline
      constexpr
      size_t getDOFs(const Geometry::Polytope& element) const
      {
        return getLHS().getDOFs(element);
      }

      inline
      constexpr
      const auto& getFiniteElementSpace() const
      {
        return getLHS().getFiniteElementSpace();
      }

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
      const Geometry::Point& getPoint() const
      {
        return m_lhs->getPoint();
      }

      inline
      Mult& setPoint(const Geometry::Point& p)
      {
        m_lhs->setPoint(p);
        return *this;
      }

      inline
      constexpr
      auto getBasis(size_t local) const
      {
        const auto& p = getPoint();
        return this->object(getLHS().getBasis(local)) * this->object(getRHS().getValue(p));
      }

      inline Mult* copy() const noexcept override
      {
        return new Mult(*this);
      }

    private:
      std::unique_ptr<LHSType> m_lhs;
      std::unique_ptr<RHSType> m_rhs;
  };

  template <class LHSDerived, class RHSDerived, class FES, ShapeFunctionSpaceType Space>
  Mult(const ShapeFunctionBase<LHSDerived, FES, Space>&, const FunctionBase<RHSDerived>&)
    -> Mult<ShapeFunctionBase<LHSDerived, FES, Space>, FunctionBase<RHSDerived>>;

  template <class LHSDerived, class RHSDerived, class FES, ShapeFunctionSpaceType Space>
  inline
  constexpr
  auto
  operator*(const ShapeFunctionBase<LHSDerived, FES, Space>& lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return Mult(lhs, rhs);
  }

  template <class LHSDerived, class FES, ShapeFunctionSpaceType Space>
  inline
  constexpr
  auto
  operator*(const ShapeFunctionBase<LHSDerived, FES, Space>& lhs, Real rhs)
  {
    return Mult(lhs, RealFunction(rhs));
  }

  template <class LHS, class Number>
  class Mult<LHS, LocalBilinearFormIntegratorBase<Number>>
    : public LocalBilinearFormIntegratorBase<decltype(std::declval<LHS>() * std::declval<Number>())>
  {
    public:
      using LHSType = LHS;

      using RHSType = LocalBilinearFormIntegratorBase<Number>;

      using ScalarType = decltype(std::declval<LHS>() * std::declval<Number>());

      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      Mult(LHSType lhs, const RHSType& rhs)
        : Parent(rhs),
          m_lhs(lhs), m_rhs(rhs.copy())
      {}

      Mult(const Mult& other)
        : Parent(other),
          m_lhs(other.m_lhs), m_rhs(other.m_rhs->copy())
      {}

      Mult(Mult&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
      {}

      Integrator::Region getRegion() const override
      {
        return getRHS().getRegion();
      }

      const LHSType& getLHS() const
      {
        return m_lhs;
      }

      const RHSType& getRHS() const
      {
        assert(m_rhs);
        return *m_rhs;
      }

      const Geometry::Polytope& getPolytope() const override
      {
        return m_rhs->getPolytope();
      }

      Mult& setPolytope(const Geometry::Polytope& polytope) override
      {
        m_rhs->setPolytope(polytope);
        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) override
      {
        return getLHS() * m_rhs->integrate(tr, te);
      }

      Mult* copy() const noexcept override
      {
        return new Mult(*this);
      }

    private:
      const LHSType m_lhs;
      std::unique_ptr<RHSType> m_rhs;
  };

  template <class LHS, class Number>
  Mult(LHS, const LocalBilinearFormIntegratorBase<Number>&)
    -> Mult<LHS, LocalBilinearFormIntegratorBase<Number>>;

  template <class LHS, class Number>
  inline
  constexpr
  auto
  operator*(LHS lhs, const LocalBilinearFormIntegratorBase<Number>& rhs)
  {
    return Mult(lhs, rhs);
  }

  template <class LHS, class Number>
  class Mult<LHS, LinearFormIntegratorBase<Number>>
    : public LinearFormIntegratorBase<decltype(std::declval<LHS>() * std::declval<Number>())>
  {
    public:
      using LHSType = LHS;

      using RHSType = LinearFormIntegratorBase<Number>;

      using ScalarType = decltype(std::declval<LHS>() * std::declval<Number>());

      using Parent = LinearFormIntegratorBase<ScalarType>;

      Mult(ScalarType lhs, const RHSType& rhs)
        : Parent(rhs),
          m_lhs(lhs), m_rhs(rhs.copy())
      {}

      Mult(const Mult& other)
        : Parent(other),
          m_lhs(other.m_lhs), m_rhs(other.m_rhs->copy())
      {}

      Mult(Mult&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
      {}

      Integrator::Region getRegion() const override
      {
        return m_rhs->getRegion();
      }

      const LHSType& getLHS() const
      {
        return m_lhs;
      }

      const RHSType& getRHS() const
      {
        assert(m_rhs);
        return *m_rhs;
      }

      const Geometry::Polytope& getPolytope() const override
      {
        return m_rhs->getPolytope();
      }

      Mult& setPolytope(const Geometry::Polytope& polytope) override
      {
        m_rhs->setPolytope(polytope);
        return *this;
      }

      ScalarType integrate(size_t local) override
      {
        return getLHS() * m_rhs->integrate(local);
      }

      Mult* copy() const noexcept override
      {
        return new Mult(*this);
      }

    private:
      const LHSType m_lhs;
      std::unique_ptr<RHSType> m_rhs;
  };

  template <class LHS, class Number>
  Mult(LHS, const LinearFormIntegratorBase<Number>&)
    -> Mult<LHS, LinearFormIntegratorBase<Number>>;

  template <class Number>
  inline
  constexpr
  auto
  operator*(Number lhs, const LinearFormIntegratorBase<Number>& rhs)
  {
    return Mult(lhs, rhs);
  }
}

#endif
