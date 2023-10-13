/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_POTENTIAL_H
#define RODIN_VARIATIONAL_POTENTIAL_H

#include "Rodin/FormLanguage/Base.h"
#include "Rodin/FormLanguage/List.h"
#include "Rodin/QF/QuadratureFormula.h"

#include "ForwardDecls.h"
#include "Mult.h"
#include "Function.h"
#include "Integral.h"
#include "ShapeFunction.h"

namespace Rodin::FormLanguage
{
  /**
   * @ingroup TraitsSpecializations
   * @brief Traits for Potential.
   */
  template <class LHSType, class RHSDerived>
  struct Traits<
    Variational::Potential<
      LHSType,
      Variational::FunctionBase<RHSDerived>>>
  {
    using LHS = LHSType;

    using RHS = Variational::FunctionBase<Variational::FunctionBase<RHSDerived>>;

    using Kernel = LHS;

    using Operand = RHS;

    using RHSRange =
      typename FormLanguage::Traits<RHS>::RangeType;

    using LHSRange =
      std::conditional_t<
        // If
        std::is_same_v<RHSRange, Scalar>,
        // Then
        Scalar,
        // Else
        std::conditional_t<
          // If
          std::is_same_v<RHSRange, Math::Vector>,
          // Then
          Math::Matrix,
          // Else
          void>>;

    using RangeType = RHSRange;
  };

  /**
   * @ingroup TraitsSpecializations
   * @brief Traits for Potential.
   */
  template <class LHSType, class RHSDerived, class FESType, Variational::ShapeFunctionSpaceType SpaceType>
  struct Traits<
    Variational::Potential<
      LHSType,
      Variational::ShapeFunctionBase<Variational::ShapeFunction<RHSDerived, FESType, SpaceType>>>>
  {
    using FES = FESType;
    static constexpr Variational::ShapeFunctionSpaceType Space = SpaceType;

    using LHS = LHSType;

    using RHS = Variational::ShapeFunctionBase<Variational::ShapeFunction<RHSDerived, FES, Space>, FES, Space>;

    using Kernel = LHS;

    using Operand = RHS;

    using RHSRange =
      typename FormLanguage::Traits<RHS>::RangeType;

    using LHSRange =
      std::conditional_t<
        // If
        std::is_same_v<RHSRange, Scalar>,
        // Then
        Scalar,
        // Else
        std::conditional_t<
          // If
          std::is_same_v<RHSRange, Math::Vector>,
          // Then
          Math::Matrix,
          // Else
          void>>;

    using RangeType = RHSRange;
  };
}

namespace Rodin::Variational
{
  /**
   * @defgroup PotentialSpecializations Potential Template Specializations
   * @brief Template specializations of the Potential class.
   * @see Potential
   */

  /**
   * @ingroup PotentialSpecializations
   */
  template <class LHSType, class RHSDerived>
  class Potential<LHSType, FunctionBase<RHSDerived>> final
    : public FunctionBase<Potential<LHSType, FunctionBase<RHSDerived>>>
  {
    public:
      using LHS = LHSType;

      using Kernel = LHS;

      using RHS = FunctionBase<RHSDerived>;

      using Operand = RHS;

      using RHSRange = typename FormLanguage::Traits<RHS>::RangeType;

      using LHSRange =
        std::conditional_t<
          // If
          std::is_same_v<RHSRange, Scalar>,
          // Then
          Scalar,
          // Else
          std::conditional_t<
            // If
            std::is_same_v<RHSRange, Math::Vector>,
            // Then
            Math::Matrix,
            // Else
            void>>;

      using Parent = FunctionBase<Potential<LHS, RHS>>;

      Potential(const Kernel& kernel, const Operand& u)
        : m_kernel(kernel), m_u(u.copy())
      {}

      Potential(const Potential& other)
        : Parent(other),
          m_kernel(other.m_kernel),
          m_u(other.m_u->copy())
      {}

      Potential(Potential&& other)
        : Parent(std::move(other)),
          m_kernel(std::move(other.m_kernel)),
          m_u(std::move(other.m_u))
      {}

      inline
      const auto& getKernel() const
      {
        return m_kernel.get();
      }

      inline
      const auto& getOperand() const
      {
        assert(m_u);
        return *m_u;
      }

      inline
      auto getValue(const Geometry::Point& p) const
      {
        const auto& kernel = getKernel();
        const auto& operand = getOperand();
        const auto& mesh = p.getPolytope().getMesh();
        if (m_qf.has_value())
        {
          if constexpr (std::is_same_v<RHSRange, Scalar>)
          {
            Scalar res = 0;
            for (auto it = mesh.getElement(); it; ++it)
            {
              const auto& polytope = *it;
              const auto& qf = m_qf.value()(polytope);
              const auto& trans = polytope.getTransformation();
              for (size_t i = 0; i < qf.getSize(); i++)
              {
                Geometry::Point y(polytope, trans, std::ref(qf.getPoint(i)));
                res += qf.getWeight(i) * p.getDistortion() * kernel(p, y) * operand(y);
              }
            }
            return res;
          }
          else if constexpr (std::is_same_v<RHSRange, Math::Vector>)
          {
            Math::Vector res;
            getValue(res, p);
            return res;
          }
          else
          {
            assert(false);
            return void();
          }
        }
        else
        {
          if constexpr (std::is_same_v<RHSRange, Scalar>)
          {
            Scalar res = 0;
            for (auto it = mesh.getElement(); it; ++it)
            {
              const auto& polytope = *it;
              const QF::QFGG qf(polytope.getGeometry());
              const auto& trans = polytope.getTransformation();
              for (size_t i = 0; i < qf.getSize(); i++)
              {
                Geometry::Point y(polytope, trans, std::ref(qf.getPoint(i)));
                res += qf.getWeight(i) * p.getDistortion() * kernel(p, y) * operand(y);
              }
            }
            return res;
          }
          else if constexpr (std::is_same_v<RHSRange, Math::Vector>)
          {
            Math::Vector res;
            getValue(res, p);
            return res;
          }
          else
          {
            assert(false);
            return void();
          }
        }
      }

      inline
      constexpr
      void getValue(Math::Vector& res, const Geometry::Point& p) const
      {
        assert(false);
        res.setConstant(NAN);
      }

      Potential& setQuadratureFormula(
          const std::function<const QF::QuadratureFormulaBase&(const Geometry::Polytope&)>& qf)
      {
        m_qf.emplace(qf);
        return *this;
      }

    private:
      std::reference_wrapper<const Kernel> m_kernel;
      std::unique_ptr<Operand> m_u;
      std::optional<
        std::function<const QF::QuadratureFormulaBase&(const Geometry::Polytope&)>> m_qf;
  };

  /**
   * @brief CTAD for Potential.
   */
  template <class LHSType, class RHSDerived>
  Potential(const LHSType&, const FunctionBase<RHSDerived>&)
    -> Potential<LHSType, FunctionBase<RHSDerived>>;

  /**
   * @ingroup PotentialSpecializations
   */
  template <class LHSType, class RHSDerived, class FESType, ShapeFunctionSpaceType SpaceType>
  class Potential<
    LHSType,
    ShapeFunctionBase<ShapeFunction<RHSDerived, FESType, SpaceType>, FESType, SpaceType>> final
    : public
      ShapeFunctionBase<
        Potential<
          LHSType,
          ShapeFunctionBase<ShapeFunction<RHSDerived, FESType, SpaceType>, FESType, SpaceType>>, FESType, SpaceType>
  {
    public:
      using FES = FESType;
      static constexpr ShapeFunctionSpaceType Space = SpaceType;

      using LHS = LHSType;

      using Kernel = LHS;

      using RHS = ShapeFunctionBase<ShapeFunction<RHSDerived, FES, SpaceType>, FES, Space>;

      using Operand = RHS;

      using Parent = ShapeFunctionBase<Potential<LHS, RHS>, FES, Space>;

      using RHSRange = typename FormLanguage::Traits<RHS>::RangeType;

      using LHSRange =
        std::conditional_t<
          // If
          std::is_same_v<RHSRange, Scalar>,
          // Then
          Scalar,
          // Else
          std::conditional_t<
            // If
            std::is_same_v<RHSRange, Math::Vector>,
            // Then
            Math::Matrix,
            // Else
            void>>;

      static_assert(
          (std::is_same_v<LHSRange, Scalar> && std::is_same_v<RHSRange, Scalar>) ||
          (std::is_same_v<LHSRange, Math::Matrix> || std::is_same_v<RHSRange, Math::Vector>));

      Potential(const Kernel& kernel, const Operand& u)
        : Parent(u.getFiniteElementSpace()),
          m_kernel(kernel), m_u(u)
      {}

      Potential(const Potential& other)
        : Parent(other),
          m_kernel(other.m_kernel), m_u(other.m_u)
      {}

      Potential(Potential&& other)
        : Parent(std::move(other)),
          m_kernel(std::move(other.m_kernel)), m_u(std::move(other.m_u))
      {}

      inline
      const LHS& getLHS() const
      {
        return getKernel();
      }

      inline
      const RHS& getRHS() const
      {
        return getOperand();
      }

      inline
      const Kernel& getKernel() const
      {
        assert(m_kernel);
        return *m_kernel;
      }

      inline
      const Operand& getOperand() const
      {
        return m_u.get();
      }

      inline
      constexpr
      const auto& getLeaf() const
      {
        return getOperand();
      }

      inline
      constexpr
      RangeShape getRangeShape() const
      {
        if constexpr (std::is_same_v<LHSRange, Scalar>)
        {
          static_assert(std::is_same_v<RHSRange, Scalar>);
          return { 1, 1 };
        }
        else if constexpr (std::is_same_v<LHSRange, Math::Matrix>)
        {
          static_assert(std::is_same_v<RHSRange, Math::Vector>);
          return getRHS().getRangeShape()[0];
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
      auto getTensorBasis(const Geometry::Point& p) const
      {
        using RangeType = typename FES::RangeType;
        const size_t d = p.getPolytope().getDimension();
        const Index i = p.getPolytope().getIndex();
        const auto& fes = this->getFiniteElementSpace();
        const auto& fe = fes.getFiniteElement(d, i);
        const auto& kernel = getKernel();
        const auto& mesh = p.getPolytope().getMesh();
        if (m_qf.has_value())
        {
          assert(false);
        }
        else
        {
          if constexpr (std::is_same_v<RangeType, Scalar>)
          {
            Scalar res = 0;
            for (auto it = mesh.getElement(); it; ++it)
            {
              const auto& polytope = *it;
              const QF::QFGG qf(polytope.getGeometry());
              const auto& trans = polytope.getTransformation();
              for (size_t i = 0; i < qf.getSize(); i++)
              {
                for (size_t local = 0; local < fe.getCount(); local++)
                {
                  Geometry::Point y(polytope, trans, std::ref(qf.getPoint(i)));
                  const auto psi = fes.getInverseMapping(polytope, fe.getBasis(local));
                  res += qf.getWeight(i) * p.getDistortion() * kernel(p, y) * psi(y);
                }
              }
            }
            return res;
          }
          else if constexpr (std::is_same_v<RangeType, Math::Vector>)
          {
            assert(false);
            return void();
          }
          else
          {
            assert(false);
            return void();
          }
        }
      }

      inline
      constexpr
      const auto& getFiniteElementSpace() const
      {
        return getOperand().getFiniteElementSpace();
      }

      Potential& setQuadratureFormula(
          const std::function<const QF::QuadratureFormulaBase&(const Geometry::Polytope&)>& qf)
      {
        m_qf.emplace(qf);
        return *this;
      }

      inline Potential* copy() const noexcept override
      {
        return new Potential(*this);
      }

    private:
      std::reference_wrapper<Kernel> m_kernel;
      std::reference_wrapper<const Operand> m_u;
      std::optional<
        std::function<const QF::QuadratureFormulaBase&(const Geometry::Polytope&)>> m_qf;
  };

  /**
   * @brief CTAD for Potential.
   */
  template <class LHSType, class RHSDerived, class FESType, ShapeFunctionSpaceType SpaceType>
  Potential(const LHSType&, const ShapeFunctionBase<ShapeFunction<RHSDerived, FESType, SpaceType>, FESType, SpaceType>&)
    -> Potential<LHSType, ShapeFunctionBase<ShapeFunction<RHSDerived, FESType, SpaceType>, FESType, SpaceType>>;
}

#endif

