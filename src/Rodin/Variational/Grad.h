/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_GRADIENT_H
#define RODIN_VARIATIONAL_GRADIENT_H

#include "ForwardDecls.h"

#include "GridFunction.h"
#include "TestFunction.h"
#include "TrialFunction.h"
#include "VectorFunction.h"

#include "Rodin/Variational/Exceptions/UndeterminedTraceDomainException.h"

namespace Rodin::FormLanguage
{
  template <class FES, class Data>
  struct Traits<Variational::Grad<Variational::GridFunction<FES, Data>>>
  {
    using FESType = FES;

    using OperandType = Variational::GridFunction<FESType, Data>;

    using RangeType = Math::Vector<typename FormLanguage::Traits<FESType>::ScalarType>;
  };

  template <class NestedDerived, class FES, Variational::ShapeFunctionSpaceType Space>
  struct Traits<
    Variational::Grad<Variational::ShapeFunction<NestedDerived, FES, Space>>>
  {
    using FESType = FES;
    static constexpr Variational::ShapeFunctionSpaceType SpaceType = Space;

    using OperandType = Variational::ShapeFunction<NestedDerived, FESType, SpaceType>;

    using RangeType = Math::Vector<typename FormLanguage::Traits<FESType>::ScalarType>;
  };
}

namespace Rodin::Variational
{
  /**
   * @defgroup GradSpecializations Grad Template Specializations
   * @brief Template specializations of the Grad class.
   * @see Grad
   */

  /**
   * @brief Base class for Grad classes.
   */
  template <class Operand, class Derived>
  class GradBase;

  /**
   * @ingroup GradSpecializations
   */
  template <class FES, class Data, class Derived>
  class GradBase<GridFunction<FES, Data>, Derived>
    : public VectorFunctionBase<
        typename FormLanguage::Traits<FES>::ScalarType, GradBase<GridFunction<FES, Data>, Derived>>
  {
    public:
      using FESType = FES;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using SpatialVectorType = Math::SpatialVector<ScalarType>;

      using OperandType = GridFunction<FESType, Data>;

      using Parent = VectorFunctionBase<ScalarType, GradBase<OperandType, Derived>>;

      GradBase(const OperandType& u)
        : m_u(u)
      {
        assert(u.getFiniteElementSpace().getVectorDimension() == 1);
      }

      /**
       * @brief Copy constructor
       */
      GradBase(const GradBase& other)
        : Parent(other),
          m_u(other.m_u)
      {}

      /**
       * @brief Move constructor
       */
      GradBase(GradBase&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u))
      {}

      constexpr
      size_t getDimension() const
      {
        return m_u.get().getFiniteElementSpace().getMesh().getSpaceDimension();
      }

      SpatialVectorType getValue(const Geometry::Point& p) const
      {
        SpatialVectorType out;
        this->getValue(out, p);
        return out;
      }

      void getValue(SpatialVectorType& out, const Geometry::Point& p) const
      {
        const auto& polytope = p.getPolytope();
        const auto& polytopeMesh = polytope.getMesh();
        const auto& gf = getOperand();
        const auto& fes = gf.getFiniteElementSpace();
        const auto& fesMesh = fes.getMesh();
        if (polytopeMesh == fesMesh)
        {
          this->interpolate(out, p);
        }
        else if (const auto inclusion = fesMesh.inclusion(p))
        {
          this->interpolate(out, *inclusion);
        }
        else if (fesMesh.isSubMesh())
        {
          const auto& submesh = fesMesh.asSubMesh();
          const auto restriction = submesh.restriction(p);
          this->interpolate(out, *restriction);
        }
        else
        {
          assert(false);
        }
      }

      /**
       * @brief Interpolation function to be overriden in Derived type.
       */
      constexpr
      void interpolate(SpatialVectorType& out, const Geometry::Point& p) const
      {
        static_cast<const Derived&>(*this).interpolate(out, p);
      }

      constexpr
      const OperandType& getOperand() const
      {
        return m_u.get();
      }

      /**
       * @brief Copy function to be overriden in Derived type.
       */
      GradBase* copy() const noexcept override
      {
        return static_cast<const Derived&>(*this).copy();
      }

    private:
      std::reference_wrapper<const OperandType> m_u;
  };

  template <class FES, class Data>
  class Grad<GridFunction<FES, Data>> final
    : public GradBase<GridFunction<FES, Data>, Grad<GridFunction<FES, Data>>>
  {
    public:
      using FESType = FES;

      using RangeType = typename FormLanguage::Traits<FESType>::RangeType;

      using ScalarType = typename FormLanguage::Traits<RangeType>::ScalarType;

      using SpatialVectorType = Math::SpatialVector<ScalarType>;

      using OperandType = GridFunction<FESType, Data>;

      using Parent = GradBase<OperandType, Grad<OperandType>>;

      Grad(const OperandType& u)
        : Parent(u)
      {}

      /**
       * @brief Copy constructor
       */
      Grad(const Grad& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor
       */
      Grad(Grad&& other)
        : Parent(std::move(other))
      {}

      void interpolate(SpatialVectorType& out, const Geometry::Point& p) const
      {
        const auto& polytope = p.getPolytope();
        const auto& d = polytope.getDimension();
        const auto& i = polytope.getIndex();
        const auto& mesh = polytope.getMesh();
        const size_t meshDim = mesh.getDimension();
        if (d == meshDim - 1) // Evaluating on a face
        {
          const auto& conn = mesh.getConnectivity();
          const auto& inc = conn.getIncidence({ meshDim - 1, meshDim }, i);
          const auto& pc = p.getPhysicalCoordinates();
          assert(inc.size() == 1 || inc.size() == 2);
          if (inc.size() == 1)
          {
            const auto& tracePolytope = mesh.getPolytope(meshDim, *inc.begin());
            const auto rc = tracePolytope->getTransformation().inverse(pc);
            const Geometry::Point np(*tracePolytope, std::cref(rc), pc);
            this->interpolate(out, np);
            return;
          }
          else
          {
            assert(inc.size() == 2);
            const auto& traceDomain = this->getTraceDomain();
            assert(traceDomain.size() > 0);
            if (traceDomain.size() == 0)
            {
              Alert::MemberFunctionException(*this, __func__)
                << "No trace domain provided: "
                << Alert::Notation::Predicate(true, "getTraceDomain().size() == 0")
                << ". Grad at an interface with no trace domain is undefined."
                << Alert::Raise;
            }
            else
            {
              for (auto& idx : inc)
              {
                const auto& tracePolytope = mesh.getPolytope(meshDim, idx);
                if (traceDomain.count(tracePolytope->getAttribute()))
                {
                  const auto rc = tracePolytope->getTransformation().inverse(pc);
                  const Geometry::Point np(*tracePolytope, std::cref(rc), pc);
                  this->interpolate(out, np);
                  return;
                }
              }
              UndeterminedTraceDomainException(
                  *this, __func__, {d, i}, traceDomain.begin(), traceDomain.end()) << Alert::Raise;
            }
            return;
          }
        }
        else // Evaluating on a cell
        {
          static thread_local SpatialVectorType s_gradient;
          static thread_local SpatialVectorType s_res;

          s_res.resize(d);
          s_res.setZero();

          assert(d == mesh.getDimension());
          const auto& gf = this->getOperand();
          const auto& fes = gf.getFiniteElementSpace();
          const auto& fe = fes.getFiniteElement(d, i);
          const auto& rc = p.getReferenceCoordinates();
          for (size_t local = 0; local < fe.getCount(); local++)
          {
            const auto basis = fe.getBasis(local);
            basis.getGradient()(s_gradient, rc);
            s_res += gf[fes.getGlobalIndex({d, i}, local)] * s_gradient;
          }
          out = p.getJacobianInverse().transpose() * s_res;
        }
      }

      Grad* copy() const noexcept override
      {
        return new Grad(*this);
      }
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for Grad of a GridFunction
   */
  template <class FES, class Data>
  Grad(const GridFunction<FES, Data>&) -> Grad<GridFunction<FES, Data>>;

  /**
   * @ingroup GradSpecializations
   * @brief Gradient of a ShapeFunction
   */
  template <class NestedDerived, class FES, ShapeFunctionSpaceType SpaceType>
  class Grad<ShapeFunction<NestedDerived, FES, SpaceType>> final
    : public ShapeFunctionBase<Grad<ShapeFunction<NestedDerived, FES, SpaceType>>>
  {
    public:
      /// Finite element space type
      using FESType = FES;
      static constexpr ShapeFunctionSpaceType Space = SpaceType;

      /// Type of scalar values in the finite element space
      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using RangeType = Math::Vector<ScalarType>;

      /// Operand type
      using OperandType = ShapeFunction<NestedDerived, FESType, Space>;

      /// Parent class
      using Parent = ShapeFunctionBase<Grad<OperandType>, FESType, Space>;

      Grad(const OperandType& u)
        : Parent(u.getFiniteElementSpace()),
          m_u(u),
          m_p(nullptr)
      {}

      Grad(const Grad& other)
        : Parent(other),
          m_u(other.m_u),
          m_p(other.m_p),
          m_gradient(other.m_gradient)
      {}

      Grad(Grad&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)),
          m_p(std::exchange(other.m_p, nullptr)),
          m_gradient(std::move(other.m_gradient))
      {}

      constexpr
      const OperandType& getOperand() const
      {
        return m_u.get();
      }

      constexpr
      const auto& getLeaf() const
      {
        return getOperand().getLeaf();
      }

      constexpr
      size_t getDOFs(const Geometry::Polytope& element) const
      {
        return getOperand().getDOFs(element);
      }

      const Geometry::Point& getPoint() const
      {
        assert(m_p);
        return *m_p;
      }

      Grad& setPoint(const Geometry::Point& p)
      {
        if (m_p == &p)
          return *this;
        m_p = &p;
        const auto& polytope = p.getPolytope();
        const auto& rc = p.getReferenceCoordinates();
        const size_t d = polytope.getDimension();
        const Index i = polytope.getIndex();
        const auto& fes = this->getFiniteElementSpace();
        const auto& fe = fes.getFiniteElement(d, i);
        const size_t count = fe.getCount();
        m_gradient.resize(count);
        for (size_t local = 0; local < count; local++)
        {
          const auto basis = fe.getBasis(local);
          basis.getGradient()(m_gradient[local], rc);
        }
        return *this;
      }

      constexpr
      auto getBasis(size_t local) const
      {
        return getPoint().getJacobianInverse().transpose() * m_gradient[local];
      }

      Grad* copy() const noexcept override
      {
        return new Grad(*this);
      }

    private:
      std::reference_wrapper<const OperandType> m_u;

      const Geometry::Point* m_p;
      std::vector<Math::SpatialVector<ScalarType>> m_gradient;
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for Grad of a ShapeFunction
   */
  template <class NestedDerived, class FES, ShapeFunctionSpaceType Space>
  Grad(const ShapeFunction<NestedDerived, FES, Space>&)
    -> Grad<ShapeFunction<NestedDerived, FES, Space>>;

}

#endif
