/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_H1_GRAD_H
#define RODIN_VARIATIONAL_H1_GRAD_H

/**
 * @file
 * @brief Gradient operator specialization for H1 (higher-order Lagrange) functions.
 *
 * For H1<K> functions, the gradient is polynomial of degree K-1 on each element:
 * @f[
 *   \nabla u|_K = \sum_{i=1}^{n} u_i \nabla \phi_i
 * @f]
 * where @f$ \phi_i @f$ are the H1<K> basis functions.
 */

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/Point.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Variational/Grad.h"
#include "Rodin/Variational/IntegrationPoint.h"
#include "Rodin/Variational/ShapeFunction.h"
#include "Rodin/Variational/Exceptions/UndeterminedTraceDomainException.h"

#include "ForwardDecls.h"

namespace Rodin::Variational
{
  /**
   * @defgroup GradSpecializations Grad Template Specializations
   * @brief Template specializations of the Grad class.
   * @see Grad
   */

  /**
   * @ingroup GradSpecializations
   * @brief Gradient of a GridFunction on H1<K> space
   */
  template <size_t K, class Scalar, class Mesh, class Data>
  class Grad<GridFunction<H1<K, Scalar, Mesh>, Data>> final
    : public GradBase<
        GridFunction<H1<K, Scalar, Mesh>, Data>,
        Grad<GridFunction<H1<K, Scalar, Mesh>, Data>>>
  {
    public:
      using FESType = H1<K, Scalar, Mesh>;
      using RangeType = typename FormLanguage::Traits<FESType>::RangeType;
      using ScalarType = typename FormLanguage::Traits<RangeType>::ScalarType;
      using SpatialVectorType = Math::SpatialVector<ScalarType>;
      using OperandType = GridFunction<FESType, Data>;
      using Parent = GradBase<OperandType, Grad<OperandType>>;

      Grad(const OperandType& u) : Parent(u) {}
      Grad(const Grad& other) : Parent(other) {}
      Grad(Grad&& other) : Parent(std::move(other)) {}

      void interpolate(SpatialVectorType& out, const IntegrationPoint& ip) const
      {
        const auto& p = ip.getPoint();
        const auto& polytope = p.getPolytope();
        const size_t d = polytope.getDimension();
        const Index  i = polytope.getIndex();

        const auto& gf  = this->getOperand();
        const auto& fes = gf.getFiniteElementSpace();
        const auto& fe  = fes.getFiniteElement(d, i);
        const auto* qf = ip.getQuadratureFormula();
        assert(qf);
        const auto& tab = fe.getTabulation(*qf);
        const auto JinvT = p.getJacobianInverse().transpose();

        SpatialVectorType ref(static_cast<std::uint8_t>(d));
        ref.setZero();

        for (size_t local = 0; local < fe.getCount(); ++local)
        {
          const auto gref = tab.getGradient(ip.getIndex(), local);
          const auto uval = gf[fes.getGlobalIndex({d, i}, local)];
          for (size_t j = 0; j < d; ++j)
            ref(static_cast<std::uint8_t>(j)) += uval * gref[j];
        }

        out = JinvT * ref;
      }

      void interpolate(SpatialVectorType& out, const Geometry::Point& p) const
      {
        const auto& polytope = p.getPolytope();
        const auto& d = polytope.getDimension();
        const auto& i = polytope.getIndex();
        const auto& mesh = polytope.getMesh();
        const size_t meshDim = mesh.getDimension();

        if (d == meshDim - 1) // face
        {
          const auto& conn = mesh.getConnectivity();
          const auto& inc = conn.getIncidence({ meshDim - 1, meshDim }, i);
          const auto& pc = p.getPhysicalCoordinates();
          assert(inc.size() == 1 || inc.size() == 2);

          if (inc.size() == 1)
          {
            const auto& tracePolytope = mesh.getPolytope(meshDim, *inc.begin());
            Math::SpatialPoint rc;
            tracePolytope->getTransformation().inverse(rc, pc);
            const Geometry::Point np(*tracePolytope, std::cref(rc), pc);
            this->interpolate(out, np);
            return;
          }

          assert(inc.size() == 2);
          const auto& traceDomain = this->getTraceDomain();
          if (traceDomain.size() == 0)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "No trace domain provided: "
              << Alert::Notation::Predicate(true, "getTraceDomain().size() == 0")
              << ". Grad at an interface with no trace domain is undefined."
              << Alert::Raise;
          }

          for (auto& idx : inc)
          {
            const auto& tracePolytope = mesh.getPolytope(meshDim, idx);
            const auto a = tracePolytope->getAttribute();
            if (a && traceDomain.count(*a))
            {
              Math::SpatialPoint rc;
              tracePolytope->getTransformation().inverse(rc, pc);
              const Geometry::Point np(*tracePolytope, std::cref(rc), pc);
              this->interpolate(out, np);
              return;
            }
          }

          UndeterminedTraceDomainException(
              *this, __func__, {d, i}, traceDomain.begin(), traceDomain.end())
            << Alert::Raise;
          return;
        }
        else // cell
        {
          out.resize(d);
          out.setZero();

          assert(d == mesh.getDimension());

          const auto& gf  = this->getOperand();
          const auto& fes = gf.getFiniteElementSpace();
          const auto& fe  = fes.getFiniteElement(d, i);
          const auto& rc  = p.getReferenceCoordinates();

          for (size_t local = 0; local < fe.getCount(); ++local)
          {
            const auto& basis = fe.getBasis(local);
            out += gf[fes.getGlobalIndex({d, i}, local)] * basis.getGradient()(rc);
          }

          out = p.getJacobianInverse().transpose() * out;
        }
      }

      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope& geom) const noexcept
      {
        const size_t k = H1Element<K, ScalarType>(geom.getGeometry()).getOrder();
        return (k == 0) ? 0 : (k - 1);
      }

      Grad* copy() const noexcept override
      {
        return new Grad(*this);
      }
  };

  /**
   * @ingroup GradSpecializations
   * @brief Gradient of a ShapeFunction on H1<K> space
   */
  template <size_t K, class NestedDerived, class Scalar, class Mesh, ShapeFunctionSpaceType SpaceType>
  class Grad<ShapeFunction<NestedDerived, H1<K, Scalar, Mesh>, SpaceType>> final
    : public ShapeFunctionBase<Grad<ShapeFunction<NestedDerived, H1<K, Scalar, Mesh>, SpaceType>>>
  {
    public:
      using FESType = H1<K, Scalar, Mesh>;
      static constexpr ShapeFunctionSpaceType Space = SpaceType;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;
      using RangeType  = Math::SpatialVector<ScalarType>;
      using OperandType = ShapeFunction<NestedDerived, FESType, Space>;
      using Parent = ShapeFunctionBase<Grad<OperandType>, FESType, Space>;

      using SpatialVectorType = Math::SpatialVector<ScalarType>;

      struct Cache
      {
        struct Key
        {
          Geometry::Polytope::Type geom = Geometry::Polytope::Type::Point;
          size_t dim = 0;
          Index cell = 0;

          const QF::QuadratureFormulaBase* qf = nullptr;
          size_t qp = 0;

          bool valid = false;

          explicit operator bool() const noexcept { return valid; }

          bool operator==(const Key& o) const noexcept
          {
            if (!valid || !o.valid)
              return false;
            return geom == o.geom
                && dim  == o.dim
                && cell == o.cell
                && qf   == o.qf
                && qp   == o.qp;
          }

          void operator=(std::initializer_list<int>) noexcept
          {
            valid = false;
            geom = Geometry::Polytope::Type::Point;
            dim = 0;
            cell = 0;
            qf = nullptr;
            qp = 0;
          }
        };

        std::vector<SpatialVectorType> grad_phys; // size = ndof, each dim = d
        Key key;
      };

      Grad(const OperandType& u)
        : Parent(u.getFiniteElementSpace()),
          m_u(u),
          m_ip(nullptr)
      {}

      Grad(const Grad& other)
        : Parent(other),
          m_u(other.m_u),
          m_ip(nullptr),
          m_cache(other.m_cache)
      {}

      Grad(Grad&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)),
          m_ip(std::exchange(other.m_ip, nullptr)),
          m_cache(std::move(other.m_cache))
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

      constexpr
      const IntegrationPoint& getIntegrationPoint() const
      {
        assert(m_ip);
        return *m_ip;
      }

      Grad& setIntegrationPoint(const IntegrationPoint& ip)
      {
        m_ip = &ip;

        const auto& p  = ip.getPoint();
        const auto* qf = ip.getQuadratureFormula();
        const size_t qp = qf ? ip.getIndex() : 0;

        const auto& poly = p.getPolytope();
        const size_t d   = poly.getDimension();
        const Index  cell = poly.getIndex();
        const auto geom = poly.getGeometry();

        typename Cache::Key key;
        key.geom  = geom;
        key.dim   = d;
        key.cell  = cell;
        key.qf    = qf;
        key.qp    = qp;
        key.valid = true;

        const bool recompute = !qf || !(m_cache.key == key);
        if (!recompute)
          return *this;

        m_cache.key = key;

        const auto& fes = this->getFiniteElementSpace();
        const auto& fe  = fes.getFiniteElement(d, cell);
        const size_t ndof = fe.getCount();

        // Ensure storage sized once.
        if (m_cache.grad_phys.size() != ndof)
          m_cache.grad_phys.resize(ndof);

        for (auto& g : m_cache.grad_phys)
        {
          if (g.size() != d)
            g.resize(d);
        }

        // Reference gradients from tabulation when integrating, otherwise
        // directly from the basis at the supplied point.
        const auto* tab = qf ? &fe.getTabulation(*qf) : nullptr;
        const auto& rc = p.getReferenceCoordinates();
        const auto JinvT = p.getJacobianInverse().transpose();

        SpatialVectorType ref(d);

        for (size_t a = 0; a < ndof; ++a)
        {
          for (size_t ii = 0; ii < d; ++ii)
            ref(ii) =
              qf
                ? tab->getGradient(qp, a)[ii]
                : fe.getBasis(a).template getDerivative<1>(ii)(rc);

          m_cache.grad_phys[a] = JinvT * ref;
        }

        return *this;
      }

      RangeType getBasis(size_t local) const
      {
        assert(m_cache.key);
        assert(local < m_cache.grad_phys.size());
        return m_cache.grad_phys[local];
      }

      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope& geom) const noexcept
      {
        const auto k = getOperand().getOrder(geom);
        if (!k.has_value())
          return std::nullopt;
        return (*k == 0) ? 0 : (*k - 1);
      }

      Grad* copy() const noexcept override
      {
        return new Grad(*this);
      }

    private:
      std::reference_wrapper<const OperandType> m_u;

      const IntegrationPoint* m_ip;
      Cache m_cache;
  };
}

#endif
