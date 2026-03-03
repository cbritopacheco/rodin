/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_H1_DIV_H
#define RODIN_VARIATIONAL_H1_DIV_H

/**
 * @file
 * @brief Divergence operator specialization for H1 vector functions.
 *
 * This file provides specialized implementations of the divergence operator
 * for H1 vector-valued GridFunctions and ShapeFunctions.
 *
 * For H1<K> elements, the divergence is computed as:
 * @f[
 *   \nabla \cdot \mathbf{u}|_K = \sum_{i=1}^{n_{dof}} \sum_{j=1}^d u_{i,j} \frac{\partial \phi_i}{\partial x_j}
 * @f]
 * where @f$ \phi_i @f$ are H1<K> basis functions and @f$ u_{i,j} @f$ are DOF values.
 *
 * For H1<K> functions, the divergence is polynomial of degree K-1 on each element.
 */

#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Variational/Div.h"
#include "Rodin/Variational/IntegrationPoint.h"

namespace Rodin::FormLanguage
{
  template <size_t K, class Scalar, class Data, class Mesh>
  struct Traits<Variational::Div<Variational::GridFunction<Variational::H1<K, Math::Vector<Scalar>, Mesh>, Data>>>
  {
    using FESType = Variational::H1<K, Math::Vector<Scalar>, Mesh>;
    using ScalarType = Scalar;
    using OperandType = Variational::GridFunction<Variational::H1<K, Math::Vector<Scalar>, Mesh>, Data>;
  };

  template <size_t K, class NestedDerived, class Scalar, class Mesh, Variational::ShapeFunctionSpaceType Space>
  struct Traits<
    Variational::Div<
      Variational::ShapeFunction<NestedDerived, Variational::H1<K, Math::Vector<Scalar>, Mesh>, Space>>>
  {
    using FESType = Variational::H1<K, Math::Vector<Scalar>, Mesh>;
    static constexpr Variational::ShapeFunctionSpaceType SpaceType = Space;
    using ScalarType = Scalar;
    using OperandType =
      Variational::ShapeFunction<NestedDerived, Variational::H1<K, Math::Vector<Scalar>, Mesh>, Space>;
  };
}

namespace Rodin::Variational
{
  /**
   * @ingroup DivSpecializations
   * @brief Divergence of an H1 vector GridFunction.
   *
   * Computes @f$ \nabla \cdot \mathbf{u} @f$ for vector-valued H1<K> functions.
   * The divergence is polynomial of degree K-1 on each element:
   * @f[
   *   (\nabla \cdot \mathbf{u})|_K = \sum_{i=1}^{n_{dof}} \mathbf{u}_i \cdot \nabla \phi_i
   * @f]
   *
   * Applications include:
   * - Incompressibility constraint: @f$ \nabla \cdot \mathbf{v} = 0 @f$
   * - Mass conservation in fluid mechanics
   * - Continuity equation
   * - Mixed formulations (Stokes, elasticity)
   */
  template <size_t K, class Scalar, class Data, class Mesh>
  class Div<GridFunction<H1<K, Math::Vector<Scalar>, Mesh>, Data>> final
    : public DivBase<
        GridFunction<H1<K, Math::Vector<Scalar>, Mesh>, Data>,
        Div<GridFunction<H1<K, Math::Vector<Scalar>, Mesh>, Data>>>
  {
    public:
      using FESType = H1<K, Math::Vector<Scalar>, Mesh>;
      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using SpatialVectorType = Math::SpatialVector<ScalarType>;

      using OperandType = GridFunction<FESType, Data>;
      using Parent = DivBase<OperandType, Div<OperandType>>;

      Div(const OperandType& u) : Parent(u) {}
      Div(const Div& other) : Parent(other) {}
      Div(Div&& other) : Parent(std::move(other)) {}

      void interpolate(ScalarType& out, const Geometry::Point& p) const
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
              << ". Div at an interface with no trace domain is undefined."
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
          assert(d == mesh.getDimension());

          const auto& gf  = this->getOperand();
          const auto& fes = gf.getFiniteElementSpace();
          decltype(auto) fe  = fes.getFiniteElement(d, i);
          const auto& rc  = p.getReferenceCoordinates();

          // div(u) = trace( Jinv^T * sum_a u_a * grad_ref(phi_a) )
          // Here u_a is vector-valued (vdim = d in typical use), but we do not
          // assume vdim==d; we take dot with mapped gradients component-wise.
          const auto JinvT = p.getJacobianInverse().transpose();

          static thread_local SpatialVectorType s_grad_phys;
          s_grad_phys.resize(d);

          out = ScalarType(0);

          const size_t vdim = fes.getVectorDimension();

          for (size_t a = 0; a < fe.getCount(); ++a)
          {
            const auto& basis = fe.getBasis(a); // vector basis
            // basis.getDerivative<1>(comp, j)(rc) exists but expensive in H1 vector element;
            // use scalar reference gradients and component structure if available is hard.
            // So: compute physical gradients component-wise via basis.getDerivative<1>(i,j).
            // However your vector H1 basis is "scalar basis in one component":
            // DerivativeFunction<1>(i,j) returns dphi/dxj if i == local%vdim else 0.
            // We can exploit that to compute divergence directly with one scalar basis gradient.

            const size_t comp = a % vdim;
            const size_t alpha = a / vdim;

            // Reference gradient of scalar basis alpha:
            // NOTE: we have access to scalar FE through H1Element<K,ScalarType>.
            const auto& feS = H1Element<K, ScalarType>(polytope.getGeometry());
            static thread_local SpatialVectorType s_ref;
            s_ref.resize(d);
            for (size_t j = 0; j < d; ++j)
              s_ref(j) = feS.getBasis(alpha).template getDerivative<1>(j)(rc);

            s_grad_phys = JinvT * s_ref;

            // divergence contribution is u_comp * dphi/dx_comp if comp < d
            // (only makes sense when vdim == d; if vdim != d, we only sum over min(vdim,d))
            if (comp < d)
            {
              const auto& u_a = gf[fes.getGlobalIndex({d, i}, a)];
              // u_a is vector coefficient for basis "component", but in your layout
              // GridFunction DOF for vector space is scalar (component value).
              // If gf[...] returns ScalarType (typical), then:
              out += u_a * s_grad_phys(comp);
            }
          }
        }
      }

      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope& geom) const noexcept
      {
        const size_t k = H1Element<K, ScalarType>(geom.getGeometry()).getOrder();
        return (k == 0) ? 0 : (k - 1);
      }

      Div* copy() const noexcept override
      {
        return new Div(*this);
      }
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for Div of an H1 GridFunction
   */
  template <size_t K, class Scalar, class Data, class Mesh>
  Div(const GridFunction<H1<K, Math::Vector<Scalar>, Mesh>, Data>&)
    -> Div<GridFunction<H1<K, Math::Vector<Scalar>, Mesh>, Data>>;

  /**
   * @ingroup DivSpecializations
   * @brief Divergence of an H1 vector ShapeFunction.
   *
   * Represents @f$ \nabla \cdot \mathbf{v} @f$ for H1 test/trial functions.
   * Used in weak formulations like:
   * @f[
   *   \int_\Omega p (\nabla \cdot \mathbf{v}) \, dx
   * @f]
   * for pressure-velocity coupling in Stokes/Navier-Stokes equations.
   */
  template <size_t K, class NestedDerived, class Number, class Mesh, ShapeFunctionSpaceType Space>
  class Div<ShapeFunction<NestedDerived, H1<K, Math::Vector<Number>, Mesh>, Space>> final
    : public ShapeFunctionBase<Div<ShapeFunction<NestedDerived, H1<K, Math::Vector<Number>, Mesh>, Space>>>
  {
    public:
      using FESType = H1<K, Math::Vector<Number>, Mesh>;
      static constexpr ShapeFunctionSpaceType SpaceType = Space;

      using OperandType = ShapeFunction<NestedDerived, FESType, SpaceType>;
      using Parent = ShapeFunctionBase<Div<OperandType>, FESType, SpaceType>;
      using ScalarType = Number;

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

        std::vector<ScalarType> div_phys; // size = ndof (vector dofs)
        Key key;
      };

      Div(const OperandType& u)
        : Parent(u.getFiniteElementSpace()),
          m_u(u),
          m_ip(nullptr)
      {}

      Div(const Div& other)
        : Parent(other),
          m_u(other.m_u),
          m_ip(nullptr),
          m_cache(other.m_cache)
      {}

      Div(Div&& other)
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
      size_t getDOFs(const Geometry::Polytope& polytope) const
      {
        return getOperand().getDOFs(polytope);
      }

      const IntegrationPoint& getIntegrationPoint() const
      {
        assert(m_ip);
        return *m_ip;
      }

      Div& setIntegrationPoint(const IntegrationPoint& ip)
      {
        m_ip = &ip;

        const auto& p  = ip.getPoint();
        const auto& qf = ip.getQuadratureFormula();
        const size_t qp = ip.getIndex();

        const auto& poly = p.getPolytope();
        const size_t d   = poly.getDimension();
        const Index  cell = poly.getIndex();
        const auto geom = poly.getGeometry();

        typename Cache::Key key;
        key.geom  = geom;
        key.dim   = d;
        key.cell  = cell;
        key.qf    = &qf;
        key.qp    = qp;
        key.valid = true;

        if (m_cache.key == key)
          return *this;

        m_cache.key = key;

        const auto& fes = this->getFiniteElementSpace();
        const size_t vdim = fes.getVectorDimension();

        // Scalar FE for reference gradients
        decltype(auto) feS = H1Element<K, ScalarType>(geom);

        const size_t nscalar = feS.getCount();
        const size_t ndof = vdim * nscalar;

        if (m_cache.div_phys.size() != ndof)
          m_cache.div_phys.resize(ndof);

        const auto& tab = feS.getTabulation(qf);
        const auto JinvT = p.getJacobianInverse().transpose();

        static thread_local SpatialVectorType s_ref;
        static thread_local SpatialVectorType s_phys;
        s_ref.resize(d);
        s_phys.resize(d);

        // For each vector dof (alpha, comp):
        // div( phi_alpha e_comp ) = d/dx_comp phi_alpha  (comp < d) else 0
        for (size_t alpha = 0; alpha < nscalar; ++alpha)
        {
          const auto gref = tab.getGradient(qp, alpha);

          for (size_t j = 0; j < d; ++j)
            s_ref(j) = gref[j];

          s_phys = JinvT * s_ref;

          for (size_t comp = 0; comp < vdim; ++comp)
          {
            const size_t local = alpha * vdim + comp;
            m_cache.div_phys[local] = (comp < d) ? s_phys(comp) : ScalarType(0);
          }
        }

        return *this;
      }

      ScalarType getBasis(size_t local) const
      {
        assert(m_cache.key);
        assert(local < m_cache.div_phys.size());
        return m_cache.div_phys[local];
      }

      constexpr
      const auto& getFiniteElementSpace() const
      {
        return getOperand().getFiniteElementSpace();
      }

      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope& geom) const noexcept
      {
        const auto k = getOperand().getOrder(geom);
        if (!k.has_value())
          return std::nullopt;
        return (*k == 0) ? 0 : (*k - 1);
      }

      Div* copy() const noexcept override
      {
        return new Div(*this);
      }

    private:
      std::reference_wrapper<const OperandType> m_u;
      const IntegrationPoint* m_ip;

      Cache m_cache;
  };

  template <size_t K, class NestedDerived, class Number, class Mesh, ShapeFunctionSpaceType Space>
  Div(const ShapeFunction<NestedDerived, H1<K, Math::Vector<Number>, Mesh>, Space>&)
    -> Div<ShapeFunction<NestedDerived, H1<K, Math::Vector<Number>, Mesh>, Space>>;
}

#endif
