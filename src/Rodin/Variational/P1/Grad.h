/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P1_GRAD_H
#define RODIN_VARIATIONAL_P1_GRAD_H

/**
 * @file
 * @brief Gradient operator specialization for P1 (piecewise linear) functions.
 *
 * For P1 functions, the gradient is constant on each element:
 * @f[
 *   \nabla u|_K = \sum_{i=1}^{n+1} u_i \nabla \phi_i
 * @f]
 * where @f$ \phi_i @f$ are the P1 basis functions and @f$ n @f$ is the spatial dimension.
 */

#include "Rodin/Geometry/Mesh.h"
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
   * @brief Base class for Grad classes.
   */
  template <class Operand, class Derived>
  class GradBase;

  template <class Scalar, class Mesh, class Data>
  class Grad<GridFunction<P1<Scalar, Mesh>, Data>> final
    : public GradBase<GridFunction<P1<Scalar, Mesh>, Data>, Grad<GridFunction<P1<Scalar, Mesh>, Data>>>
  {
    public:
      using FESType = P1<Scalar, Mesh>;

      using RangeType = Math::SpatialVector<Scalar>;

      using ScalarType = typename FormLanguage::Traits<RangeType>::ScalarType;

      using SpatialVectorType = Math::SpatialVector<ScalarType>;

      using OperandType = GridFunction<FESType, Data>;

      using Parent = GradBase<OperandType, Grad<OperandType>>;

      /**
       * @brief Constructs the gradient of a P1 function @f$ u @f$.
       * @param[in] u P1 GridFunction (piecewise linear)
       *
       * @note The gradient is constant on each element for P1 functions.
       */
      Grad(const OperandType& u)
        : Parent(u)
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Grad object to copy
       */
      Grad(const Grad& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Grad object to move from
       */
      Grad(Grad&& other)
        : Parent(std::move(other))
      {}

      void interpolate(SpatialVectorType& out, const IntegrationPoint& ip) const
      {
        interpolate(out, ip.getPoint());
      }

      /**
       * @brief Interpolates the gradient at a given point.
       * @param[out] out Output spatial vector for gradient
       * @param[in] p Point at which to evaluate gradient
       *
       * Computes @f$ \nabla u(p) @f$ using P1 basis function gradients.
       * Handles evaluation on faces by projecting to adjacent cells.
       */
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
            Math::SpatialPoint rc;
            tracePolytope->getTransformation().inverse(rc, pc);
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
                const Optional<Geometry::Attribute> a = tracePolytope->getAttribute();
                if (!a || !traceDomain.contains(*a)) // or traceDomain.count(*a)
                  continue;
                Math::SpatialPoint rc;
                tracePolytope->getTransformation().inverse(rc, pc);
                const Geometry::Point np(*tracePolytope, std::cref(rc), pc);
                this->interpolate(out, np);
                return;
              }
              UndeterminedTraceDomainException(
                  *this, __func__, {d, i}, traceDomain.begin(), traceDomain.end()) << Alert::Raise;
            }
            return;
          }
        }
        else // Evaluating on a cell
        {
          SpatialVectorType res(d);
          res.setZero();

          assert(d == mesh.getDimension());
          const auto& gf = this->getOperand();
          const auto& fes = gf.getFiniteElementSpace();
          const auto& fe = fes.getFiniteElement(d, i);
          const auto& rc = p.getReferenceCoordinates();
          for (size_t local = 0; local < fe.getCount(); local++)
          {
            const auto& basis = fe.getBasis(local);
            basis.getGradient()(rc);
            res += gf[fes.getGlobalIndex({d, i}, local)] * basis.getGradient()(rc);
          }
          out = p.getJacobianInverse().transpose() * res;
        }
      }

      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope& geom) const noexcept
      {
        const size_t k = P1Element<ScalarType>(geom.getGeometry()).getOrder();
        return (k == 0) ? 0 : (k - 1);
      }

      Grad* copy() const noexcept override
      {
        return new Grad(*this);
      }
  };

  /**
   * @ingroup GradSpecializations
   * @brief Gradient of a ShapeFunction
   */
  template <class NestedDerived, class Scalar, class Mesh, ShapeFunctionSpaceType SpaceType>
  class Grad<ShapeFunction<NestedDerived, P1<Scalar, Mesh>, SpaceType>> final
    : public ShapeFunctionBase<Grad<ShapeFunction<NestedDerived, P1<Scalar, Mesh>, SpaceType>>>
  {
    public:
      using FESType = P1<Scalar, Mesh>;
      static constexpr ShapeFunctionSpaceType Space = SpaceType;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using SpatialVectorType = Math::SpatialVector<ScalarType>;

      using RangeType = Math::SpatialVector<ScalarType>;

      using OperandType = ShapeFunction<NestedDerived, FESType, Space>;

      using Parent = ShapeFunctionBase<Grad<OperandType>, FESType, Space>;

      struct Cache
      {
        struct CellKey
        {
          const void* mesh = nullptr;
          size_t d = 0;
          Index i = 0;
          Geometry::Polytope::Type geom = Geometry::Polytope::Type::Point;
          int transOrder = 1;
          bool valid = false;

          explicit operator bool() const noexcept { return valid; }

          bool operator==(const CellKey& o) const noexcept
          {
            if (!valid || !o.valid)
              return false;
            return mesh == o.mesh && d == o.d && i == o.i
                && geom == o.geom && transOrder == o.transOrder;
          }

          void operator=(std::initializer_list<int>) noexcept
          {
            valid = false;
            mesh = nullptr;
            d = 0;
            i = 0;
            geom = Geometry::Polytope::Type::Point;
            transOrder = 1;
          }
        };

        struct QpKey
        {
          const QF::QuadratureFormulaBase* qf = nullptr;
          size_t qp = 0;
          bool valid = false;

          explicit operator bool() const noexcept { return valid; }

          bool operator==(const QpKey& o) const noexcept
          {
            if (!valid || !o.valid)
              return false;
            return qf == o.qf && qp == o.qp;
          }

          void operator=(std::initializer_list<int>) noexcept
          {
            valid = false;
            qf = nullptr;
            qp = 0;
          }
        };

        // Cached physical gradients ∇_x φ_a (one per scalar basis function)
        std::vector<SpatialVectorType> grad;

        CellKey cellKey;
        QpKey qpKey;
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
        // Gradient has same number of scalar DOFs as the operand basis count.
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

        const auto& pt   = ip.getPoint();
        const auto& poly = pt.getPolytope();
        const auto& mesh = poly.getMesh();

        const size_t d    = poly.getDimension();
        const Index  i    = poly.getIndex();
        const auto   geom = poly.getGeometry();

        const size_t meshDim = mesh.getDimension();

        // If the integration point is on a face (d < meshDim), redirect the
        // gradient computation to the adjacent cell.  This ensures that
        // getBasis() returns the full spatial gradient ∇φ (from the cell
        // basis), not the lower-dimensional tangential derivative of the face
        // parametrization — which is what is needed for boundary integrals
        // such as the Nitsche / SIPG consistency and symmetry terms.
        if (d < meshDim)
        {
          const auto& conn = mesh.getConnectivity();
          const auto& inc  = conn.getIncidence({d, d + 1}, i);
          assert(!inc.empty());

          // Select adjacent cell.  For boundary faces there is only one
          // neighbour; for interior faces use the trace domain (if set).
          // Note: getTraceDomain() is only available on FunctionBase-derived
          // classes (e.g. Grad<GridFunction<...>>).  For ShapeFunction-based
          // Grad objects (which inherit from ShapeFunctionBase), no trace
          // domain is available, so we fall back to the first adjacent cell.
          Index cellIdx = inc[0];
          if (inc.size() > 1)
          {
            if constexpr (requires { this->getTraceDomain(); })
            {
              const auto& traceDomain = this->getTraceDomain();
              if (!traceDomain.empty())
              {
                for (size_t k = 0; k < inc.size(); ++k)
                {
                  const auto cellPoly = mesh.getPolytope(d + 1, inc[k]);
                  const Optional<Geometry::Attribute> attr = cellPoly->getAttribute();
                  if (attr && traceDomain.contains(*attr))
                  {
                    cellIdx = inc[k];
                    break;
                  }
                }
              }
            }
          }

          const auto cellPoly     = mesh.getPolytope(d + 1, cellIdx);
          const auto cellGeom     = cellPoly->getGeometry();
          const int  cellTransOrd = cellPoly->getTransformation().getOrder();

          // ---- cell key (keyed on the adjacent cell, not the face)
          typename Cache::CellKey ckey;
          ckey.mesh       = static_cast<const void*>(&mesh);
          ckey.d          = d + 1;       // cell dimension (= meshDim)
          ckey.i          = cellIdx;
          ckey.geom       = cellGeom;
          ckey.transOrder = cellTransOrd;
          ckey.valid      = true;

          const bool cell_changed = !(m_cache.cellKey == ckey);

          const size_t sdim   = d + 1;
          const size_t faceNv = Geometry::Polytope::Traits(geom).getVertexCount();
          const P1Element<ScalarType> cellFe(cellGeom);
          const size_t cellNv = cellFe.getCount();

          if (cell_changed)
          {
            m_cache.cellKey = ckey;
            m_cache.qpKey   = {};
            m_cache.grad.resize(faceNv);
            for (auto& gvec : m_cache.grad)
            {
              gvec.resize(sdim);
              gvec.setZero();
            }
          }

          // For P1 with affine maps the gradient is constant within the cell,
          // so a single computation per cell suffices.  For higher-order maps
          // (transOrder > 1) or tensor-product geometries recompute per QP.
          const bool tensor_ref =
            (cellGeom == Geometry::Polytope::Type::Quadrilateral) ||
            (cellGeom == Geometry::Polytope::Type::Wedge)         ||
            (cellGeom == Geometry::Polytope::Type::Hexahedron);
          const bool needs_qp = (cellTransOrd > 1) || tensor_ref;

          typename Cache::QpKey qkey;
          if (needs_qp)
          {
            qkey.qf    = &ip.getQuadratureFormula();
            qkey.qp    = ip.getIndex();
            qkey.valid = true;
          }
          else
          {
            qkey.qf    = nullptr;
            qkey.qp    = 0;
            qkey.valid = true;
          }

          const bool qp_changed = !(m_cache.qpKey == qkey);
          if (cell_changed || qp_changed)
          {
            m_cache.qpKey = qkey;

            // Map the face physical point to cell reference coordinates.
            const auto& pc = pt.getPhysicalCoordinates();
            Math::SpatialPoint cellRc;
            cellPoly->getTransformation().inverse(cellRc, pc);

            // Build a temporary cell point to obtain the cell Jacobian.
            const Geometry::Point cellPt(*cellPoly, std::cref(cellRc), pc);
            const auto JcellInvT = cellPt.getJacobianInverse().transpose();

            // Compute physical gradients for every cell DOF.
            std::vector<SpatialVectorType> cellGrad(cellNv, SpatialVectorType(sdim));
            for (size_t a = 0; a < cellNv; ++a)
            {
              Math::SpatialVector<ScalarType> ghat(sdim);
              for (size_t k = 0; k < sdim; ++k)
                ghat(k) = cellFe.getBasis(a).template getDerivative<1>(k)(cellRc);
              cellGrad[a] = JcellInvT * ghat;
            }

            // Build face-local → cell-local DOF correspondence by matching
            // vertex indices and fill the cached face-ordered gradients.
            const auto& faceVerts = conn.getPolytope(d, i);
            const auto& cellVerts = conn.getPolytope(d + 1, cellIdx);
            for (size_t k = 0; k < faceNv; ++k)
            {
              const Index v = faceVerts(k);
              bool found = false;
              for (size_t j = 0; j < cellNv; ++j)
              {
                if (cellVerts(j) == v)
                {
                  m_cache.grad[k] = cellGrad[j];
                  found = true;
                  break;
                }
              }
              assert(found); // every face vertex must belong to its adjacent cell
              (void)found;
            }
          }

          return *this;
        }

        // ---- Cell evaluation (d == meshDim): original logic -----------------
        const int transOrder = poly.getTransformation().getOrder();

        // ---- cell key: allocate/size once per cell
        typename Cache::CellKey ckey;
        ckey.mesh = static_cast<const void*>(&mesh);
        ckey.d = d;
        ckey.i = i;
        ckey.geom = geom;
        ckey.transOrder = transOrder;
        ckey.valid = true;

        const bool cell_changed = !(m_cache.cellKey == ckey);
        if (cell_changed)
        {
          m_cache.cellKey = ckey;
          m_cache.qpKey = {}; // invalidate qp cache

          const size_t nv = Geometry::Polytope::Traits(geom).getVertexCount();
          m_cache.grad.resize(nv);
          for (auto& gvec : m_cache.grad)
          {
            gvec.resize(d);
            gvec.setZero();
          }
        }

        // ---- decide if gradients depend on quadrature point
        const bool tensor_ref =
          (geom == Geometry::Polytope::Type::Quadrilateral) ||
          (geom == Geometry::Polytope::Type::Wedge) ||
          (geom == Geometry::Polytope::Type::Hexahedron);

        const bool needs_qp = (transOrder > 1) || tensor_ref;

        typename Cache::QpKey qkey;
        if (needs_qp)
        {
          qkey.qf = &ip.getQuadratureFormula();
          qkey.qp = ip.getIndex();
          qkey.valid = true;
        }
        else
        {
          // one state per cell
          qkey.qf = nullptr;
          qkey.qp = 0;
          qkey.valid = true;
        }

        const bool qp_changed = !(m_cache.qpKey == qkey);
        if (cell_changed || qp_changed)
        {
          m_cache.qpKey = qkey;

          const P1Element<ScalarType> fe(geom);
          const size_t nv = fe.getCount();

          const auto& qf = ip.getQuadratureFormula();
          const size_t qp = ip.getIndex();
          const auto& rc = qf.getPoint(qp);

          // J^{-T} at this integration point (constant for affine maps)
          const auto JinvT = pt.getJacobianInverse().transpose();

          // Compute physical gradients: ∇_x φ_a = J^{-T} ∇_hat φ_a
          for (size_t a = 0; a < nv; ++a)
          {
            // Reference gradient (size d). Build without using GradientFunction()
            // to avoid constructing thread_local vectors repeatedly.
            Math::SpatialVector<ScalarType> ghat(d);
            for (size_t k = 0; k < d; ++k)
              ghat(k) = fe.getBasis(a).template getDerivative<1>(k)(rc);

            m_cache.grad[a] = JinvT * ghat;
          }
        }

        return *this;
      }

      constexpr
      auto getBasis(size_t local) const
      {
        assert(m_cache.cellKey);
        assert(local < m_cache.grad.size());
        return m_cache.grad[local].getData().head(m_cache.cellKey.d);
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
