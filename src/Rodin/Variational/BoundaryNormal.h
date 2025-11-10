/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file BoundaryNormal.h
 * @brief Outward unit normal vector on domain boundaries.
 *
 * This file defines the BoundaryNormal class, which represents the outward
 * unit normal vector field on the boundary of a domain. Boundary normals are
 * essential for imposing natural boundary conditions and computing boundary
 * integrals in finite element formulations.
 *
 * ## Mathematical Foundation
 * For a domain @f$ \Omega \subset \mathbb{R}^d @f$ with boundary @f$ \partial\Omega @f$,
 * the outward unit normal @f$ \mathbf{n} @f$ satisfies:
 * @f[
 *   \|\mathbf{n}\| = 1, \quad \mathbf{n} \cdot \mathbf{t} = 0
 * @f]
 * where @f$ \mathbf{t} @f$ is any tangent vector to @f$ \partial\Omega @f$.
 *
 * ## Applications
 * - **Neumann boundary conditions**: @f$ \nabla u \cdot \mathbf{n} = g @f$
 * - **Robin boundary conditions**: @f$ \alpha u + \beta \nabla u \cdot \mathbf{n} = g @f$
 * - **Flux calculations**: @f$ \int_{\partial\Omega} (\mathbf{v} \cdot \mathbf{n}) \, ds @f$
 * - **Integration by parts**: Boundary terms involving normal derivatives
 * - **DG methods**: Numerical flux terms with normal components
 *
 * ## Properties
 * - Orientation: Points outward from the domain
 * - Uniqueness: Well-defined except at corners/edges (where limits differ)
 * - Computation: Derived from surface geometry and element transformations
 *
 * ## Usage Example
 * ```cpp
 * Mesh mesh = /* ... */;
 * BoundaryNormal n(mesh);
 * 
 * // Neumann BC: ∫_Γ g(n·∇u) v ds
 * auto neumann = BoundaryIntegral(g * Dot(n, Grad(u)), v).on(boundary_attr);
 * 
 * // Normal derivative
 * auto normal_deriv = Dot(Grad(u), n);  // ∂u/∂n
 * ```
 */
#ifndef RODIN_VARIATIONAL_BOUNDARYNORMAL_H
#define RODIN_VARIATIONAL_BOUNDARYNORMAL_H

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/SubMesh.h"
#include "Rodin/Geometry/PolytopeTransformation.h"
#include "Rodin/Math/Common.h"
#include "Rodin/Variational/Exceptions/UndeterminedTraceDomainException.h"

#include "ForwardDecls.h"
#include "VectorFunction.h"

namespace Rodin::Variational
{
  /**
   * @ingroup RodinVariational
   * @brief Outward unit normal vector on domain boundaries.
   *
   * BoundaryNormal represents the outward-pointing unit normal vector field
   * on the boundary of a mesh. It is computed from the mesh geometry and
   * element transformations.
   */
  class BoundaryNormal final : public VectorFunctionBase<Real, BoundaryNormal>
  {
    public:
      using ScalarType = Real;

      using VectorType = Math::SpatialVector<ScalarType>;

      using Parent = VectorFunctionBase<ScalarType, BoundaryNormal>;

      using Parent::traceOf;

      using Parent::operator();

      /**
       * @brief Constructs the outward unit normal.
       */
      BoundaryNormal(const Geometry::MeshBase& mesh)
        : m_sdim(mesh.getSpaceDimension()),
          m_mesh(mesh)
      {
        assert(m_sdim > 0);
      }

      BoundaryNormal(const BoundaryNormal& other)
        : Parent(other),
          m_sdim(other.m_sdim),
          m_mesh(other.m_mesh)
      {}

      BoundaryNormal(BoundaryNormal&& other)
        : Parent(std::move(other)),
          m_sdim(std::move(other.m_sdim)),
          m_mesh(std::move(other.m_mesh))
      {}

      constexpr
      size_t getDimension() const
      {
        return m_sdim;
      }

      void interpolate(VectorType& res, const Geometry::Point& p) const
      {
        const auto& polytope = p.getPolytope();
        const auto& d = polytope.getDimension();
        const auto& i = polytope.getIndex();
        const auto& jacobian = p.getJacobian();
        const auto& mesh = polytope.getMesh();
        const size_t meshDim = mesh.getDimension();
        res.resize(m_sdim);
        if (d == meshDim - 2) // Evaluating on a codimension-2 element
        {
          const auto& conn = mesh.getConnectivity();
          const auto& inc = conn.getIncidence({ meshDim - 2, meshDim - 1 }, i);
          const auto& pc = p.getPhysicalCoordinates();
          assert(inc.size() == 1 || inc.size() == 2);
          if (inc.size() == 1)
          {
            const auto& tracePolytope = mesh.getPolytope(meshDim - 1, *inc.begin());
            VectorType rc;
            tracePolytope->getTransformation().inverse(rc, pc);
            const Geometry::Point np(*tracePolytope, std::cref(rc), pc);
            interpolate(res, np);
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
                << ". BoundaryNormal at a codimension-2 element with no trace domain is undefined."
                << Alert::Raise;
            }
            else
            {
              for (auto& idx : inc)
              {
                const auto& tracePolytope = mesh.getPolytope(meshDim - 1, idx);
                if (traceDomain.count(tracePolytope->getAttribute()))
                {
                  VectorType rc;
                  tracePolytope->getTransformation().inverse(rc, pc);
                  const Geometry::Point np(*tracePolytope, std::cref(rc), pc);
                  interpolate(res, np);
                  return;
                }
              }
              UndeterminedTraceDomainException(
                  *this, __func__, {d, i}, traceDomain.begin(), traceDomain.end()) << Alert::Raise;
            }
            return;
          }
        }
        else
        {
          assert(d == mesh.getDimension() - 1);
          assert(mesh.isBoundary(i));
          if (jacobian.rows() == 2)
          {
            res << jacobian(1, 0), -jacobian(0, 0);
          }
          else if (jacobian.rows() == 3)
          {
            res <<
              jacobian(1, 0) * jacobian(2, 1) - jacobian(2, 0) * jacobian(1, 1),
              jacobian(2, 0) * jacobian(0, 1) - jacobian(0, 0) * jacobian(2, 1),
              jacobian(0, 0) * jacobian(1, 1) - jacobian(1, 0) * jacobian(0, 1);
          }
          else
          {
            assert(false);
            res.setConstant(NAN);
            return;
          }

          const auto& incidence = mesh.getConnectivity().getIncidence({ d, d + 1 }, i);
          assert(incidence.size() == 1);
          auto pit = mesh.getPolytope(d + 1, *incidence.begin());
          Integer ori = -1;
          for (auto vit = pit->getVertex(); vit; ++vit)
          {
            const auto v = vit->getCoordinates() - polytope.getVertex()->getCoordinates();
            if (res.dot(v) < 0)
            {
              ori *= -1;
              break;
            }
          }
          res *= ori;
          res.normalize();
        }
      }

      decltype(auto) getValue(const Geometry::Point& p) const
      {
        static thread_local VectorType s_out;
        const auto& polytope = p.getPolytope();
        const auto& polytopeMesh = polytope.getMesh();
        if (polytopeMesh == m_mesh.get())
        {
          this->interpolate(s_out, p);
        }
        else if (const auto inclusion = m_mesh.get().inclusion(p))
        {
          this->interpolate(s_out, *inclusion);
        }
        else if (m_mesh.get().isSubMesh())
        {
          const auto& submesh = m_mesh.get().asSubMesh();
          const auto restriction = submesh.restriction(p);
          this->interpolate(s_out, *restriction);
        }
        else
        {
          s_out.setConstant(Math::nan<ScalarType>());
          assert(false);
        }
        return s_out;
      }

      BoundaryNormal* copy() const noexcept override
      {
        return new BoundaryNormal(*this);
      }

    private:
      const size_t m_sdim;
      std::reference_wrapper<const Geometry::MeshBase> m_mesh;
  };
}

#endif

