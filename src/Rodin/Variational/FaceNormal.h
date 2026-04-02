/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file FaceNormal.h
 * @brief Unit normal vector on mesh faces for DG and face-based formulations.
 *
 * This file defines the FaceNormal class, which represents a unit normal
 * vector on codimension-one mesh entities (faces in 3D, edges in 2D). Face
 * normals are fundamental in face-based finite element formulations, notably
 * in Discontinuous Galerkin (DG) methods where numerical fluxes, jumps, and
 * averages depend on a consistent face orientation.
 *
 * ## Mathematical Foundation
 *
 * Let @f$ F @f$ be a codimension-one polytope of the mesh. A face normal is a
 * vector field @f$ \mathbf{n}_F @f$ such that
 * @f[
 *   \|\mathbf{n}_F\| = 1.
 * @f]
 * Geometrically, the normal is orthogonal to the tangent space of the face.
 *
 * For an interior face shared by two adjacent cells @f$ K_1 @f$ and
 * @f$ K_2 @f$, the normal is defined only up to sign unless one chooses which
 * adjacent cell defines the outward direction. In trace-based formulations,
 * this side is determined by the chosen trace domain.
 *
 * For a boundary face, there is a unique adjacent cell @f$ K @f$, and the
 * normal is oriented automatically outward from that cell.
 *
 * ## Orientation Convention
 *
 * The orientation is chosen as follows:
 *
 * - **Boundary face (one adjacent cell)**:
 *   the normal is oriented automatically outward from the unique incident cell.
 *
 * - **Interior face (two adjacent cells)**:
 *   the normal is oriented outward from the unique adjacent cell whose
 *   attribute belongs to the trace domain.
 *
 * Consequently, on interior faces the normal is trace-dependent, while on
 * boundary faces it is determined geometrically from the unique incident cell.
 *
 * ## Trace Semantics
 *
 * The FaceNormal object supports both boundary and interior codimension-one
 * entities, but the meaning differs slightly:
 *
 * - On a **boundary face**, no trace domain is required, since there is only
 *   one adjacent cell and the outward direction is unambiguous.
 *
 * - On an **interior face**, a trace domain is required to select which side
 *   of the interface defines the outward direction. If no trace domain is
 *   provided, or if none of the adjacent cells matches it, the normal is
 *   considered undefined and an exception is raised.
 *
 * This behavior makes FaceNormal suitable for DG interface terms where one
 * needs a normal consistent with a chosen side of the interface.
 *
 * ## Applications in DG Methods
 *
 * Typical uses include:
 *
 * - **Numerical fluxes**:
 *   @f$ \hat{f}(u^+, u^-, \mathbf{n}) @f$
 * - **Jump terms**:
 *   @f$ [\![u]\!] = u^+ \mathbf{n}^+ + u^- \mathbf{n}^- @f$
 * - **Average-normal couplings**:
 *   @f$ \{\!\{\nabla u\}\!\} \cdot \mathbf{n} @f$
 * - **Penalty terms**:
 *   @f$ \frac{\sigma}{h} [\![u]\!] \cdot [\![v]\!] @f$
 *
 * ## Difference from BoundaryNormal
 *
 * - **FaceNormal**:
 *   defined on all codimension-one mesh entities; on interior faces its
 *   orientation may depend on the trace domain.
 *
 * - **BoundaryNormal**:
 *   defined only on the boundary of the domain and always oriented outward.
 *
 * ## Manifold Assumption
 *
 * This implementation assumes that a codimension-one mesh entity is incident
 * to either:
 *
 * - exactly one cell (boundary face), or
 * - exactly two cells (interior face).
 *
 * Non-manifold situations with more than two incident cells are not supported.
 *
 * ## Usage Example
 * ```cpp
 * FaceNormal n(mesh);
 *
 * // DG numerical flux on interfaces
 * auto flux = InterfaceIntegral(Average(u) * Dot(n.traceOf(omega), Jump(v)));
 *
 * // Normal component of a vector field on a boundary region
 * auto un = Dot(velocity, n);
 * ```
 *
 * @see BoundaryNormal, InterfaceIntegral, Jump, Average
 */
#ifndef RODIN_VARIATIONAL_FACENORMAL_H
#define RODIN_VARIATIONAL_FACENORMAL_H

#include <Eigen/Geometry>

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Variational/Exceptions/UndeterminedTraceDomainException.h"

#include "ForwardDecls.h"
#include "VectorFunction.h"

namespace Rodin::Variational
{
  /**
   * @ingroup RodinVariational
   * @brief Unit normal vector on codimension-one mesh entities.
   *
   * FaceNormal represents a unit normal vector on mesh faces and is intended
   * for face integral formulations such as DG methods. On boundary faces, the
   * normal is oriented automatically outward from the unique adjacent cell. On
   * interior faces, the orientation is determined by the selected trace domain.
   */
  class FaceNormal : public VectorFunctionBase<Real, FaceNormal>
  {
    public:
      using ScalarType = Real;
      using RangeType = Math::Vector<ScalarType>;
      using SpatialVectorType = Math::SpatialVector<ScalarType>;
      using Parent = VectorFunctionBase<ScalarType, FaceNormal>;

      using Parent::traceOf;

      /**
       * @brief Constructs a face normal field on the given mesh.
       *
       * @param mesh Underlying mesh supporting codimension-one entities.
       */
      FaceNormal(const Geometry::MeshBase& mesh)
        : m_sdim(mesh.getSpaceDimension())
      {
        assert(m_sdim > 0);
      }

      FaceNormal(const FaceNormal& other)
        : Parent(other),
          m_sdim(other.m_sdim)
      {}

      FaceNormal(FaceNormal&& other)
        : Parent(std::move(other)),
          m_sdim(std::move(other.m_sdim))
      {}

      constexpr
      size_t getDimension() const
      {
        return m_sdim;
      }

      decltype(auto) getValue(const Geometry::Point& p) const
      {
        static thread_local RangeType s_res;

        const auto& polytope = p.getPolytope();
        const auto& vs = polytope.getVertices();
        const auto  d = polytope.getDimension();
        const auto  i = polytope.getIndex();
        const auto& mesh = polytope.getMesh();
        const auto& jacobian = p.getJacobian();

        assert(d == mesh.getDimension() - 1);

        SpatialVectorType res;
        res.resize(m_sdim);

        if (jacobian.rows() == 2)
        {
          res[0] = jacobian(1, 0);
          res[1] = -jacobian(0, 0);
        }
        else if (jacobian.rows() == 3)
        {
          if (jacobian.cols() == 1)
          {
            const Index v1 = vs[0];
            const Index v2 = vs[1];

            const Math::SpatialVector<ScalarType> a =
              mesh.getVertexCoordinates(v1) - mesh.getVertexCoordinates(v2);

            Math::SpatialVector<ScalarType> n(3);
            n[0] = jacobian(1, 0);
            n[1] = -jacobian(0, 0);
            n[2] = jacobian(2, 0);
            n = n.cross(a);
            n.normalize();

            res = n.cross(a) + n * (n.dot(a));
          }
          else if (jacobian.cols() == 2)
          {
            res[0] =
              jacobian(1, 0) * jacobian(2, 1) - jacobian(2, 0) * jacobian(1, 1);
            res[1] =
              jacobian(2, 0) * jacobian(0, 1) - jacobian(0, 0) * jacobian(2, 1);
            res[2] =
              jacobian(0, 0) * jacobian(1, 1) - jacobian(1, 0) * jacobian(0, 1);
          }
          else
          {
            assert(false);
            res.setConstant(Math::nan<ScalarType>());
          }
        }
        else
        {
          assert(false);
          res.setConstant(Math::nan<ScalarType>());
        }

        const auto& incidence = mesh.getConnectivity().getIncidence({ d, d + 1 }, i);
        assert(incidence.size() == 1 || incidence.size() == 2);

        const auto orientFromCell =
          [&](const Geometry::Polytope& cellPoly)
          {
            const auto xf = p.getCoordinates();

            const auto rc =
              Geometry::Polytope::Traits(cellPoly.getGeometry()).getCentroid();
            Geometry::Point pc(cellPoly, rc);
            const auto xc = pc.getCoordinates();

            if (res.dot(xc - xf) > 0)
              res *= ScalarType(-1);

            res.normalize();
          };

        if (incidence.size() == 1)
        {
          const Index cell = incidence[0];
          const auto pit = mesh.getPolytope(d + 1, cell);
          assert(pit);
          orientFromCell(*pit);
        }
        else
        {
          const auto& traceDomain = getTraceDomain();
          if (traceDomain.size() == 0)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "No trace domain provided: "
              << Alert::Notation::Predicate(true, "getTraceDomain().size() == 0")
              << ". FaceNormal on an interior face with two adjacent cells is "
                 "undefined without a trace domain."
              << Alert::Raise;
          }

          bool matched = false;
          for (const Index cell : incidence)
          {
            const auto pit = mesh.getPolytope(d + 1, cell);
            assert(pit);

            const Optional<Geometry::Attribute> ca = pit->getAttribute();
            if (!ca || !traceDomain.contains(*ca))
              continue;

            orientFromCell(*pit);
            matched = true;
            break;
          }

          if (!matched)
          {
            UndeterminedTraceDomainException(
              *this, __func__, { d, i }, traceDomain.begin(), traceDomain.end()).raise();
          }
        }

        s_res = res.getData().head(static_cast<Eigen::Index>(m_sdim));
        return s_res;
      }

      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope&) const noexcept
      {
        return 0;
      }

      FaceNormal* copy() const noexcept override
      {
        return new FaceNormal(*this);
      }

    private:
      const size_t m_sdim;
  };
}

#endif
