/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file FaceNormal.h
 * @brief Normal vector on mesh faces for DG methods.
 *
 * This file defines the FaceNormal class, which represents the unit normal
 * vector on mesh faces (both interior and boundary). Face normals are crucial
 * for Discontinuous Galerkin (DG) methods where numerical fluxes depend on
 * the normal direction at element interfaces.
 *
 * ## Mathematical Foundation
 * For a face @f$ F @f$ shared by elements @f$ K^+ @f$ and @f$ K^- @f$, the
 * face normal @f$ \mathbf{n} @f$ satisfies:
 * @f[
 *   \|\mathbf{n}\| = 1
 * @f]
 *
 * ## Orientation Convention
 * The normal direction is chosen according to a consistent orientation:
 * - For interior faces: points from @f$ K^+ @f$ to @f$ K^- @f$
 * - For boundary faces: points outward from the domain
 *
 * ## Applications in DG Methods
 * - **Numerical fluxes**: @f$ \hat{f}(\mathbf{u}^+, \mathbf{u}^-, \mathbf{n}) @f$
 * - **Jump terms**: @f$ [\![u]\!] = u^+ \mathbf{n}^+ + u^- \mathbf{n}^- @f$
 * - **Average terms**: @f$ \{\!\{\nabla u\}\!\} \cdot \mathbf{n} @f$
 * - **Penalty terms**: @f$ \frac{\sigma}{h} [\![u]\!] \cdot [\![v]\!] @f$
 *
 * ## Difference from BoundaryNormal
 * - **FaceNormal**: Defined on all faces (interior + boundary)
 * - **BoundaryNormal**: Only defined on domain boundary
 *
 * ## Usage Example
 * ```cpp
 * FaceNormal n(mesh);
 * 
 * // DG numerical flux on interior faces
 * auto flux = InterfaceIntegral(Average(u) * Dot(n, Jump(v)));
 * 
 * // Normal component of vector field
 * auto normal_component = Dot(velocity, n);
 * ```
 *
 * @see BoundaryNormal, InterfaceIntegral, Jump, Average
 */
#ifndef RODIN_VARIATIONAL_FACENORMAL_H
#define RODIN_VARIATIONAL_FACENORMAL_H

#include <Eigen/Geometry>

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Variational/Exceptions/UndeterminedTraceDomainException.h"

#include "ForwardDecls.h"
#include "VectorFunction.h"

namespace Rodin::Variational
{
  /**
   * @ingroup RodinVariational
   * @brief Unit normal vector on mesh faces.
   *
   * FaceNormal represents the unit normal vector on mesh faces, essential for
   * DG methods and face integral computations. The normal orientation follows
   * a consistent convention across the mesh.
   */
  class FaceNormal : public VectorFunctionBase<Real, FaceNormal>
  {
    public:
      using ScalarType = Real;

      using SpatialVectorType = Math::SpatialVector<ScalarType>;

      using Parent = VectorFunctionBase<ScalarType, FaceNormal>;

      using Parent::traceOf;

      /**
       * @brief Constructs the outward unit on a face.
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
        static thread_local SpatialVectorType s_res;
        const auto& polytope = p.getPolytope();
        const auto& vs = polytope.getVertices();
        const auto  d  = polytope.getDimension();
        const auto  i  = polytope.getIndex();
        const auto& mesh = polytope.getMesh();
        assert(d == mesh.getDimension() - 1);
        const auto& jacobian = p.getJacobian();

        s_res.resize(m_sdim);
        if (jacobian.rows() == 2)
        {
          s_res << jacobian(1,0), -jacobian(0,0);
        }
        else if (jacobian.rows() == 3)
        {
          if (jacobian.cols() == 1)
          {
            const Index v1 = vs[0];
            const Index v2 = vs[1];
            Eigen::Vector3<ScalarType> a =
                mesh.getVertexCoordinates(v1) - mesh.getVertexCoordinates(v2);
            Eigen::Vector3<ScalarType> n;
            n << jacobian(1, 0), -jacobian(0, 0), jacobian(2, 0);
            n = n.cross(a);
            n.stableNormalize();
            s_res = n.cross(a) + n * (n.dot(a));
          }
          else if (jacobian.cols() == 2)
          {
            s_res <<
              jacobian(1, 0) * jacobian(2, 1) - jacobian(2, 0) * jacobian(1, 1),
              jacobian(2, 0) * jacobian(0, 1) - jacobian(0, 0) * jacobian(2, 1),
              jacobian(0, 0) * jacobian(1, 1) - jacobian(1, 0) * jacobian(0, 1);
          }
        }
        else
        {
          assert(false);
          s_res.setConstant(Math::nan<ScalarType>());
        }

        const auto& incidence = mesh.getConnectivity().getIncidence({d, d+1}, i);
        assert(incidence.size() == 1 || incidence.size() == 2);

        const auto& traceDomain = getTraceDomain();
        if (traceDomain.size() == 0)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "No trace domain provided: "
            << Alert::Notation::Predicate(true, "getTraceDomain().size() == 0")
            << ". FaceNormal at an interface with no trace domain is undefined."
            << Alert::Raise;
        }
        else
        {
          bool matched = false;
          for (const Index cell : incidence)
          {
            auto pit = mesh.getPolytope(d + 1, cell);
            if (traceDomain.contains(pit->getAttribute()))
            {
              Integer ori = -1;
              for (auto vit = pit->getVertex(); vit; ++vit)
              {
                const auto v = vit->getCoordinates() - polytope.getVertex()->getCoordinates();
                if (s_res.dot(v) < 0)
                {
                  ori *= -1;
                  break;
                }
              }
              s_res *= ori;
              s_res.normalize();
              matched = true;
              break;
            }
          }

          if (!matched)
          {
            UndeterminedTraceDomainException(
              *this, __func__, {d, i}, traceDomain.begin(), traceDomain.end()).raise();
          }
        }
        return s_res;
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

