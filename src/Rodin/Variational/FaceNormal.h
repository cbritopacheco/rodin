/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FACENORMAL_H
#define RODIN_VARIATIONAL_FACENORMAL_H

/**
 * @file
 * @brief Outward unit normal vector on mesh faces.
 *
 * This file defines the outward unit normal vector @f$ \mathbf{n} @f$ on interior
 * and boundary faces of a mesh. For interior faces, the orientation depends on
 * the choice of "positive" and "negative" sides of the face.
 *
 * ## Mathematical Definition
 * The face normal @f$ \mathbf{n} @f$ is perpendicular to the face and satisfies:
 * @f[
 *   \|\mathbf{n}\| = 1
 * @f]
 *
 * ## Usage
 * Face normals are essential in:
 * - Discontinuous Galerkin (DG) methods
 * - Interior penalty methods
 * - Flux computations across element interfaces
 */

#include <Eigen/Geometry>

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Variational/Exceptions/UndeterminedTraceDomainException.h"

#include "ForwardDecls.h"
#include "VectorFunction.h"

namespace Rodin::Variational
{
  /**
   * @brief Outward unit normal on a face.
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

