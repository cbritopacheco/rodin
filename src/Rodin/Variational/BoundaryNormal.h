#ifndef RODIN_VARIATIONAL_BOUNDARYNORMAL_H
#define RODIN_VARIATIONAL_BOUNDARYNORMAL_H

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/SubMesh.h"
#include "Rodin/Geometry/PolytopeTransformation.h"
#include "Rodin/Variational/Exceptions/UndeterminedTraceDomainException.h"

#include "ForwardDecls.h"
#include "VectorFunction.h"

namespace Rodin::Variational
{
  /**
   * @brief Outward unit normal.
   */
  class BoundaryNormal final : public VectorFunctionBase<BoundaryNormal>
  {
    public:
      using Parent = VectorFunctionBase<BoundaryNormal>;

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

      inline
      constexpr
      size_t getDimension() const
      {
        return m_sdim;
      }

      inline
      constexpr
      BoundaryNormal& traceOf(Geometry::Attribute attr)
      {
        Parent::traceOf(attr);
        return *this;
      }

      inline
      constexpr
      BoundaryNormal& traceOf(const FlatSet<Geometry::Attribute>& attrs)
      {
        Parent::traceOf(attrs);
        return *this;
      }

      inline
      Math::SpatialVector getValue(const Geometry::Point& p) const
      {
        Math::SpatialVector res;
        getValue(res, p);
        return res;
      }

      void interpolate(Math::SpatialVector& res, const Geometry::Point& p) const
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
            const Math::SpatialVector rc = tracePolytope->getTransformation().inverse(pc);
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
                  const Math::SpatialVector rc = tracePolytope->getTransformation().inverse(pc);
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

      void getValue(Math::SpatialVector& out, const Geometry::Point& p) const
      {
        const auto& polytope = p.getPolytope();
        const auto& polytopeMesh = polytope.getMesh();
        if (polytopeMesh == m_mesh.get())
        {
          interpolate(out, p);
        }
        else
        {
          if (polytopeMesh.isSubMesh())
          {
            const auto& submesh = polytopeMesh.asSubMesh();
            assert(submesh.getParent() == fes.getMesh());
            interpolate(out, submesh.inclusion(p));
          }
          else if (m_mesh.get().isSubMesh())
          {
            const auto& submesh = m_mesh.get().asSubMesh();
            assert(submesh.getParent() == polytopeMesh);
            interpolate(out, submesh.restriction(p));
          }
          else
          {
            assert(false);
            out.setConstant(NAN);
          }
        }
      }

      inline BoundaryNormal* copy() const noexcept override
      {
        return new BoundaryNormal(*this);
      }

    private:
      const size_t m_sdim;
      std::reference_wrapper<const Geometry::MeshBase> m_mesh;
  };
}

#endif

