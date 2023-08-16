#ifndef RODIN_VARIATIONAL_BOUNDARYNORMAL_H
#define RODIN_VARIATIONAL_BOUNDARYNORMAL_H

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/PolytopeTransformation.h"

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
        : m_sdim(mesh.getSpaceDimension())
      {
        assert(m_sdim > 0);
      }

      BoundaryNormal(const BoundaryNormal& other)
        : Parent(other),
          m_sdim(other.m_sdim)
      {}

      BoundaryNormal(BoundaryNormal&& other)
        : Parent(std::move(other)),
          m_sdim(std::move(other.m_sdim))
      {}

      inline
      constexpr
      size_t getDimension() const
      {
        return m_sdim;
      }

      inline
      Math::SpatialVector getValue(const Geometry::Point& p) const
      {
        Math::SpatialVector res;
        getValue(res, p);
        return res;
      }

      void getValue(Math::SpatialVector& res, const Geometry::Point& p) const
      {
        const auto& polytope = p.getPolytope();
        const auto& d = polytope.getDimension();
        const auto& i = polytope.getIndex();
        const auto& mesh = polytope.getMesh();
        assert(d == mesh.getDimension() - 1);
        assert(mesh.isBoundary(i));
        const auto& jacobian = p.getJacobian();
        res.resize(m_sdim);
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
        }

        const auto& incidence = mesh.getConnectivity().getIncidence({ d, d + 1 }, i);
        assert(incidence.size() == 1);
        auto pit = mesh.getPolytope(d + 1, *incidence.begin());
        for (auto vit = pit->getVertex(); vit; ++vit)
        {
          const auto v = vit->getCoordinates() - polytope.getVertex()->getCoordinates();
          if (res.dot(v) > 0)
          {
            res *= -1;
            break;
          }
        }
        res.normalize();
      }

      inline BoundaryNormal* copy() const noexcept override
      {
        return new BoundaryNormal(*this);
      }

    private:
      const size_t m_sdim;
  };
}

#endif

