#ifndef RODIN_VARIATIONAL_FACENORMAL_H
#define RODIN_VARIATIONAL_FACENORMAL_H

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/PolytopeTransformation.h"
#include "Rodin/Variational/Exceptions/UndeterminedTraceDomainException.h"

#include "ForwardDecls.h"
#include "VectorFunction.h"

namespace Rodin::Variational
{
  /**
   * @brief Outward unit normal on a face.
   */
  class FaceNormal : public VectorFunctionBase<FaceNormal>
  {
    public:
      using Parent = VectorFunctionBase<FaceNormal>;

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
          for (const Index i : incidence)
          {
            auto pit = mesh.getPolytope(d + 1, i);
            if (traceDomain.contains(pit->getAttribute()))
            {
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
              return;
            }
          }

          UndeterminedTraceDomainException(
              *this, __func__, {d, i}, traceDomain.begin(), traceDomain.end()).raise();
        }
      }

      inline
      constexpr
      FaceNormal& traceOf(Geometry::Attribute attr)
      {
        Parent::traceOf(attr);
        return *this;
      }

      inline
      constexpr
      FaceNormal& traceOf(const FlatSet<Geometry::Attribute>& attrs)
      {
        Parent::traceOf(attrs);
        return *this;
      }

      inline FaceNormal* copy() const noexcept override
      {
        return new FaceNormal(*this);
      }

    private:
      const size_t m_sdim;
  };
}

#endif

