#ifndef RODIN_VARIATIONAL_NORMAL_H
#define RODIN_VARIATIONAL_NORMAL_H

#include "Rodin/Geometry/Mesh.h"

#include "ForwardDecls.h"
#include "VectorFunction.h"

namespace Rodin::Variational
{
  /**
   * @brief Outward unit normal.
   */
  class Normal : public VectorFunctionBase<Normal>
  {
    public:
      using Parent = VectorFunctionBase<Normal>;

      /**
       * @brief Constructs the outward unit normal.
       */
      constexpr
      Normal(size_t dimension)
        : m_dimension(dimension)
      {
        assert(dimension > 0);
      }

      constexpr
      Normal(const Normal& other)
        : Parent(other),
          m_dimension(other.m_dimension)
      {}

      constexpr
      Normal(Normal&& other)
        : Parent(std::move(other)),
          m_dimension(other.m_dimension)
      {}

      inline
      constexpr
      size_t getDimension() const
      {
        return m_dimension;
      }

      auto getValue(const Geometry::Point& p) const
      {
        const auto& simplex = p.getSimplex();
        const auto& mesh = simplex.getMesh();
        auto& trans = simplex.getTransformation();
        const auto& jac = trans.Jacobian();
        const auto* jacdata = jac.Data();
        Math::Vector value(m_dimension);

        if (jac.Height() == 2)
        {
          value(0) =  jacdata[1];
          value(1) = -jacdata[0];
        }
        else
        {
          value(0) = jacdata[1] * jacdata[5] - jacdata[2] * jacdata[4];
          value(1) = jacdata[2] * jacdata[3] - jacdata[0] * jacdata[5];
          value(2) = jacdata[0] * jacdata[4] - jacdata[1] * jacdata[3];
        }

        const Scalar norm = value.norm();
        value /= norm;
        assert(norm >= 0.0);
        assert(std::isfinite(norm));

        if ((mesh.getDimension() + 1 == mesh.getSpaceDimension()) &&
            simplex.getDimension() == mesh.getDimension())
        {
          // We are on an element of a d-mesh in (d + 1)-space.
          return value;
        }
        else if ((mesh.getDimension() == mesh.getSpaceDimension()) &&
            simplex.getDimension() == mesh.getDimension() - 1)
        {
          // Or we are on a face of a d-mesh in d-space
          if (mesh.isBoundary(simplex.getIndex()))
          {
            return value;
          }
          else
          {
            int el1 = -1;
            int el2 = -1;
            const auto& meshHandle = simplex.getMesh().getHandle();
            meshHandle.GetFaceElements(simplex.getIndex(), &el1, &el2);
            if (el1 >= 0 && getTraceDomain() == meshHandle.GetAttribute(el1))
            {
              return value;
            }
            else if (el2 >= 0 && getTraceDomain() == meshHandle.GetAttribute(el2))
            {
              value = -1.0 * value;
              return value;
            }
            else
            {
              assert(false);
              value = Math::Vector::Zero(m_dimension);
              return value;
            }
          }
        }
        else
        {
          assert(false);
          value = Math::Vector::Zero(m_dimension);
          return value;
        }
        return value;
      }

    private:
      const size_t m_dimension;
  };
}

#endif
