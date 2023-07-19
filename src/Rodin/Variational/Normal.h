#ifndef RODIN_VARIATIONAL_NORMAL_H
#define RODIN_VARIATIONAL_NORMAL_H

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/SimplexTransformation.h"

#include "ForwardDecls.h"
#include "VectorFunction.h"

namespace Rodin::Variational
{
  /**
   * @brief Outward unit normal.
   */
  class Normal final : public VectorFunctionBase<Normal>
  {
    public:
      using Parent = VectorFunctionBase<Normal>;

      /**
       * @brief Constructs the outward unit normal.
       */
      Normal(const Geometry::MeshBase& surface)
        : m_dimension(surface.getSpaceDimension())
      {
        assert(m_dimension > 0);
      }

      Normal(const Normal& other)
        : Parent(other),
          m_dimension(other.m_dimension)
      {}

      Normal(Normal&& other)
        : Parent(std::move(other)),
          m_dimension(std::move(other.m_dimension))
      {}

      inline
      constexpr
      size_t getDimension() const
      {
        return m_dimension;
      }

      Math::Vector getValue(const Geometry::Point& p) const
      {
        assert(p.getPolytope().getMesh().isSurface());
        const auto& jacobian = p.getJacobian();
        Math::Vector value(m_dimension);
        if (jacobian.rows() == 2)
        {
          value(0) =  jacobian(1, 0);
          value(1) = -jacobian(0, 0);
        }
        else if (jacobian.rows() == 3)
        {
          value(0) = jacobian(1, 0) * jacobian(2, 1) - jacobian(2, 0) * jacobian(1, 1);
          value(1) = jacobian(2, 0) * jacobian(0, 1) - jacobian(0, 0) * jacobian(2, 1);
          value(2) = jacobian(0, 0) * jacobian(1, 1) - jacobian(1, 0) * jacobian(0, 1);
        }
        else
        {
          assert(false);
          value.setConstant(NAN);
        }
        return value.normalized();
      }

      inline Normal* copy() const noexcept override
      {
        return new Normal(*this);
      }

    private:
      const size_t m_dimension;
  };
}

#endif
