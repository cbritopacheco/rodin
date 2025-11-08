#ifndef RODIN_VARIATIONAL_P1_DERIVATIVE_H
#define RODIN_VARIATIONAL_P1_DERIVATIVE_H

#include "Rodin/Math/Vector.h"
#include "Rodin/Variational/Derivative.h"

namespace Rodin::FormLanguage
{
  template <class Range, class Data, class Mesh>
  struct Traits<Variational::Derivative<Variational::GridFunction<Variational::P1<Range, Mesh>, Data>>>
  {
    using FESType = Variational::P1<Range, Mesh>;

    using OperandType = Variational::GridFunction<FESType, Data>;

    using RangeType = Range;
  };

  template <class NestedDerived, class Range, class Mesh, Variational::ShapeFunctionSpaceType Space>
  struct Traits<
    Variational::Derivative<
      Variational::ShapeFunction<NestedDerived, Variational::P1<Range, Mesh>, Space>>>
  {
    using FESType = Variational::P1<Range, Mesh>;
    static constexpr Variational::ShapeFunctionSpaceType SpaceType = Space;

    using OperandType = Variational::ShapeFunction<NestedDerived, FESType, SpaceType>;

    using RangeType = Range;
  };
}

namespace Rodin::Variational
{
  /**
   * @ingroup DerivativeSpecializations
   * @brief Derivative of a P1 GridFunction
   */
  template <class Range, class Data, class Mesh>
  class Derivative<GridFunction<P1<Range, Mesh>, Data>> final
    : public DerivativeBase<GridFunction<P1<Range, Mesh>, Data>, Derivative<GridFunction<P1<Range, Mesh>, Data>>>
  {
    public:
      using FESType = P1<Range, Mesh>;

      using RangeType = Range;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using SpatialVectorType = Math::SpatialVector<ScalarType>;

      using OperandType = GridFunction<FESType, Data>;

      using Parent = DerivativeBase<OperandType, Derivative<OperandType>>;

      /**
       * @brief Constructs the derivative of an @f$ \mathbb{P}^1 @f$ function
       * @f$ u @f$.
       * @param[in] u P1 GridFunction
       */
      Derivative(size_t i, const OperandType& u)
        : Parent(u),
          m_i(i)
      {}

      /**
       * @brief Copy constructor
       */
      Derivative(const Derivative& other)
        : Parent(other),
          m_i(other.m_i)
      {}

      /**
       * @brief Move constructor
       */
      Derivative(Derivative&& other)
        : Parent(std::move(other)),
          m_i(std::move(other.m_i))
      {}

      void interpolate(ScalarType& out, const Geometry::Point& p) const
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
                << ". Derivative at an interface with no trace domain is undefined."
                << Alert::Raise;
            }
            else
            {
              for (auto& idx : inc)
              {
                const auto& tracePolytope = mesh.getPolytope(meshDim, idx);
                if (traceDomain.count(tracePolytope->getAttribute()))
                {
                  Math::SpatialPoint rc;
                  tracePolytope->getTransformation().inverse(rc, pc);
                  const Geometry::Point np(*tracePolytope, std::cref(rc), pc);
                  this->interpolate(out, np);
                  return;
                }
              }
              UndeterminedTraceDomainException(
                  *this, __func__, {d, i}, traceDomain.begin(), traceDomain.end()) << Alert::Raise;
            }
            return;
          }
        }
        else // Evaluating on a cell
        {
          assert(d == mesh.getDimension());
          const auto& gf = this->getOperand();
          const auto& fes = gf.getFiniteElementSpace();
          const auto& fe = fes.getFiniteElement(d, i);
          const auto& rc = p.getReferenceCoordinates();
          SpatialVectorType grad(d);
          SpatialVectorType res(d);
          res.setZero();
          for (size_t local = 0; local < fe.getCount(); local++)
          {
            fe.getGradient(local)(grad, rc);
            res += gf.getValue(fes.getGlobalIndex({d, i}, local)) * grad;
          }
          out = (p.getJacobianInverse().transpose() * res).eval().coeff(m_i);
        }
      }

      Derivative* copy() const noexcept override
      {
        return new Derivative(*this);
      }

    private:
      size_t m_i;
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for Derivative of a P1 GridFunction
   */
  template <class Range, class Data, class Mesh>
  Derivative(size_t, const GridFunction<P1<Range, Mesh>, Data>&)
    -> Derivative<GridFunction<P1<Range, Mesh>, Data>>;
}

#endif
