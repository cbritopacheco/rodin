#include "Rodin/Alert/MemberFunctionException.h"

#include "GrundmannMoller.h"

namespace Rodin::QF
{
  /**
   * @internal
   */
  boost::multi_array<size_t, 2> GrundmannMoller::initSizes()
  {
    const auto extents = boost::extents[RODIN_MAXIMAL_SPACE_DIMENSION + 1][RODIN_QF_GRUNDMANNMOLLER_MAX_S + 1];
    boost::multi_array<size_t, 2> res(extents);
    for (size_t n = 0; n <= RODIN_MAXIMAL_SPACE_DIMENSION; n++)
    {
      for (size_t s = 0; s <= RODIN_QF_GRUNDMANNMOLLER_MAX_S; s++)
        res[n][s] = Math::binom(n + s + 1, n + 1);
    }
    return res;
  }

  /**
   * @internal
   */
  boost::multi_array<Math::Vector, 2> GrundmannMoller::initWeights()
  {
    const auto extents = boost::extents[RODIN_MAXIMAL_SPACE_DIMENSION + 1][RODIN_QF_GRUNDMANNMOLLER_MAX_S + 1];
    boost::multi_array<Math::Vector, 2> res(extents);
    for (size_t n = 0; n <= RODIN_MAXIMAL_SPACE_DIMENSION; n++)
    {
      // Handle case n == 0
      if (n == 0)
      {
        for (size_t s = 0; s <= RODIN_QF_GRUNDMANNMOLLER_MAX_S; s++)
          res[0][s] = Math::Vector{{1}};
        continue;
      }

      // Handle case n > 0
      Array<size_t> beta(n), sums(n);
      for (size_t s = 0; s <= RODIN_QF_GRUNDMANNMOLLER_MAX_S; s++)
      {
        const size_t d = 2 * s + 1;
        auto& weights = res[n][s];
        weights.resize(s_sizes[n][s]);
        size_t pt = 0;
        for (size_t i = 0; i <= s; i++)
        {
          Scalar w = Math::pow(2, -2 * Integer(s));
          w *= Math::pow(Integer(d + n - 2 * i), Integer(d));
          w /= Math::factorial(i) * Math::factorial(d + n - i);
          if (i % 2)
            w = -w;
          size_t k = s - i;
          beta.setZero();
          sums.setZero();
          bool doneBeta = false;
          while (true)
          {
            weights.coeffRef(pt++) = w;
            int j = 0;
            while (sums[j] == k)
            {
              beta[j++] = 0;
              assert(j >= 0);
              if (static_cast<size_t>(j) == n)
              {
                doneBeta = true;
                break;
              }
            }
            if (doneBeta)
              break;
            beta[j]++;
            sums[j]++;
            for (j--; j >= 0; j--)
              sums[j] = sums[j + 1];
          }
        }
      }
    }
    return res;
  }

  boost::multi_array<std::vector<Math::SpatialVector>, 2> GrundmannMoller::initPoints()
  {
    const auto extents = boost::extents[RODIN_MAXIMAL_SPACE_DIMENSION + 1][RODIN_QF_GRUNDMANNMOLLER_MAX_S + 1];
    boost::multi_array<std::vector<Math::SpatialVector>, 2> res(extents);
    for (size_t n = 0; n <= RODIN_MAXIMAL_SPACE_DIMENSION; n++)
    {
      // Handle case n == 0
      if (n == 0)
      {
        for (size_t s = 0; s <= RODIN_QF_GRUNDMANNMOLLER_MAX_S; s++)
          res[0][s] = { Math::SpatialVector{{0}} };
        continue;
      }

      Array<size_t> beta(n), sums(n);
      for (size_t s = 0; s <= RODIN_QF_GRUNDMANNMOLLER_MAX_S; s++)
      {
        const size_t d = 2 * s + 1;
        auto& ps = res[n][s];
        ps.resize(s_sizes[n][s]);
        size_t pt = 0;
        for (size_t i = 0; i <= s; i++)
        {
          size_t k = s - i;
          beta.setZero();
          sums.setZero();
          bool doneBeta = false;
          while (true)
          {
            ps[pt++] = (2 * beta + 1).cast<Scalar>() / (d + n - 2 * i);
            int j = 0;
            while (sums[j] == k)
            {
              beta[j++] = 0;
              assert(j >= 0);
              if (static_cast<size_t>(j) == n)
              {
                doneBeta = true;
                break;
              }
            }
            if (doneBeta)
              break;
            beta[j]++;
            sums[j]++;
            for (j--; j >= 0; j--)
              sums[j] = sums[j + 1];
          }
        }
      }
    }
    return res;
  }

  boost::multi_array<size_t, 2> GrundmannMoller::s_sizes = GrundmannMoller::initSizes();

  boost::multi_array<Math::Vector, 2> GrundmannMoller::s_weights = GrundmannMoller::initWeights();

  boost::multi_array<std::vector<Math::SpatialVector>, 2> GrundmannMoller::s_points = GrundmannMoller::initPoints();

  GrundmannMoller::GrundmannMoller(size_t s, Geometry::Polytope::Type geom)
    : Parent(geom),
      m_s(s),
      m_n(Geometry::Polytope::getGeometryDimension(geom)),
      m_order(2 * m_s + 1)
  {
    assert(m_order <= RODIN_QF_GRUNDMANNMOLLER_MAX_ORDER);
    if (!Geometry::Polytope::isSimplex(geom))
    {
      Alert::MemberFunctionException(*this, __func__)
        << "Polytope::Type " << Alert::Notation::Print(geom) << " is not a simplex."
        << Alert::Raise;
    }
    assert(Geometry::Polytope::isSimplex(geom));
  }

  size_t GrundmannMoller::getSize() const
  {
    return s_sizes[m_n][m_s];
  }

  Scalar GrundmannMoller::getWeight(size_t i) const
  {
    return s_weights[m_n][m_s].coeff(i);
  }

  const Math::SpatialVector& GrundmannMoller::getPoint(size_t i) const
  {
    return s_points[m_n][m_s][i];
  }
}

