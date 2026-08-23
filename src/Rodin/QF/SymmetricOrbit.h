/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_QF_SYMMETRICORBIT_H
#define RODIN_QF_SYMMETRICORBIT_H

/**
 * @file
 * @brief Defines SymmetricOrbit, the barycentric symmetry orbit of a
 * fully symmetric simplex quadrature rule.
 */

#include <array>
#include <vector>
#include <algorithm>

#include "Rodin/Math.h"
#include "Rodin/Geometry/Polytope.h"

namespace Rodin::QF
{
  /**
   * @ingroup RodinQuadrature
   * @brief One symmetry orbit of a fully symmetric simplex quadrature rule.
   *
   * A rule on a simplex is *fully symmetric* when its point set is invariant
   * under every permutation of the barycentric coordinates. Such a rule is
   * therefore a union of orbits, each orbit being the set of distinct
   * permutations of a single barycentric tuple, and every point of an orbit
   * carrying the same weight.
   *
   * This is the representation used by the rule families of Xiao and Gimbutas
   * @cite xiao2010numerical and of Witherden and Vincent
   * @cite witherden2015identification: an orbit is stored once, as a tuple and
   * a weight, and expanded on construction. It is compact, it makes the
   * symmetry structural rather than incidental, and it reduces the unknowns of
   * the moment equations from @f$ d \cdot n + n @f$ to a handful of orbit
   * parameters.
   *
   * The orbit classes conventionally named @f$ S_3 @f$, @f$ S_{21}(a) @f$ and
   * @f$ S_{111}(a,b) @f$ on the triangle, and @f$ S_4 @f$, @f$ S_{31}(a) @f$,
   * @f$ S_{22}(a) @f$, @f$ S_{211}(a,b) @f$, @f$ S_{1111}(a,b,c) @f$ on the
   * tetrahedron, are not enumerated here. They are exactly the multiplicity
   * patterns of the stored tuple, so expanding the distinct permutations
   * recovers each class and its cardinality without a separate taxonomy.
   *
   * @see expand
   */
  class SymmetricOrbit
  {
    public:
      /// @brief Barycentric coordinates; @f$ d + 1 @f$ entries for a
      /// @f$ d @f$-simplex.
      using Barycentric = std::vector<Real>;

      /// @brief Values along the tensor direction of a product element.
      ///
      /// Empty on a simplex. On the wedge, which is a triangle crossed with
      /// @f$ [0,1] @f$, the symmetry group is @f$ S_3 @f$ on the barycentric
      /// coordinates times the reflection @f$ z \mapsto 1-z @f$, so an orbit
      /// carries either the single mid-plane value @f$ \{1/2\} @f$ or the
      /// reflected pair @f$ \{c, 1-c\} @f$.
      using Tensor = std::vector<Real>;

      /**
       * @brief Constructs an orbit from a barycentric tuple and its weight.
       * @param barycentric Representative tuple, summing to one.
       * @param weight Weight carried by *each* point of the orbit, not by the
       * orbit as a whole.
       */
      SymmetricOrbit(Barycentric barycentric, Real weight)
        : m_barycentric(std::move(barycentric)),
          m_weight(weight)
      {}

      /**
       * @brief Constructs a product orbit on a tensor element.
       * @param barycentric Representative tuple on the base simplex.
       * @param tensor Values along the tensor direction.
       * @param weight Weight carried by each point of the orbit.
       *
       * The orbit is the Cartesian product of the barycentric permutations
       * with @p tensor, so its cardinality is the product of the two.
       */
      SymmetricOrbit(Barycentric barycentric, Tensor tensor, Real weight)
        : m_barycentric(std::move(barycentric)),
          m_tensor(std::move(tensor)),
          m_weight(weight)
      {}

      /// @brief The tensor-direction values; empty on a simplex.
      const Tensor& getTensor() const
      {
        return m_tensor;
      }

      /// @brief The representative barycentric tuple.
      const Barycentric& getBarycentric() const
      {
        return m_barycentric;
      }

      /// @brief The weight carried by each point of the orbit.
      Real getWeight() const
      {
        return m_weight;
      }

      /**
       * @brief The distinct permutations of the representative tuple.
       *
       * The cardinality is the multinomial coefficient of the multiplicity
       * pattern: 1, 3 or 6 on the triangle, and 1, 4, 6, 12 or 24 on the
       * tetrahedron, recovering the conventional orbit classes.
       */
      std::vector<Barycentric> expand() const
      {
        Barycentric sorted = m_barycentric;
        std::sort(sorted.begin(), sorted.end());
        std::vector<Barycentric> out;
        do
        {
          out.push_back(sorted);
        } while (std::next_permutation(sorted.begin(), sorted.end()));
        return out;
      }

      /// @brief Number of points contributed by this orbit: the number of
      /// distinct barycentric permutations times the number of tensor values.
      size_t getSize() const
      {
        return expand().size() * std::max<size_t>(m_tensor.size(), 1);
      }

      /**
       * @brief Maps a barycentric tuple to reference coordinates of @p g.
       *
       * The image is @f$ \sum_i \lambda_i v_i @f$ over the reference vertices
       * published by Geometry::Polytope::Traits, so the mapping follows the
       * library's reference element rather than restating it.
       */
      static Math::SpatialVector<Real> toReference(
        Geometry::Polytope::Type g, const Barycentric& barycentric)
      {
        const Geometry::Polytope::Traits traits(g);
        const size_t d = traits.getDimension();
        assert(barycentric.size() == traits.getVertexCount());
        Math::SpatialVector<Real> x;
        x.resize(static_cast<Eigen::Index>(d));
        x.setZero();
        for (size_t i = 0; i < barycentric.size(); ++i)
        {
          const auto& v = traits.getVertex(i);
          for (size_t k = 0; k < d; ++k)
          {
            x[static_cast<Eigen::Index>(k)] +=
              barycentric[i] * v[static_cast<Eigen::Index>(k)];
          }
        }
        return x;
      }

      /**
       * @brief Every point of the orbit, in reference coordinates of @p g.
       *
       * On a simplex this is the barycentric expansion mapped by toReference.
       * On a product element the barycentric part is mapped through the base
       * simplex --- the triangle, for the reference wedge of
       * Geometry::Polytope::Traits, which is the unit triangle crossed with
       * @f$ [0,1] @f$ --- and each tensor value is appended as the last
       * coordinate.
       */
      std::vector<Math::SpatialVector<Real>> expandPoints(
        Geometry::Polytope::Type g) const
      {
        std::vector<Math::SpatialVector<Real>> out;
        if (m_tensor.empty())
        {
          for (const auto& b : expand())
            out.push_back(toReference(g, b));
          return out;
        }

        const size_t d = Geometry::Polytope::Traits(g).getDimension();
        const auto base = baseSimplex(m_barycentric.size());
        for (const auto& b : expand())
        {
          const auto bx = toReference(base, b);
          for (const auto t : m_tensor)
          {
            Math::SpatialVector<Real> x;
            x.resize(static_cast<Eigen::Index>(d));
            for (Eigen::Index k = 0; k < bx.size(); ++k)
              x[k] = bx[k];
            x[static_cast<Eigen::Index>(d - 1)] = t;
            out.push_back(std::move(x));
          }
        }
        return out;
      }

      /**
       * @brief Orbit classes of the square cross-section of a pyramid.
       *
       * The reference pyramid of Geometry::Polytope::Traits has its apex over
       * the base corner @f$ (0,0,0) @f$, so it carries almost no symmetry: only
       * the swap @f$ x \leftrightarrow y @f$. It is however the image, under a
       * shear of unit determinant, of the pyramid whose apex sits over the base
       * centre, and that one has the full dihedral symmetry @f$ D_4 @f$ of the
       * square.
       *
       * Orbits are therefore built in centred coordinates
       * @f$ (u,v) \in [-\tfrac{1-z}{2}, \tfrac{1-z}{2}]^2 @f$ and sheared onto
       * the reference element by
       * @f$ x = u + \tfrac{1-z}{2},\ y = v + \tfrac{1-z}{2} @f$. Exactness is
       * affine invariant and the determinant is one, so weights transfer
       * unchanged. Symmetry in the reference coordinates is lost, which is
       * expected: it was a device for constructing the rule, not a property
       * required of it.
       */
      enum class PyramidClass
      {
        Centre,    ///< @f$ (0,0) @f$, one point.
        Axis,      ///< @f$ (\pm a, 0), (0, \pm a) @f$, four points.
        Diagonal,  ///< @f$ (\pm a, \pm a) @f$, four points.
        General    ///< the eight images of @f$ (a,b) @f$ under @f$ D_4 @f$.
      };

      /// @brief Number of points an orbit of @p c contributes.
      static size_t pyramidClassSize(PyramidClass c)
      {
        switch (c)
        {
          case PyramidClass::Centre:   return 1;
          case PyramidClass::Axis:     return 4;
          case PyramidClass::Diagonal: return 4;
          default:                     return 8;
        }
      }

      /// @brief Free parameters of @p c, besides its height and its weight.
      static size_t pyramidClassParameters(PyramidClass c)
      {
        switch (c)
        {
          case PyramidClass::Centre:   return 0;
          case PyramidClass::General:  return 2;
          default:                     return 1;
        }
      }

      /**
       * @brief The points of a pyramid orbit, in reference coordinates.
       * @param c Orbit class.
       * @param alpha First shape parameter, in @f$ (-1,1) @f$, scaled by the
       * half-width of the cross-section at height @p z.
       * @param beta Second shape parameter, used by PyramidClass::General.
       * @param z Height of the orbit.
       */
      static std::vector<Math::SpatialVector<Real>> expandPyramid(
        PyramidClass c, Real alpha, Real beta, Real z)
      {
        const Real half = (1 - z) / 2;
        std::vector<std::pair<Real, Real>> uv;
        switch (c)
        {
          case PyramidClass::Centre:
            uv = {{0, 0}};
            break;
          case PyramidClass::Axis:
            uv = {{alpha, 0}, {-alpha, 0}, {0, alpha}, {0, -alpha}};
            break;
          case PyramidClass::Diagonal:
            uv = {{alpha, alpha}, {-alpha, alpha},
                  {alpha, -alpha}, {-alpha, -alpha}};
            break;
          default:
            uv = {{alpha, beta}, {-alpha, beta}, {alpha, -beta}, {-alpha, -beta},
                  {beta, alpha}, {-beta, alpha}, {beta, -alpha}, {-beta, -alpha}};
            break;
        }
        std::vector<Math::SpatialVector<Real>> out;
        out.reserve(uv.size());
        for (const auto& [u, v] : uv)
        {
          Math::SpatialVector<Real> p;
          p.resize(3);
          p[0] = u * half + half;   // shear onto the corner-apex element
          p[1] = v * half + half;
          p[2] = z;
          out.push_back(std::move(p));
        }
        return out;
      }

      /// @brief The simplex carrying @p vertices barycentric coordinates.
      static Geometry::Polytope::Type baseSimplex(size_t vertices)
      {
        switch (vertices)
        {
          case 2:
            return Geometry::Polytope::Type::Segment;
          case 3:
            return Geometry::Polytope::Type::Triangle;
          default:
            return Geometry::Polytope::Type::Tetrahedron;
        }
      }

    private:
      Barycentric m_barycentric;
      Tensor m_tensor;
      Real m_weight;
  };
}

#endif
