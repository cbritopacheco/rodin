/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_H1_FEKETE_H
#define RODIN_VARIATIONAL_H1_FEKETE_H

#include <array>
#include <cstddef>

#include "Rodin/Types.h"
#include "Rodin/Math/Vector.h"

#include "WarpBlend.h"

namespace Rodin::Variational
{
  /**
   * @brief Cached Fekete-type nodes on the reference triangle for degree K.
   *
   * The nodes are constructed as:
   *  - start from equispaced barycentric nodes on the reference triangle
   *    with vertices (0,0), (1,0), (0,1),
   *  - apply the warp–blend algorithm to move them towards Fekete positions.
   *
   * For a given polynomial degree K, the total number of nodes is
   * \f$ (K+1)(K+2)/2 \f$. The nodes are computed once per template
   * instantiation and cached in a static array.
   *
   * @tparam K Polynomial degree on the triangle.
   */
  template <size_t K>
  class FeketeTriangle
  {
    public:
      /// Total number of nodes: (K+1)(K+2)/2.
      static constexpr size_t Count = (K + 1) * (K + 2) / 2;

      /// Return cached nodes as a std::array.
      static const std::array<Math::SpatialPoint, Count>& getNodes()
      {
        static const std::array<Math::SpatialPoint, Count> s_nodes = compute();
        return s_nodes;
      }

    private:
      /// Build equispaced nodes and apply warp–blend once.
      static constexpr std::array<Math::SpatialPoint, Count> compute()
      {
        std::array<Math::SpatialPoint, Count> nodes{};
        size_t idx = 0;

        if constexpr (K == 0)
        {
          // Single vertex node; choose (0,0) by convention.
          nodes[0] = Math::SpatialPoint{{0.0, 0.0}};
          return nodes;
        }

        // Equispaced nodes on the reference triangle (0,0)-(1,0)-(0,1)
        for (size_t j = 0; j <= K; ++j)
        {
          for (size_t i = 0; i <= K - j; ++i, ++idx)
          {
            const Real s = static_cast<Real>(i) / static_cast<Real>(K);
            const Real t = static_cast<Real>(j) / static_cast<Real>(K);
            nodes[idx] = Math::SpatialPoint{{s, t}};
          }
        }

        // Move equispaced nodes toward Fekete positions using warp–blend
        WarpBlendTriangle<K>::template apply<Count>(nodes);

        return nodes;
      }
  };

  /**
   * @brief Cached Fekete-type nodes on the reference tetrahedron for degree K.
   *
   * Reference tetrahedron:
   *   (0,0,0), (1,0,0), (0,1,0), (0,0,1).
   *
   * Construction:
   *  - start from equispaced barycentric nodes (integer lattice with
   *    i + j + k <= K, normalized by K),
   *  - apply warp–blend algorithm (WarpBlendTetrahedron<K>) to move them
   *    towards Fekete-type positions.
   *
   * For a given polynomial degree K, the number of nodes is
   * @f$ (K + 1)(K + 2)(K + 3) / 6 @f$.
   *
   * Nodes are computed once per template instantiation and cached.
   *
   * @tparam K Polynomial degree on the tetrahedron.
   */
  template <size_t K>
  class FeketeTetrahedron
  {
    public:
      /// Total number of nodes: (K+1)(K+2)(K+3)/6.
      static constexpr size_t Count = (K + 1) * (K + 2) * (K + 3) / 6;

      /// Return cached nodes as a std::array.
      static const std::array<Math::SpatialPoint, Count>& getNodes()
      {
        static const std::array<Math::SpatialPoint, Count> s_nodes = compute();
        return s_nodes;
      }

    private:
      /// Build equispaced nodes and apply warp–blend once.
      static constexpr std::array<Math::SpatialPoint, Count> compute()
      {
        std::array<Math::SpatialPoint, Count> nodes{};
        size_t idx = 0;

        if constexpr (K == 0)
        {
          // Single vertex node; choose (0,0,0) by convention.
          nodes[0] = Math::SpatialPoint{{0.0, 0.0, 0.0}};
          return nodes;
        }

        // Equispaced nodes on the reference tetrahedron:
        // vertices (0,0,0), (1,0,0), (0,1,0), (0,0,1)
        for (size_t k = 0; k <= K; ++k)
        {
          for (size_t j = 0; j <= K - k; ++j)
          {
            for (size_t i = 0; i <= K - j - k; ++i, ++idx)
            {
              const Real r = static_cast<Real>(i) / static_cast<Real>(K);
              const Real s = static_cast<Real>(j) / static_cast<Real>(K);
              const Real t = static_cast<Real>(k) / static_cast<Real>(K);
              nodes[idx] = Math::SpatialPoint{{r, s, t}};
            }
          }
        }

        // Move equispaced nodes toward Fekete positions using warp–blend
        WarpBlendTetrahedron<K>::template apply<Count>(nodes);

        return nodes;
      }
  };
}

#endif
