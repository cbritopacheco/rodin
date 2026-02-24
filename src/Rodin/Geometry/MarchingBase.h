/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_MARCHINGBASE_H
#define RODIN_GEOMETRY_MARCHINGBASE_H

#include <functional>

#include "Rodin/Geometry/Types.h"
#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Variational/P1/P1.h"

namespace Rodin::Geometry
{
  /**
   * @brief Base class for marching-based discretizers driven by a level-set.
   *
   * This class provides common configuration/state for algorithms that
   * discretize a domain or interface described implicitly by a scalar field
   * (typically a level-set function) defined on a mesh.
   *
   * The input is a P1 grid function @f$\phi@f$ on a mesh @f$\mathcal{T}@f$.
   * A derived class implements a specific marching strategy (e.g. marching
   * tetrahedra) to produce an output mesh where:
   * - cells may be split according to the sign of @f$\phi@f$,
   * - regions may be relabeled depending on which side of the interface they
   *   belong to,
   * - interface entities (faces/edges) may be marked with a dedicated attribute.
   *
   * Splitting/relabeling is controlled by a per-dimension map (SplitMap) that
   * associates an input attribute to either:
   * - a Split policy (negative/positive relabel), or
   * - a NoSplit policy (prevent splitting where that attribute appears).
   *
   * @tparam Params Specialization parameters; this header provides a
   * specialization for:
   *   Variational::GridFunction<Variational::P1<Real, Mesh>, Data>.
   */
  template <class ... Params>
  class MarchingBase;

  /**
   * @brief Specialization for marching algorithms operating on a P1 grid function.
   *
   * The marching is driven by a scalar P1 grid function on a mesh. The derived
   * class is responsible for implementing discretize().
   *
   * @tparam Mesh Mesh type hosting the P1 space.
   * @tparam Data Underlying storage/backend of the GridFunction.
   */
  template <class Mesh, class Data>
  class MarchingBase<Variational::GridFunction<Variational::P1<Real, Mesh>, Data>>
  {
    public:
      /// Mesh type associated with the underlying P1 finite element space.
      using MeshType = Mesh;

      /// Finite element space type (P1 scalar field).
      using FESType = Variational::P1<Real, Mesh>;

      /// Grid function type driving the marching procedure.
      using GridFunctionType = Variational::GridFunction<FESType, Data>;

      /**
       * @brief Tag type representing a "do not split" policy.
       *
       * When an entity/cell is associated with an attribute configured as NoSplit
       * for a given topological dimension, the derived algorithm should avoid
       * performing a cut/split across that entity/cell in that dimension.
       */
      struct NoSplitT {};

      /// Singleton value used to specify the NoSplit policy.
      static constexpr NoSplitT NoSplit{};

      /**
       * @brief Split policy for an input attribute.
       *
       * When splitting is performed, the output entities/cells may be relabeled
       * based on the side of the level set:
       * - @ref negative: attribute assigned to the negative side (phi < 0)
       * - @ref positive: attribute assigned to the positive side (phi > 0)
       *
       * The exact semantics of "side" and how labels are applied are defined
       * by the derived class, but the intended use is region relabeling.
       */
      struct Split
      {
        /// Attribute to assign to elements/entities classified on the negative side.
        Attribute negative;

        /// Attribute to assign to elements/entities classified on the positive side.
        Attribute positive;
      };

      /**
       * @brief Per-attribute splitting configuration map.
       *
       * For a fixed topological dimension @f$d@f$, this map associates an input
       * attribute to either:
       * - Split: allow splitting and relabel output with the configured pair, or
       * - NoSplit: prevent splitting for that attribute in that dimension.
       *
       * Notes:
       * - The map is stored per dimension (see m_split).
       * - The default behavior when an attribute is absent from the map is left
       *   to the derived class (common choices: "split and keep label", or
       *   "do not split unless specified").
       */
      using SplitMap = FlatMap<Attribute, std::variant<Split, NoSplitT>>;

      /**
       * @brief Construct from a P1 grid function.
       *
       * Initializes one SplitMap per topological dimension of the mesh
       * (0..dim). The maps are initially empty.
       *
       * @param gf Input grid function (typically a level-set).
       */
      MarchingBase(const GridFunctionType& gf)
        : m_gf(gf)
      {
        m_split.resize(gf.getFiniteElementSpace().getMesh().getDimension() + 1);
      }

      /**
       * @brief Set the splitting policy for an attribute at a given dimension.
       *
       * @param d Topological dimension (0..mesh_dim).
       * @param attr Input attribute to configure.
       * @param value Either Split{neg,pos} or NoSplit.
       * @return *this
       */
      MarchingBase& setSplit(size_t d, const Attribute& attr, const std::variant<Split, NoSplitT>& value)
      {
        assert(d < m_split.size());
        m_split[d][attr] = value;
        return *this;
      }

      /**
       * @brief Convenience wrapper to configure a Split policy.
       *
       * @param d Topological dimension (0..mesh_dim).
       * @param attr Input attribute to configure.
       * @param split Split policy (negative/positive relabel).
       * @return *this
       */
      MarchingBase& split(size_t d, Attribute attr, const Split& split)
      {
        return this->setSplit(d, attr, split);
      }

      /**
       * @brief Convenience wrapper to configure a NoSplit policy.
       *
       * @param d Topological dimension (0..mesh_dim).
       * @param attr Input attribute to configure as "do not split".
       * @return *this
       */
      MarchingBase& noSplit(size_t d, Attribute attr)
      {
        return this->setSplit(d, attr, NoSplit);
      }

      /**
       * @brief Access the SplitMap for a given topological dimension.
       *
       * @param d Topological dimension (0..mesh_dim).
       * @return The corresponding SplitMap.
       */
      const SplitMap& getSplitMap(size_t d) const
      {
        assert(d < m_split.size());
        return m_split[d];
      }

      /**
       * @brief Set the attribute used to mark the extracted interface.
       *
       * The derived algorithm may mark interface entities (typically faces in 3D,
       * edges in 2D) with this attribute when they lie on, or separate, the
       * negative and positive side.
       *
       * If unset (nullopt), the derived algorithm should not assign a dedicated
       * interface attribute.
       *
       * @param attr Optional interface attribute.
       * @return *this
       */
      MarchingBase& setInterface(size_t d, const Optional<Attribute>& attr)
      {
        m_interface[d] = attr;
        return *this;
      }

      /**
       * @brief Get the configured interface attribute.
       *
       * @return Optional interface attribute (may be nullopt).
       */
      const Optional<Attribute>& getInterface(size_t d) const
      {
        return m_interface[d];
      }

      /**
       * @brief Get the input grid function driving the marching.
       *
       * @return Reference to the input GridFunction.
       */
      const GridFunctionType& getGridFunction() const
      {
        return m_gf.get();
      }

      /**
       * @brief Perform the discretization.
       *
       * The derived class implements a specific marching strategy and returns
       * a new mesh representing the discretized geometry (and possibly relabeled
       * regions/interface) induced by the input grid function.
       *
       * @return Output mesh generated by the marching algorithm.
       */
      virtual MeshType discretize() const = 0;

    private:
      /// Reference to the input grid function (level-set).
      std::reference_wrapper<const GridFunctionType> m_gf;

      /// Per-dimension split configuration maps (indexed by topological dimension).
      std::vector<SplitMap> m_split;

      /// Optional attribute used to mark interface entities in the output.
      std::array<Optional<Attribute>, 4> m_interface;
  };
}

#endif
