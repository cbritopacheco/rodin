/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_LEVELSETDISCRETIZERBASE_H
#define RODIN_GEOMETRY_LEVELSETDISCRETIZERBASE_H

#include <array>
#include <cassert>
#include <functional>
#include <variant>
#include <vector>

#include "Rodin/Geometry/Types.h"
#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Variational/P1/P1.h"

namespace Rodin::Geometry
{
  /**
   * @brief Base class for discretizers driven by a level-set.
   *
   * This class stores the common configuration/state for algorithms that
   * discretize a domain or interface described implicitly by a scalar field
   * (typically a level-set function) defined on a mesh.
   *
   * Attribute transfer is controlled by a per-dimension map (SplitMap) that
   * associates an input attribute to either:
   * - a Split policy, meaning "split geometrically and relabel by side", or
   * - a PreserveAttribute policy, meaning "split geometrically and keep the
   *   original input attribute".
   *
   * Geometry is handled entirely by derived classes. The policies stored here
   * affect only how attributes are transferred to the output entities/cells.
   *
   * @tparam Params Specialization parameters; this header provides a
   * specialization for:
   *   Variational::GridFunction<Variational::P1<Real, Mesh>, Data>.
   */
  template <class... Params>
  class LevelSetDiscretizerBase;

  /**
   * @brief Specialization for level-set discretization algorithms operating on a P1 grid function.
   *
   * @tparam Mesh Mesh type hosting the P1 space.
   * @tparam Data Underlying storage/backend of the GridFunction.
   */
  template <class Mesh, class Data>
  class LevelSetDiscretizerBase<Variational::GridFunction<Variational::P1<Real, Mesh>, Data>>
  {
    public:
      /// Mesh type associated with the underlying P1 finite element space.
      using MeshType = Mesh;

      /// Finite element space type (P1 scalar field).
      using FESType = Variational::P1<Real, Mesh>;

      /// Grid function type used as the input level-set.
      using GridFunctionType = Variational::GridFunction<FESType, Data>;

      /**
       * @brief Tag type representing the policy "preserve the original attribute".
       *
       * This policy means that if the derived discretizer geometrically splits an
       * entity/cell carrying a given input attribute, the descendants should keep
       * that original attribute instead of being relabeled by side.
       */
      struct PreserveAttributeT {};

      /// Singleton value used to specify the PreserveAttribute policy.
      static constexpr PreserveAttributeT PreserveAttribute{};

      /**
       * @brief Side-based relabeling policy.
       *
       * When splitting is performed, descendants may be relabeled according to
       * the sign of the level-set:
       * - @ref negative: attribute assigned to the negative side (phi < 0)
       * - @ref positive: attribute assigned to the positive side (phi > 0)
       *
       * The exact geometric meaning of "side" is defined by the derived class,
       * but the intended use is region/interface-aware relabeling.
       */
      struct Split
      {
        /// Attribute to assign to descendants classified on the negative side.
        Attribute negative;

        /// Attribute to assign to descendants classified on the positive side.
        Attribute positive;
      };

      /**
       * @brief Per-attribute attribute-transfer configuration map.
       *
       * For a fixed topological dimension @f$d@f$, this map associates an input
       * attribute to either:
       * - Split: split geometrically and relabel by side, or
       * - PreserveAttribute: split geometrically and keep the original attribute.
       *
       * Notes:
       * - The map is stored per dimension (see @ref m_split).
       * - If an attribute is absent from the map, the derived class may choose
       *   its own default behavior. A common choice is to preserve the original
       *   attribute.
       */
      using SplitMap = FlatMap<Attribute, std::variant<Split, PreserveAttributeT>>;

      /**
       * @brief Construct from a P1 grid function.
       *
       * Initializes one SplitMap per topological dimension of the mesh
       * (0..dim). The maps are initially empty.
       *
       * @param gf Input grid function (typically a level-set).
       */
      LevelSetDiscretizerBase(const GridFunctionType& gf)
        : m_gf(gf)
      {
        const size_t dim = gf.getFiniteElementSpace().getMesh().getDimension();
        m_split.resize(dim + 1);
        m_interface.fill({});
      }

      /**
       * @brief Set the attribute-transfer policy for an input attribute at a given dimension.
       *
       * @param d Topological dimension (0..mesh_dim).
       * @param attr Input attribute to configure.
       * @param value Either Split{neg,pos} or PreserveAttribute.
       * @return *this
       */
      LevelSetDiscretizerBase&
      setSplit(size_t d, const Attribute& attr, const std::variant<Split, PreserveAttributeT>& value)
      {
        assert(d < m_split.size());
        m_split[d][attr] = value;
        return *this;
      }

      /**
       * @brief Convenience wrapper to configure side-based relabeling.
       *
       * @param d Topological dimension (0..mesh_dim).
       * @param attr Input attribute to configure.
       * @param split Split policy (negative/positive relabel).
       * @return *this
       */
      LevelSetDiscretizerBase& split(size_t d, Attribute attr, const Split& split)
      {
        return this->setSplit(d, attr, split);
      }

      /**
       * @brief Convenience wrapper to configure attribute preservation.
       *
       * This means that descendants of entities/cells carrying @p attr are still
       * allowed to be split geometrically by the derived algorithm, but their
       * output attribute should remain equal to @p attr instead of being relabeled
       * by side.
       *
       * @param d Topological dimension (0..mesh_dim).
       * @param attr Input attribute to configure for preservation.
       * @return *this
       */
      LevelSetDiscretizerBase& preserve(size_t d, Attribute attr)
      {
        return this->setSplit(d, attr, PreserveAttribute);
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
       * @param d Topological dimension of the interface entities to mark.
       * @param attr Optional interface attribute.
       * @return *this
       */
      LevelSetDiscretizerBase& setInterface(size_t d, const Optional<Attribute>& attr)
      {
        assert(d < m_interface.size());
        m_interface[d] = attr;
        return *this;
      }

      /**
       * @brief Get the configured interface attribute for a given dimension.
       *
       * @param d Topological dimension.
       * @return Optional interface attribute (may be nullopt).
       */
      const Optional<Attribute>& getInterface(size_t d) const
      {
        assert(d < m_interface.size());
        return m_interface[d];
      }

      /**
       * @brief Get the input grid function driving the discretization.
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
       * Must be implemented by the derived class.
       *
       * @return Output mesh.
       */
      virtual MeshType discretize() const = 0;

    private:
      /// Reference to the input grid function (level-set).
      std::reference_wrapper<const GridFunctionType> m_gf;

      /// Per-dimension attribute-transfer configuration maps.
      std::vector<SplitMap> m_split;

      /// Optional attribute used to mark interface entities in the output.
      std::array<Optional<Attribute>, 4> m_interface;
  };
}

#endif
