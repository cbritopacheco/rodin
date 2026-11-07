/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KELVINBALL_SEWEDOUTPUT_H
#define RODIN_KELVINBALL_SEWEDOUTPUT_H

#include <array>
#include <span>
#include <vector>

#include <Rodin/IO/XDMF.h>

#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

namespace KelvinBall
{
  using namespace Rodin;
  using namespace Rodin::Geometry;
  using namespace Rodin::Variational;

  /**
   * @brief Reconstructs the complete design from its fundamental chamber.
   *
   * The chamber is copied by the 24 proper rotations of the cube. Coincident
   * vertices are identified, and chamber fields are transported according to
   * their scalar, vector, or load-index transformation law.
   *
   * Architecture:
   * 1. construct the proper octahedral rotations;
   * 2. identify rotated copies of each chamber vertex;
   * 3. assemble the full-domain visualization mesh;
   * 4. sample fields lazily, averaging all chamber representatives.
   */
  class SewedOutput
  {
    public:
      using Mesh = Geometry::Mesh<Context::Local>;

      SewedOutput(const Mesh& chamber, const FlatSet<Attribute>& boundaryAttributes = {},
        Real tolerance = 1e-10);

      /// @brief A sewn mesh is owned explicitly and is not implicitly copied.
      SewedOutput(const SewedOutput& other) = delete;

      /// @brief Copy assignment of owned output geometry is disabled.
      SewedOutput& operator=(const SewedOutput& other) = delete;

      ~SewedOutput() = default;

      using VelocitySpace = P1<Math::SpatialVector<Real>, Mesh>;
      using PressureSpace = P1<Real, Mesh>;
      using Velocity = GridFunction<VelocitySpace, Math::Vector<Real>>;
      using Pressure = GridFunction<PressureSpace, Math::Vector<Real>>;

      /**
       * @brief Attaches the sewn Stokes families and their free-motion field.
       * @param grid Output grid, whose callbacks must be cleared before this object expires.
       * @param translations Chamber translational velocity family.
       * @param rotations Chamber rotational velocity family.
       * @param translationPressures Chamber translational pressure family.
       * @param rotationPressures Chamber rotational pressure family.
       * @param translation Free translational velocity.
       * @param angular Free angular velocity.
       *
       * Fields are sampled lazily from chamber coefficients, one attribute at
       * a time. No sewn finite element space or full state family is allocated.
       * This object and the chamber states must remain alive through the write;
       * clear the grid callbacks before their destruction.
       */
      void addStates(IO::XDMF::Grid& grid, const std::array<Velocity*, 3>& translations,
        const std::array<Velocity*, 3>& rotations,
        const std::array<Pressure*, 3>& translationPressures,
        const std::array<Pressure*, 3>& rotationPressures,
        const Math::SpatialVector<Real>& translation,
        const Math::SpatialVector<Real>& angular);

      const Mesh& getMesh() const;

      Mesh& getMesh();

      static const std::vector<Math::SpatialMatrix<Real>>& getCubeRotations();

      /**
       * @brief Averages a chamber scalar at a sewn vertex.
       * @param input Chamber field.
       * @param vertex Sewn vertex index.
       */
      template <class Input>
      Real getScalar(const Input& input, Index vertex) const
      {
        Real value = 0;
        const auto sources = getSources(vertex);
        for (const auto& source : sources)
        {
          const auto dofs = input.getFiniteElementSpace().getDOFs(0, source.vertex);
          value += input.getData()(dofs(0));
        }
        return value / sources.size();
      }

      /**
       * @brief Rotates and averages a chamber vector at a sewn vertex.
       * @param input Chamber field.
       * @param vertex Sewn vertex index.
       */
      template <class Input>
      Math::SpatialVector<Real> getVector(const Input& input, Index vertex) const
      {
        Math::SpatialVector<Real> value(3);
        value.setZero();
        const auto sources = getSources(vertex);
        for (const auto& source : sources)
        {
          const auto dofs = input.getFiniteElementSpace().getDOFs(0, source.vertex);
          Math::SpatialVector<Real> chamberValue(3);
          for (size_t component = 0; component < 3; ++component)
            chamberValue(component) = input.getData()(dofs(component));
          value += m_rotations[source.rotation] * chamberValue;
        }
        value /= sources.size();
        return value;
      }

      /**
       * @brief Rotates both the vector and load indices of a state family.
       * @param family Three chamber load responses.
       * @param load Sewn load index.
       * @param vertex Sewn vertex index.
       */
      template <class Family>
      Math::SpatialVector<Real> getVectorLoad(
        const Family& family, size_t load, Index vertex) const
      {
        Math::SpatialVector<Real> value(3);
        value.setZero();
        const auto sources = getSources(vertex);
        for (const auto& source : sources)
        {
          const auto& rotation = m_rotations[source.rotation];
          for (size_t sourceLoad = 0; sourceLoad < 3; ++sourceLoad)
          {
            const auto dofs =
              family[sourceLoad]->getFiniteElementSpace().getDOFs(0, source.vertex);
            Math::SpatialVector<Real> chamberValue(3);
            for (size_t component = 0; component < 3; ++component)
              chamberValue(component) = family[sourceLoad]->getData()(dofs(component));
            value += rotation(load, sourceLoad) * rotation * chamberValue;
          }
        }
        value /= sources.size();
        return value;
      }

      /**
       * @brief Rotates the load index and averages a pressure family.
       * @param family Three chamber pressures.
       * @param load Sewn load index.
       * @param vertex Sewn vertex index.
       */
      template <class Family>
      Real getScalarLoad(const Family& family, size_t load, Index vertex) const
      {
        Real value = 0;
        const auto sources = getSources(vertex);
        for (const auto& source : sources)
        {
          const auto& rotation = m_rotations[source.rotation];
          for (size_t sourceLoad = 0; sourceLoad < 3; ++sourceLoad)
          {
            const auto dofs =
              family[sourceLoad]->getFiniteElementSpace().getDOFs(0, source.vertex);
            value += rotation(load, sourceLoad) * family[sourceLoad]->getData()(dofs(0));
          }
        }
        return value / sources.size();
      }

    private:
      struct Source
      {
          size_t rotation;
          Index vertex;
      };

      /// @brief Borrowed representatives, in the original rotation/vertex order.
      std::span<const Source> getSources(Index vertex) const
      {
        return {m_sources.data() + m_sourceOffsets[vertex],
          m_sourceOffsets[vertex + 1] - m_sourceOffsets[vertex]};
      }

      Mesh m_mesh;
      std::vector<Source> m_sources;
      std::vector<Index> m_sourceOffsets;
      const std::vector<Math::SpatialMatrix<Real>>& m_rotations;
  };
}

#endif
