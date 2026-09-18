/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef KELVIN_BALL_SEWED_OUTPUT_H
#define KELVIN_BALL_SEWED_OUTPUT_H

#include <vector>

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
   * 3. assemble the conforming full-domain mesh;
   * 4. transfer fields by averaging all chamber representatives.
   */
  class SewedOutput
  {
    public:
      using Mesh = Geometry::Mesh<Context::Local>;

      SewedOutput(const Mesh& chamber,
        const FlatSet<Attribute>& boundaryAttributes = {}, Real tolerance = 1e-10);

      const Mesh& getMesh() const;

      Mesh& getMesh();

      const std::vector<Math::SpatialMatrix<Real>>& getRotations() const;

      static const std::vector<Math::SpatialMatrix<Real>>& getCubeRotations();

      template <class Output, class Input>
      void setScalar(Output& output, const Input& input) const
      {
        for (Index vertex = 0; vertex < m_mesh.getVertexCount(); ++vertex)
        {
          const auto targetDOFs = output.getFiniteElementSpace().getDOFs(0, vertex);
          Real value = 0;
          for (const auto& source : m_sources[vertex])
          {
            const auto sourceDOFs =
              input.getFiniteElementSpace().getDOFs(0, source.vertex);
            value += input.getData()(sourceDOFs(0));
          }
          output.getData()(targetDOFs(0)) = value / m_sources[vertex].size();
        }
      }

      template <class Output, class Input>
      void setVector(Output& output, const Input& input) const
      {
        for (Index vertex = 0; vertex < m_mesh.getVertexCount(); ++vertex)
        {
          const auto targetDOFs = output.getFiniteElementSpace().getDOFs(0, vertex);
          Math::SpatialVector<Real> value(3);
          value.setZero();
          for (const auto& source : m_sources[vertex])
          {
            const auto sourceDOFs =
              input.getFiniteElementSpace().getDOFs(0, source.vertex);
            Math::SpatialVector<Real> chamberValue(3);
            for (size_t component = 0; component < 3; ++component)
              chamberValue(component) = input.getData()(sourceDOFs(component));
            value += m_rotations[source.rotation] * chamberValue;
          }
          value /= m_sources[vertex].size();
          for (size_t component = 0; component < 3; ++component)
            output.getData()(targetDOFs(component)) = value(component);
        }
      }

      template <class Output, class Family>
      void setVectorLoad(Output& output, const Family& family, size_t load) const
      {
        for (Index vertex = 0; vertex < m_mesh.getVertexCount(); ++vertex)
        {
          const auto targetDOFs = output.getFiniteElementSpace().getDOFs(0, vertex);
          Math::SpatialVector<Real> value(3);
          value.setZero();
          for (const auto& source : m_sources[vertex])
          {
            const auto& rotation = m_rotations[source.rotation];
            for (size_t sourceLoad = 0; sourceLoad < 3; ++sourceLoad)
            {
              const auto sourceDOFs =
                family[sourceLoad]->getFiniteElementSpace().getDOFs(0, source.vertex);
              Math::SpatialVector<Real> chamberValue(3);
              for (size_t component = 0; component < 3; ++component)
                chamberValue(component) =
                  family[sourceLoad]->getData()(sourceDOFs(component));
              value += rotation(load, sourceLoad) * rotation * chamberValue;
            }
          }
          value /= m_sources[vertex].size();
          for (size_t component = 0; component < 3; ++component)
            output.getData()(targetDOFs(component)) = value(component);
        }
      }

      template <class Output, class Family>
      void setScalarLoad(Output& output, const Family& family, size_t load) const
      {
        for (Index vertex = 0; vertex < m_mesh.getVertexCount(); ++vertex)
        {
          const auto targetDOFs = output.getFiniteElementSpace().getDOFs(0, vertex);
          Real value = 0;
          for (const auto& source : m_sources[vertex])
          {
            const auto& rotation = m_rotations[source.rotation];
            for (size_t sourceLoad = 0; sourceLoad < 3; ++sourceLoad)
            {
              const auto sourceDOFs =
                family[sourceLoad]->getFiniteElementSpace().getDOFs(0, source.vertex);
              value += rotation(load, sourceLoad) *
                family[sourceLoad]->getData()(sourceDOFs(0));
            }
          }
          output.getData()(targetDOFs(0)) = value / m_sources[vertex].size();
        }
      }

    private:
      struct Source
      {
          size_t rotation;
          Index vertex;
      };

      Mesh m_mesh;
      std::vector<std::vector<Source>> m_sources;
      std::vector<Math::SpatialMatrix<Real>> m_rotations;
  };
}

#endif
