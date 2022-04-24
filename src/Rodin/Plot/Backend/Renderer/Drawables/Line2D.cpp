/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <vector>
#include <numeric>

#include <Magnum/Mesh.h>
#include <Magnum/Math/Color.h>
#include <Magnum/GL/Renderer.h>
#include <Magnum/Math/Vector2.h>
#include <Magnum/GL/DefaultFramebuffer.h>
#include <Magnum/SceneGraph/Camera.h>

#include <Corrade/Containers/ArrayViewStl.h>

#include "Rodin/Plot/Geometry/Line2D.h"

#include "Line2D.h"

namespace Rodin::Plot::Backend::Renderer::Drawables
{
  Line2D::Line2D(
      Object2D& object,
      DrawableGroup2D* group,
      const Eigen::ArrayX<float>& xData,
      const Eigen::ArrayX<float>& yData,
      const Magnum::Color3& color,
      std::variant<LineStyle, DashTuple> lineStyle,
      float lineWidth,
      unsigned int lineSmoothness
      )
    : BaseDrawable2D(object, group),
      m_xData(xData), m_yData(yData),
      m_color(color),
      m_lineStyle(lineStyle),
      m_lineWidth(lineWidth),
      m_lineSmoothness(lineSmoothness)
  {
    assert(xData.size() == yData.size());
  }

  Line2D& Line2D::setColor(const Magnum::Color4& color)
  {
    m_color = color;
    return *this;
  }

  Line2D& Line2D::setLineWidth(float lineWidth)
  {
    m_lineWidth = lineWidth;
    return *this;
  }

  void Line2D::computeMesh(Magnum::SceneGraph::Camera2D& camera)
  {
    using LineSegment2D = Geometry::LineSegment2D<float>;

    struct Vertex {
      Magnum::Math::Vector2<float>  position;
    };

    assert(m_xData.size() == m_yData.size());

    long vertexCount = m_xData.size();
    std::vector<Vertex> vertices;
    vertices.reserve(vertexCount);

    std::vector<unsigned int> index;
    index.reserve(vertexCount);

    unsigned int idx = 0;
    for (long i = 0; i < vertexCount - 1; i++)
    {
      LineSegment2D segment(
          {m_xData[i], m_yData[i]},
          {m_xData[i + 1], m_yData[i + 1]}
          );

      Magnum::Math::Vector2<float> normal = segment.direction().perpendicular().normalized();

      // Use the scaled normal to ensure same thickness regardless of scale
      Magnum::Math::Vector2<float> scaledNormal
        = m_lineWidth / 2.0f * (camera.projectionMatrix().inverted() * (
            Magnum::Math::Vector3<float>(
              Magnum::Math::Vector2<float>(normal), 0))).xy();
      scaledNormal.x() /= Magnum::GL::defaultFramebuffer.viewport().sizeX();
      scaledNormal.y() /= Magnum::GL::defaultFramebuffer.viewport().sizeY();


      // TODO: We just draw quads for now, until probably the Core::Geometry
      // module is done so we can devise a more elegant solution
      vertices.push_back(Vertex{
          Magnum::Math::Vector2<float>(segment.start()
              ) - scaledNormal
          });
      vertices.push_back(Vertex{
              Magnum::Math::Vector2<float>(segment.end()
                  ) - scaledNormal
          });
      vertices.push_back(Vertex{
              Magnum::Math::Vector2<float>(segment.end()
                  ) + scaledNormal
          });
      vertices.push_back(Vertex{
              Magnum::Math::Vector2<float>(segment.start()
                  ) + scaledNormal
          });
      index.insert(index.end(), {idx, idx + 1, idx + 2});
      index.insert(index.end(), {idx + 2, idx + 3, idx});
      idx += 4;
    }

    Magnum::GL::Mesh mesh;
    Magnum::GL::Buffer vertexBuffer;
    vertexBuffer.setData(vertices);

    Magnum::GL::Buffer indices;
    indices.setData(index);

    mesh.setPrimitive(Magnum::GL::MeshPrimitive::Triangles)
        .setCount(vertexBuffer.size())
        .addVertexBuffer(
            std::move(vertexBuffer),
            0,
            Magnum::Shaders::FlatGL2D::Position{})
        .setIndexBuffer(
            std::move(indices),
            0,
            Magnum::MeshIndexType::UnsignedInt);

    m_mesh  = std::move(mesh);
  }

  void Line2D::draw(
      const Magnum::Matrix3& transformationMatrix,
      Magnum::SceneGraph::Camera2D& camera)
  {
    computeMesh(camera);

    Magnum::GL::Renderer::enable(Magnum::GL::Renderer::Feature::Multisampling);
    Magnum::GL::Renderer::enable(Magnum::GL::Renderer::Feature::Blending);
    Magnum::GL::Renderer::setBlendFunction(
        Magnum::GL::Renderer::BlendFunction::SourceAlpha,
        Magnum::GL::Renderer::BlendFunction::OneMinusSourceAlpha);

    // Draw to the framebuffer
    m_flatShader.setColor(m_color)
                .setTransformationProjectionMatrix(
                  camera.projectionMatrix() * transformationMatrix)
                .draw(m_mesh);
  }
}
