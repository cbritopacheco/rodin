/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ATRIUM_CAVITYVOLUME_H
#define RODIN_ATRIUM_CAVITYVOLUME_H

#include <cmath>
#include <vector>
#include <boost/mpi/collectives.hpp>
#include <Rodin/Alert.h>
#include <Rodin/MPI/Geometry/Mesh.h>
#include <Rodin/MPI/Geometry/SubMesh.h>
#include <Rodin/Variational/FaceNormal.h>
#include <Rodin/Variational/Jacobian.h>
#include <Rodin/Variational/P1/P1Element.h>

namespace Rodin::Examples::Heart
{
  /**
   * @brief Cavity volume bounded by a moving inner wall and fixed opening caps.
   *
   * ## Architecture
   *
   * The common mesh's fluid region supplies the reference volume. Only owned
   * solid interface faces contribute the boundary-volume change. Parent maps
   * identify the interface independently of face labels. Fixed inlet/outlet
   * caps contribute no change. No displacement extension into the fluid is used.
   *
   * For the cavity-outward normal, the divergence theorem gives
   * @f$ V(u)=V_0+\frac13\int_{\Gamma_0}
   * [(X+u)\cdot\operatorname{cof}(I+\nabla u)N-X\cdot N]\,dS @f$.
   * Face quadrature is embedded into the incident solid reference cell before
   * evaluating displacement and its gradient. Tangential derivatives realize the
   * deformed area vector without extending displacement into the fluid. Volumes
   * use the mesh length unit cubed.
   */
  class CavityVolume
  {
    public:
      using MeshType = Geometry::Mesh<Context::MPI>;

      /**
       * @brief Collectively measures the common reference cavity.
       * @param parent Common mesh, which must outlive this object.
       * @param fluidMaterial Label identifying the filled reference cavity.
       */
      CavityVolume(const MeshType& parent, Geometry::Attribute fluidMaterial);

      /// @brief Returns the retained reference volume in source length units cubed.
      Real getReferenceVolume() const noexcept;

      /**
       * @brief Collectively integrates owned cells, excluding overlap replicas.
       * @param mesh Distributed volume mesh.
       * @param material Optional cell label to select; empty selects all cells.
       * @returns Global volume in mesh length units cubed.
       */
      static Real measure(
        const MeshType& mesh, Optional<Geometry::Attribute> material = {});

      /**
       * @brief Collectively measures the cavity enclosed by the deformed solid wall.
       * @tparam Field Vector displacement field on a submesh of the common parent.
       * @param displacement Displacement in the parent's length unit, fixed on cap rims.
       * @returns Global enclosed volume in source length units cubed.
       */
      template <class Field>
      Real evaluate(const Field& displacement) const
      {
        const auto& mesh = displacement.getFiniteElementSpace().getMesh();
        if (!mesh.isSubMesh() || &mesh.asSubMesh().getParent() != &m_parent)
        {
          Alert::Exception()
            << "Cavity displacement requires a submesh of the common parent"
            << Alert::Raise;
        }
        const auto& map = mesh.asSubMesh().getPolytopeMap(2).left;
        Variational::FaceNormal normal(mesh);
        auto gradient = Variational::Jacobian(displacement);
        Real change = 0;
        size_t faces = 0;
        for (auto face = mesh.getFace(); face; ++face)
        {
          const Index index = face->getIndex();
          if (!mesh.getShard().isOwned(2, index) || !mesh.getShard().isBoundary(index))
            continue;
          bool interface = false;
          for (Index cell :
            m_parent.getConnectivity().getIncidence({2, 3}, map.at(index)))
            interface = interface || m_parent.getAttribute(3, cell) == m_fluidMaterial;
          if (!interface)
            continue;
          ++faces;
          const auto& incidence = mesh.getConnectivity().getIncidence({2, 3}, index);
          assert(incidence.size() == 1);
          const auto cell = mesh.getCell(incidence[0]);
          const Geometry::Polytope::Traits cellTraits(cell->getGeometry());
          const Variational::P1Element<Real> faceElement(face->getGeometry());
          std::vector<Math::SpatialPoint> referenceVertices;
          for (Index vertex : face->getVertices())
          {
            size_t local = 0;
            while (
              local < cell->getVertices().size() && cell->getVertices()[local] != vertex)
              ++local;
            assert(local < cell->getVertices().size());
            referenceVertices.push_back(cellTraits.getVertex(local));
          }
          const auto& formula =
            QF::QuadratureFormula::get(QuadratureOrder, face->getGeometry());
          const auto& quadrature = face->getQuadrature(formula);
          for (size_t q = 0; q < quadrature.getSize(); ++q)
          {
            const auto& point = quadrature.getPoint(q);
            Math::SpatialPoint reference(3);
            reference.setZero();
            for (size_t j = 0; j < referenceVertices.size(); ++j)
            {
              reference += faceElement.getBasis(j)(point.getReferenceCoordinates()) *
                referenceVertices[j];
            }
            const Geometry::Point cellPoint(*cell, reference);
            const auto& jacobian = point.getJacobian();
            const Math::SpatialMatrix<Real> deformed =
              jacobian + gradient.getValue(cellPoint) * jacobian;
            const Math::SpatialVector<Real> first = jacobian.col(0),
                                            second = jacobian.col(1);
            const Math::SpatialVector<Real> currentFirst = deformed.col(0),
                                            currentSecond = deformed.col(1);
            const Math::SpatialVector<Real> referenceArea = first.cross(second);
            const Math::SpatialVector<Real> currentArea =
              currentFirst.cross(currentSecond);
            const Real orientation =
              referenceArea.dot(normal.getValue(point)) > 0 ? -1 : 1;
            const Math::SpatialVector<Real> position = point.getCoordinates();
            const Math::SpatialVector<Real> current = position + displacement(cellPoint);
            change += formula.getWeight(q) * orientation *
              (current.dot(currentArea) - position.dot(referenceArea)) / 3;
          }
        }
        const auto& comm = m_parent.getContext().getCommunicator();
        const Real volume =
          m_referenceVolume + boost::mpi::all_reduce(comm, change, std::plus<Real>());
        if (boost::mpi::all_reduce(comm, faces, std::plus<size_t>()) == 0)
        {
          Alert::Exception() << "The solid has no interface with the reference cavity"
                             << Alert::Raise;
        }
        if (!std::isfinite(volume) || volume <= 0)
        {
          Alert::Exception() << "Invalid enclosed atrial volume " << volume
                             << Alert::Raise;
        }
        return volume;
      }

    private:
      // Exact for quadratic displacement on affine triangular interface faces.
      static constexpr size_t QuadratureOrder = 4;
      const MeshType& m_parent;
      Geometry::Attribute m_fluidMaterial;
      Real m_referenceVolume;
  };
}
#endif
