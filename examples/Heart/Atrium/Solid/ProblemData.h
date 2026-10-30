/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ATRIUM_SOLID_PROBLEMDATA_H
#define RODIN_ATRIUM_SOLID_PROBLEMDATA_H

#include "Parameters.h"
#include <array>
#include <limits>
#include <string>
#include <set>
#include <vector>
#include <Eigen/Core>
#include <Rodin/Geometry.h>
#include <Rodin/MPI.h>
#include <Rodin/PETSc.h>
#include <Rodin/Solid/Constitutive/DemirayHolzapfelGasserOgden.h>
#include <Rodin/Variational.h>

namespace Rodin::Examples::Heart::Atrium
{
  using MeshType = Geometry::Mesh<Context::MPI>;
  using SubMeshType = Geometry::SubMesh<Context::MPI>;
  using VectorSpace = Variational::P1<Math::SpatialVector<Real>, MeshType>;
  using ScalarSpace = Variational::P1<Real, MeshType>;
  using CellVectorSpace = Variational::P0<Math::SpatialVector<Real>, MeshType>;
  using CellScalarSpace = Variational::P0<Real, MeshType>;
  using VectorField = PETSc::Variational::GridFunction<VectorSpace>;
  using ScalarField = PETSc::Variational::GridFunction<ScalarSpace>;
  using CellVectorField = PETSc::Variational::GridFunction<CellVectorSpace>;
  using CellScalarField = PETSc::Variational::GridFunction<CellScalarSpace>;

  /// @brief Scalar diagnostics from a completed solid state, in mm/kPa/s units.
  struct DiagnosticState
  {
    Real maxDisplacement, activationMax;
    size_t newtonIterations;
    Real residual, minimumJacobian, maxVelocity, maxAcceleration;
    Real minimumExtension, maximumExtension, maximumStress;
  };

  /**
   * @brief Validates a programmatic configuration.
   * @param cfg Material, wave, and time parameters.
   * @returns Validated configuration.
   */
  SolidParameters validateParameters(const SolidParameters& cfg);

  /**
   * @brief Computes the minimum six-node wedge determinant on its reference cell.
   * @param q Physical wedge vertex coordinates.
   * @returns Minimum Jacobian determinant.
   */
  Real minimumJacobian(const std::array<Vec, 6>& q);

  /**
   * @brief Prescribed geometry, fiber frame, and activation data for the solid model.
   *
   * ## Architecture
   *
   * Geometry inspection identifies the basal rim. Harmonic transmural and
   * landmark coordinates supply a pointwise wall-tangent fiber frame. The fitted plane defines activation;
   * the fields and their spaces are retained for constitutive input and output.
   * The supporting distributed submesh and configuration must outlive this object.
   * Geometry reductions and the retained harmonic-coordinate solve are collective.
   */
  struct ProblemData
  {
    MeshType& mesh;
    const SolidParameters& cfg;
    ScalarSpace scalarSpace;
    CellVectorSpace fiberSpace;
    CellScalarSpace diagnosticSpace;
    ScalarField fiberCoordinate, transmuralCoordinate, waveCoordinate;
    CellVectorField fiberField, normalField, longitudinalField;
    CellScalarField fallbackField;
    std::vector<Vec> directions, normals, longitudinal;
    std::vector<Real> fallback;
    std::set<Index> baseAllVertices;
    Vec baseOrigin = Vec::Zero(), waveDirection = Vec::Zero();
    size_t rimFaces = 0, apexPatchFaces = 0, fallbackCells = 0;
    Real minimumReferenceJacobian = std::numeric_limits<Real>::infinity();
    Real coordinateResidual = 0, transmuralResidual = 0;
    Real endocardialTangencyError = 0, epicardialTangencyError = 0;
    Real minimumProjection = 0, maximumProjection = 0;
    Real endTime = 0;

    /**
     * @brief Constructs the prescribed data on a loaded wall.
     * @param wall Supporting wall mesh.
     * @param options Prescribed input configuration.
     */
    ProblemData(MeshType& wall, const SolidParameters& options);

    /**
     * @brief Evaluates the material fiber direction.
     * @param point Material point.
     * @returns Reference fiber direction.
     */
    Math::SpatialVector<Real> fiberAt(const Geometry::Point& point) const;

    /**
     * @brief Evaluates the prescribed activation.
     * @param point Material point.
     * @param time Physical time in seconds.
     * @returns Dimensionless activation.
     */
    Real activationAt(const Geometry::Point& point, Real time) const;

    /**
     * @brief Evaluates the electrical source for the active internal equations.
     * @param point Reference material point.
     * @param time Physical time in seconds.
     * @returns Electrical drive in inverse seconds, negative during relaxation.
     */
    Real electricalAt(const Geometry::Point& point, Real time) const;

    /**
     * @brief Collectively checks each owned wedge over its entire reference cell.
     * @param displacement Current nodal displacement in mm.
     * @returns Minimum scaled determinant in mm, or no value for nonfinite geometry.
     */
    Optional<Real> minimumDeformedJacobian(const VectorField& displacement) const;

    private:
      using CoordinateTrial = PETSc::Variational::TrialFunction<ScalarField, ScalarSpace>;
      using CoordinateTest = PETSc::Variational::TestFunction<ScalarSpace>;
      using CoordinateProblem =
        Variational::Problem<PETSc::Math::LinearSystem, CoordinateTrial, CoordinateTest>;
      CoordinateTrial m_coordinate;
      CoordinateTest m_coordinateTest;
      CoordinateProblem m_laplace;
      Solver::KSP m_coordinateSolver;

      /// @brief Collectively validates the wall and identifies its basal plane.
      void inspectGeometry();

      /// @brief Pointwise orthonormal reference material frame.
      struct Frame
      {
        Vec normal, longitudinal, fiber;
        Real fallback = 0;
      };

      /**
       * @brief Evaluates the frame from both coordinate gradients.
       * @param point Reference material point in a wall cell.
       * @returns Unit normal, projected longitudinal direction, and fiber direction.
       */
      Frame frameAt(const Geometry::Point& point) const;

      /// @brief Solves the transmural coordinate with endocardial and epicardial data.
      void generateTransmuralCoordinate();

      /// @brief Constructs the harmonic fiber frame and its diagnostic fields.
      void generateFibers();

      /**
       * @brief Projects a reference point onto the activation direction.
       * @param point Reference material point.
       * @returns Signed distance from the fitted basal plane in mm.
       */
      Real waveCoordinateAt(const Geometry::Point& point) const;
  };
}

#endif
