/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ATRIUM_SOLID_PROBLEMDATA_H
#define RODIN_ATRIUM_SOLID_PROBLEMDATA_H

#include <array>
#include <limits>
#include <string>
#include <set>
#include <vector>
#include <Eigen/Core>
#include <Rodin/Geometry.h>
#include <Rodin/Solid/Constitutive/DemirayGasser.h>
#include <Rodin/Variational.h>

namespace Rodin::Examples::Heart::AtriumContraction
{
  using Vec = Eigen::Matrix<Real, 3, 1>;
  using VectorSpace = Variational::P1<Math::SpatialVector<Real>>;
  using ScalarSpace = Variational::P1<Real>;
  using CellVectorSpace = Variational::P0<Math::SpatialVector<Real>>;
  using CellScalarSpace = Variational::P0<Real>;
  using VectorField = Variational::GridFunction<VectorSpace, Math::Vector<Real>>;
  using ScalarField = Variational::GridFunction<ScalarSpace, Math::Vector<Real>>;
  using CellVectorField = Variational::GridFunction<CellVectorSpace, Math::Vector<Real>>;
  using CellScalarField = Variational::GridFunction<CellScalarSpace, Math::Vector<Real>>;

  struct Options
  {
    std::string mesh =
      "resources/examples/Heart/Atrium/AtriumWall_Repaired_Wedges_1mm.mesh";
    std::string output = "output/atrium-contraction/Contraction";
    size_t steps = 24;
    size_t stopAfter = 0;
    Real timeStep = 0; // Zero uses the wave-duration time step.
    Geometry::Attribute baseLabel = 6;
    Vec appendageApex = Vec(46.39638367, 187.44307667, -205.03317767);
    Real apexRadius = 3;
    bool fibersOnly = false;
    Real activation = 20; // Positive electrical peak, in 1/s.
    Real relaxation = 10; // Magnitude of the resting electrical drive, in 1/s.
    Real speed = 200;
    Real width = 6;
    Rodin::Solid::DemirayGasser::Parameters passive;
    Real viscosity = 0.2; // Tissue viscosity, in kPa s.
    Real activeStiffness = 20; // Serial elastic stiffness, in kPa.
    Real activeDamping = 3; // Contractile dashpot, in kPa s.
    Real crossBridgeStiffness = 100; // Cross-bridge stiffness source, in kPa.
    Real contractility = 20; // Cross-bridge stress source, in kPa.
    Real destructionRate = 0.4; // Dimensionless strain-rate destruction factor.
    Real density = 1060; // Reference mass density in kg/m^3.
  };

  /**
   * @brief Parses the example configuration.
   * @param argc Argument count.
   * @param argv Argument values.
   * @returns Validated input configuration.
   */
  Options readOptions(int argc, char** argv);

  /**
   * @brief Validates a programmatic configuration.
   * @param cfg Material, wave, and time parameters.
   * @returns Validated configuration.
   */
  Options validateOptions(const Options& cfg);

  /**
   * @brief Loads the wedge wall and computes boundary connectivity.
   * @param cfg Input configuration.
   * @returns Wall mesh.
   */
  Geometry::LocalMesh loadWall(const Options& cfg);

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
   * Geometry inspection identifies the basal rim and extrusion normals. A scalar
   * Laplace problem supplies the fiber frame. The fitted plane defines activation;
   * the fields and their spaces are retained for constitutive input and output.
   * The supporting mesh and configuration must outlive this object.
   */
  struct ProblemData
  {
    Geometry::LocalMesh& mesh;
    const Options& cfg;
    ScalarSpace scalarSpace;
    CellVectorSpace fiberSpace;
    CellScalarSpace diagnosticSpace;
    ScalarField fiberCoordinate, waveCoordinate;
    CellVectorField fiberField, normalField, longitudinalField;
    CellScalarField fallbackField;
    std::vector<Vec> directions, normals, longitudinal;
    std::vector<Real> fallback;
    std::set<Index> baseAllVertices;
    Vec baseOrigin = Vec::Zero(), waveDirection = Vec::Zero();
    size_t rimFaces = 0, apexPatchFaces = 0, fallbackCells = 0;
    Real minimumReferenceJacobian = std::numeric_limits<Real>::infinity();
    Real coordinateResidual = 0, minimumProjection = 0, maximumProjection = 0;
    Real endTime = 0;

    /**
     * @brief Constructs the prescribed data on a loaded wall.
     * @param wall Supporting wall mesh.
     * @param options Prescribed input configuration.
     */
    ProblemData(Geometry::LocalMesh& wall, const Options& options);

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
     * @brief Bounds electrical drive over a wedge for an admissible Newton predictor.
     * @param cell Wedge whose reference vertices define the projection interval.
     * @param time Current physical time in seconds.
     * @returns Maximum drive in inverse seconds; residuals still use quadrature values.
     */
    Real maximumElectrical(const Geometry::Polytope& cell, Real time) const;

    /**
     * @brief Checks each deformed wedge over its entire reference cell.
     * @param displacement Current nodal displacement in mm.
     * @returns Minimum scaled determinant in mm, or no value for nonfinite geometry.
     */
    Optional<Real> minimumDeformedJacobian(const VectorField& displacement) const;

    private:
      /// @brief Validates the wall and identifies its basal plane and normals.
      void inspectGeometry();

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
