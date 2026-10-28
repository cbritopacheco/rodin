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
#include <Rodin/MPI.h>
#include <Rodin/PETSc.h>
#include <Rodin/Solid/Constitutive/DemirayHolzapfelGasserOgden.h>
#include <Rodin/Variational.h>

namespace Rodin::Examples::Heart::AtriumContraction
{
  using Vec = Eigen::Matrix<Real, 3, 1>;
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

  struct Options
  {
    std::string output = "output/atrium-contraction/Contraction";

    /// @brief Writes a standalone component CSV; the coordinator owns the run CSV.
    bool writeCSV = true;
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
    /// @brief Literature-based atrial passive coefficients; see the solid README.
    Rodin::Solid::DemirayHolzapfelGasserOgden::Parameters passive{
      .matrixScale = 2.92,
      .matrixExponent = 5.6,
      .fiberScale = 11.84,
      .fiberExponent = 17.95,
      .volumePenalty = 500
    };
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
   * @param meshPath Common parent mesh path selected by the command line.
   * @returns Validated input configuration.
   */
  Options readOptions(int argc, char** argv, std::string& meshPath);

  /**
   * @brief Validates a programmatic configuration.
   * @param cfg Material, wave, and time parameters.
   * @returns Validated configuration.
   */
  Options validateOptions(const Options& cfg);

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
   * The supporting distributed submesh and configuration must outlive this object.
   * Geometry reductions and the retained harmonic-coordinate solve are collective.
   */
  struct ProblemData
  {
    MeshType& mesh;
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
    ProblemData(MeshType& wall, const Options& options);

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

      /// @brief Collectively validates the wall and identifies its basal plane and normals.
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
