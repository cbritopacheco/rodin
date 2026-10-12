/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_SWIFT_ADAPT_H
#define RODIN_ADAPTATION_SWIFT_ADAPT_H

#include "Rodin/Context/ForwardDecls.h"
#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Variational/P1.h"
#include "ForwardDecls.h"
#include "Problem.h"

namespace Rodin::Adaptation::SWIFT
{
  /**
   * @brief Mesh-fitting workflow selected by mesh context.
   *
   * Local meshes use Eigen storage. MPI meshes require a PETSc implementation;
   * that specialization is not yet available and construction is disabled.
   * No context is silently converted to a local mesh.
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref Adapt<Mesh, Context::Local> | In-place P1 fitting on a local mesh. |
   *
   * @tparam Mesh Mesh type, including its execution context.
   * @tparam ContextType Mesh execution context.
   */
  template <class Mesh, class ContextType>
  class Adapt
  {
    public:
      /// @brief Whether this context has an implemented fitting workflow.
      static constexpr bool isSupported = false;

      /**
       * @brief Unsupported contexts have no local fallback.
       * @param[in] mesh Mesh requiring an unimplemented backend specialization.
       */
      Adapt(Mesh& mesh) = delete;
  };

  /**
   * @brief Owns the local P1 fitting setup and displaces the supplied mesh.
   * @tparam Mesh Local mesh type.
   *
   * # Architecture
   *
   * The workflow owns a vector P1 space, trial and test functions, and a
   * @ref Problem. Each call to execute() starts from zero displacement on the
   * current mesh, computes a fit, and applies the accepted displacement only
   * when the sampled quality budget is satisfied. Valid best-effort results
   * are applied even if the geometric target was not reached.
   *
   * The mesh must outlive this object. Connectivity and classification must
   * remain unchanged; this operation neither remeshes nor classifies facets.
   * Only full-dimensional affine simplicial geometry is supported. Curved geometry requires
   * an explicitly chosen space and @ref Problem, followed by an appropriate
   * transformation update rather than vertex displacement.
   *
   * # Usage
   * @code{.cpp}
   * Adaptation::SWIFT::Adapt adapt(mesh);
   * Adaptation::SWIFT::Parameters parameters;
   * parameters.model.h = h; // Fixed background mesh scale.
   * adapt.setParameters(parameters).setInterfaceAttribute(interface);
   * const auto report = adapt.execute(phi, gradPhi); // Updates mesh in place.
   * @endcode
   *
   */
  template <class Mesh>
  class Adapt<Mesh, Context::Local> final
  {
    public:
      /// @brief This specialization uses the local Eigen variational backend.
      static constexpr bool isSupported = true;
      /// @brief Vector P1 space on the local mesh.
      using Space =
        Variational::P1<Math::SpatialVector<Real>, Geometry::Mesh<Context::Local>>;
      /// @brief Trial function owning displacement.
      using TrialFunction = decltype(Variational::TrialFunction(std::declval<Space&>()));
      /// @brief Test function on the displacement space.
      using TestFunction = decltype(Variational::TestFunction(std::declval<Space&>()));
      /// @brief Fitting problem using the owned trial and test functions.
      using ProblemType = Problem<TrialFunction, TestFunction>;

      /**
       * @brief Constructs the fitting setup without modifying the mesh.
       * @param[in,out] mesh Affine simplicial mesh to fit in place.
       */
      explicit Adapt(Mesh& mesh)
        : m_mesh(mesh),
          m_space(mesh, mesh.getDimension()),
          m_trial(m_space),
          m_test(m_space),
          m_problem(m_trial, m_test)
      {
        if (mesh.getSpaceDimension() != mesh.getDimension())
          Alert::Exception() << "SWIFT::Adapt requires a full-dimensional mesh."
                             << Alert::Raise;
        for (auto cell = mesh.getCell(); !cell.end(); ++cell)
        {
          const auto geometry = cell->getGeometry();
          if ((geometry != Geometry::Polytope::Type::Segment &&
                geometry != Geometry::Polytope::Type::Triangle &&
                geometry != Geometry::Polytope::Type::Tetrahedron) ||
            mesh.getPolytopeTransformation(mesh.getDimension(), cell->getIndex())
                .getOrder() != 1)
            Alert::Exception() << "SWIFT::Adapt requires affine simplicial geometry. "
                               << "Use SWIFT::Problem for curved geometry."
                               << Alert::Raise;
        }
      }

      /**
       * @brief Copying an owned adaptation workflow is disabled.
       * @param other Workflow that cannot be copied.
       */
      Adapt(const Adapt& other) = delete;
      /**
       * @brief Copy assignment of an owned adaptation workflow is disabled.
       * @param other Workflow that cannot be assigned.
       * @returns No value, since assignment is deleted.
       */
      Adapt& operator=(const Adapt& other) = delete;

      /**
       * @brief Configures the owned fitting problem.
       * @param[in] parameters Model, convergence and solver controls.
       * @returns This workflow.
       */
      Adapt& setParameters(const Parameters& parameters)
      {
        m_problem.setParameters(parameters);
        return *this;
      }

      /**
       * @brief Selects the preclassified interface facets.
       * @param[in] attribute Interface facet attribute.
       * @returns This workflow.
       */
      Adapt& setInterfaceAttribute(Geometry::Attribute attribute)
      {
        m_problem.setInterfaceAttribute(attribute);
        return *this;
      }

      /**
       * @brief Accesses metric, boundary and monitor configuration.
       * @returns The owned fitting problem.
       */
      ProblemType& getProblem()
      {
        return m_problem;
      }
      /**
       * @brief Inspects the owned fitting problem.
       * @returns The read-only owned problem.
       */
      const ProblemType& getProblem() const
      {
        return m_problem;
      }
      /**
       * @brief Accesses the trial function for additional variational terms.
       * @returns The owned trial function.
       */
      TrialFunction& getTrialFunction()
      {
        return m_trial;
      }
      /**
       * @brief Accesses the test function for additional variational terms.
       * @returns The owned test function.
       */
      TestFunction& getTestFunction()
      {
        return m_test;
      }
      /**
       * @brief Inspects the current controls.
       * @returns The owned problem parameters.
       */
      const Parameters& getParameters() const
      {
        return m_problem.getParameters();
      }
      /**
       * @brief Inspects the last fit.
       * @returns The last fitting report.
       */
      const Report& getReport() const
      {
        return m_problem.getReport();
      }

      /**
       * @brief Fits the current interface and applies a quality-valid displacement.
       * @param[in] phi Target level-set function evaluated in physical coordinates.
       * @param[in] gradient Target gradient evaluated in physical coordinates.
       * @returns Diagnostics relative to the mesh at the start of this call.
       */
      template <class Phi, class Gradient>
      Report execute(const Variational::RealFunctionBase<Phi>& phi,
        const Variational::VectorFunctionBase<Real, Gradient>& gradient)
      {
        m_trial.getSolution().getData().setZero();
        const auto report = m_problem.solve(phi, gradient);
        if (report.qualityBudgetSatisfied)
          m_mesh.get().displace(m_trial.getSolution());
        return report;
      }

    private:
      std::reference_wrapper<Mesh> m_mesh;
      Space m_space;
      TrialFunction m_trial;
      TestFunction m_test;
      ProblemType m_problem;
  };

  /**
   * @brief Selects the workflow from the supplied mesh context.
   * @param[in,out] mesh Mesh to fit in place.
   */
  template <class Mesh>
  Adapt(Mesh& mesh) -> Adapt<Mesh>;
}

#endif
