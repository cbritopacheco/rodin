/*
 * Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE or copy at https://www.boost.org/LICENSE_1_0.txt)
 */
#include "FSI.h"
#include "../Fluid/FlowIntegral.h"
#include "../Log.h"
#include <functional>
#include <unordered_map>
#include <iomanip>

namespace Rodin::Examples::Heart::Atrium
{
  using namespace Geometry;
  using namespace Variational;

  namespace
  {
    // Directional derivative coefficients of F.v + M:grad(v) + P.q + Q.grad(q).
    // Gradients use spatial coordinates; the direction is material mesh motion.
    struct ShapeCoefficients
    {
      Math::SpatialVector<Real> value = Math::SpatialVector<Real>(3),
                                pressureGradient = Math::SpatialVector<Real>(3);
      Math::SpatialMatrix<Real> gradient{3, 3};
      Math::SpatialVector<Real> supgResidual{3}, supgBase{3}, supgTransport{3};
      Real pressure = 0;
    };

    // Copy a field vector into one distributed block without assuming nodal ordering.
    inline void copyBlock(::Vec field, ::Vec mixed, PetscInt offset)
    {
      PetscInt begin = 0, end = 0;
      PetscErrorCode ierr = VecGetOwnershipRange(field, &begin, &end);
      assert(ierr == PETSC_SUCCESS);
      ::IS indices = nullptr;
      ierr = ISCreateStride(
        PetscObjectComm((PetscObject)mixed), end - begin, offset + begin, 1, &indices);
      assert(ierr == PETSC_SUCCESS);
      ::Vec block = nullptr;
      ierr = VecGetSubVector(mixed, indices, &block);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecCopy(field, block);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecRestoreSubVector(mixed, indices, &block);
      assert(ierr == PETSC_SUCCESS);
      ierr = ISDestroy(&indices);
      assert(ierr == PETSC_SUCCESS);
    }

    struct InterfaceNewtonContext
    {
      std::function<bool(::Vec, ::Vec)> residual;
      std::function<void(::Vec, ::Vec)> multiply;
      std::function<void(::Vec)> linearize;
      std::function<void()> observe;
      std::string failure;
      size_t evaluations = 0, products = 0, rejectedTrials = 0;
      bool root = false;

      static PetscErrorCode function(::SNES snes, ::Vec x, ::Vec f, void* data)
      {
        auto& c = *static_cast<InterfaceNewtonContext*>(data);
        try
        {
          if (!c.residual(x, f))
          {
            ++c.rejectedTrials;
            if (c.root)
            {
              Log("FSI/Newton trial", false, Log::Status::Retry)
                << "rejected: " << c.failure << Alert::Raise;
            }
            const PetscErrorCode ierr = SNESSetFunctionDomainError(snes);
            assert(ierr == PETSC_SUCCESS);
            return ierr;
          }
          return PETSC_SUCCESS;
        }
        catch (const std::exception& e)
        {
          c.failure = e.what();
          return PETSC_ERR_LIB;
        }
      }
      static PetscErrorCode jacobian(::SNES, ::Vec x, ::Mat, ::Mat, void* data)
      {
        auto& c = *static_cast<InterfaceNewtonContext*>(data);
        try
        {
          c.linearize(x);
          return PETSC_SUCCESS;
        }
        catch (const std::exception& e)
        {
          c.failure = e.what();
          return PETSC_ERR_LIB;
        }
      }
      static PetscErrorCode mult(::Mat matrix, ::Vec x, ::Vec y)
      {
        void* data = nullptr;
        PetscErrorCode ierr = MatShellGetContext(matrix, &data);
        assert(ierr == PETSC_SUCCESS);
        auto& c = *static_cast<InterfaceNewtonContext*>(data);
        try
        {
          ++c.products;
          c.multiply(x, y);
          return PETSC_SUCCESS;
        }
        catch (const std::exception& e)
        {
          c.failure = e.what();
          return PETSC_ERR_LIB;
        }
      }
      static PetscErrorCode monitor(
        ::SNES snes, PetscInt iteration, PetscReal norm, void* data)
      {
        const auto& c = *static_cast<InterfaceNewtonContext*>(data);
        PetscInt linear = 0;
        const PetscErrorCode ierr = SNESGetLinearSolveIterations(snes, &linear);
        assert(ierr == PETSC_SUCCESS);
        if (c.root)
        {
          Log("FSI/Newton") << "iteration=" << iteration
                            << " residual=" << std::setprecision(16) << norm
                            << " evaluations=" << c.evaluations
                            << " tangent_products=" << c.products
                            << " linear_iterations=" << linear << Alert::Raise;
        }
        if (c.observe)
          c.observe();
        return PETSC_SUCCESS;
      }
    };
  }

  void FSI::solveNewton()
  {
    constexpr PetscInt MaximumSolidIterations = 80, MaximumSolidEvaluations = 10000;
    constexpr PetscInt MaximumInterfaceIterations = 12,
                       MaximumInterfaceEvaluations = 1000;
    constexpr PetscInt MaximumLinearIterations = 200;
    constexpr Real LinearRelativeTolerance = 1e-8, LinearAbsoluteTolerance = 1e-14;
    constexpr Real LinearDivergenceTolerance = 1e5;
    const Real residualAtol = m_parameters.absoluteTolerance;
    m_solid.m_newton->setTolerances(m_parameters.solidAbsoluteTolerance, 0, 0,
      MaximumSolidIterations, MaximumSolidEvaluations);
    InterfaceNewtonContext context;
    size_t stateCacheHits = 0;
    double residualSeconds = 0, solidSensitivitySeconds = 0, aleSensitivitySeconds = 0;
    double fluidDerivativeSeconds = 0, fluidSensitivitySeconds = 0,
           tangentSetupSeconds = 0;
    size_t sensitivityActions = 0;
    context.root = m_fluid.isRoot();
    const MPI_Comm communicator =
      PetscObjectComm((PetscObject)m_solid.m_robinData.sync().getData());
    m_solid.m_displacement.sync();
    m_solid.m_extension.sync();
    m_solid.m_gamma.sync();
    m_solid.m_beta.sync();
    VectorField seedD(m_solid.m_displacement);
    CellScalarField seedE(m_solid.m_extension), seedGamma(m_solid.m_gamma),
      seedBeta(m_solid.m_beta);
    VectorField deltaG(m_solid.m_space), deltaT(m_solid.m_space), deltaD(m_solid.m_space);
    VectorField deltaImportedU(m_solid.m_space), deltaRate(m_solid.m_space);
    Fluid::VectorGridFunctionType direction(m_fluid.m_vh),
      deltaFluidTraction(m_fluid.m_vh);
    Fluid::VectorGridFunctionType deltaVelocity(m_fluid.m_vh), baseWall(m_fluid.m_vh);
    ::Vec xParts[2] = {nullptr, nullptr}, rParts[2] = {nullptr, nullptr};
    PetscErrorCode ierr;
    for (int i = 0; i < 2; ++i)
    {
      ierr = VecDuplicate(m_solid.m_robinData.getData(), &xParts[i]);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecDuplicate(m_solid.m_robinData.getData(), &rParts[i]);
      assert(ierr == PETSC_SUCCESS);
    }
    ierr = VecCopy(m_solid.m_robinData.getData(), xParts[0]);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecCopy(m_importedTraction.sync().getData(), xParts[1]);
    assert(ierr == PETSC_SUCCESS);
    ::Vec x = nullptr, residual = nullptr;
    ierr = VecCreateNest(communicator, 2, nullptr, xParts, &x);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecCreateNest(communicator, 2, nullptr, rParts, &residual);
    assert(ierr == PETSC_SUCCESS);
    ::Vec solidRhs = nullptr, solidDelta = nullptr, fluidRhs = nullptr,
          fluidDelta = nullptr;

    ::Vec evaluatedInput = nullptr, evaluatedResidual = nullptr;
    ierr = VecDuplicate(x, &evaluatedInput);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecDuplicate(residual, &evaluatedResidual);
    assert(ierr == PETSC_SUCCESS);
    bool evaluatedStateValid = false;
    const auto evaluate = [&](::Vec input, ::Vec output) {
      // Cache only the most recently evaluated point: the current bulk fields
      // and geometry must still belong to it. No tolerance-based comparison.
      PetscBool same = PETSC_FALSE;
      if (evaluatedStateValid)
      {
        ierr = VecEqual(input, evaluatedInput, &same);
        assert(ierr == PETSC_SUCCESS);
      }
      if (same)
      {
        ierr = VecCopy(evaluatedResidual, output);
        assert(ierr == PETSC_SUCCESS);
        ++stateCacheHits;
        return true;
      }
      evaluatedStateValid = false;
      ierr = VecCopy(input, evaluatedInput);
      assert(ierr == PETSC_SUCCESS);
      ++context.evaluations;
      const double residualStart = MPI_Wtime();
      ::Vec g = nullptr, t = nullptr, rg = nullptr, rt = nullptr;
      ierr = VecNestGetSubVec(input, 0, &g);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecNestGetSubVec(input, 1, &t);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecNestGetSubVec(output, 0, &rg);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecNestGetSubVec(output, 1, &rt);
      assert(ierr == PETSC_SUCCESS);
      m_solid.m_displacement = seedD;
      m_solid.m_extension = seedE;
      m_solid.m_gamma = seedGamma;
      m_solid.m_beta = seedBeta;
      m_solid.m_robinData.setData(g);
      m_importedTraction.setData(t);
      transfer(m_solidToFluid, m_importedTraction, m_fluid.m_wallTraction, 1e3);
      transfer(m_solidToFluid, m_fluid.m_wallTraction, m_importedTraction, 1e-3,
        SCATTER_REVERSE);
      m_solid.solveStep();
      moveFluidMesh();
      if (!m_fluid.solveFlow())
      {
        context.failure =
          "fluid velocity guard: speed=" + std::to_string(m_fluid.m_speed) +
          " m/s, limit=" + std::to_string(m_fluid.m_cfg.maxVelocity) + " m/s";
        residualSeconds += MPI_Wtime() - residualStart;
        return false;
      }
      updateTraction();
      ierr = VecWAXPY(rg, -1, m_solid.m_robinData.sync().getData(), g);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecWAXPY(rt, -1, m_trialTraction.sync().getData(), t);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecCopy(output, evaluatedResidual);
      assert(ierr == PETSC_SUCCESS);
      evaluatedStateValid = true;
      residualSeconds += MPI_Wtime() - residualStart;
      context.failure.clear();
      return true;
    };
    context.residual = evaluate;
    if (!evaluate(x, residual))
    {
      Alert::Exception() << "Invalid initial interface state: " << context.failure
                         << Alert::Raise;
    }
    ierr = VecDuplicate(m_solid.m_mechanics.getLinearSystem().getSolution(), &solidRhs);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecDuplicate(solidRhs, &solidDelta);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecDuplicate(m_fluid.m_flow.getLinearSystem().getSolution(), &fluidRhs);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecDuplicate(fluidRhs, &fluidDelta);
    assert(ierr == PETSC_SUCCESS);

    // Material shape derivative of the assembled SUPG/PSPG Oseen residual.
    // delta grad(f)=-grad(f)G, delta volume=tr(G) volume, G=grad(theta).
    struct CoefficientSample
    {
      Math::SpatialVector<Real> reference;
      ShapeCoefficients value;
    };
    std::unordered_map<Index, std::vector<CoefficientSample>> coefficientCache;
    const auto coefficients = [&](const Point& point) -> ShapeCoefficients {
      auto& samples = coefficientCache[point.getPolytope().getIndex()];
      const auto& reference = point.getReferenceCoordinates();
      for (const auto& sample : samples)
      {
        if ((sample.reference - reference).squaredNorm() == 0)
          return sample.value;
      }
      using V = Math::SpatialVector<Real>;
      using M = Math::SpatialMatrix<Real>;
      const auto& cfg = m_fluid.m_cfg;
      const Real rho = cfg.rho, dt = m_fluid.m_run.timeStep;
      const V u = m_fluid.m_u.getSolution()(point), old = m_fluid.m_uOld(point);
      const V a = old - m_fluid.m_meshVelocity(point), theta = direction(point),
              da = -theta / dt;
      const M U = Jacobian(m_fluid.m_u.getSolution()).getValue(point);
      const M O = Jacobian(m_fluid.m_uOld).getValue(point);
      const M W = Jacobian(m_fluid.m_meshVelocity).getValue(point);
      const M G = Jacobian(direction).getValue(point);
      const Real dilation = G.trace();
      const M dU = -U * G, dO = -O * G, dW = G / dt - W * G;
      const M D = .5 * (U + U.transpose()), dD = .5 * (dU + dU.transpose());
      const M oldD = .5 * (O + O.transpose()), oldDotD = .5 * (dO + dO.transpose());
      const Real pressure = m_fluid.m_p.getSolution()(point);
      const V gradP = Grad(m_fluid.m_p.getSolution()).getValue(point);
      const Real ratio = m_fluid.m_massRatios.at(point.getPolytope().getIndex());
      const Real b = O.trace() - 2 * W.trace(), db = dO.trace() - 2 * dW.trace();
      const V F = (rho / dt) * (u - ratio * old) + rho * U * a + (.5 * rho * b) * u;
      const V dF = (rho / dt) * (ratio * dilation * old) + rho * (dU * a + U * da) +
        (.5 * rho * db) * u;
      const V rm = F + gradP, drm = dF - G.transpose() * gradP;
      const auto& cy = cfg.viscosity;
      const Real gamma = std::sqrt(
        cy.gammaRegularization * cy.gammaRegularization + 2 * oldD.squaredNorm());
      Real contraction = 0;
      for (Index i = 0; i < 3; ++i)
      {
        for (Index j = 0; j < 3; ++j)
          contraction += oldD(i, j) * oldDotD(i, j);
      }
      const Real dGamma = 2 * contraction / gamma;
      const Real mu = m_fluid.viscosityAt(point);
      const Real power = std::pow(cy.lambda * gamma, cy.yasuda);
      const Real dMu = (cy.mu0 - cy.muInf) * (cy.n - 1) * cy.lambda *
        std::pow(cy.lambda * gamma, cy.yasuda - 1) *
        std::pow(1 + power, (cy.n - 1) / cy.yasuda - 1) * dGamma;
      const Real h = m_fluid.cellSize(point), dh = h * dilation / 3;
      const Real speed = a.norm(), dSpeed = speed > 0 ? a.dot(da) / speed : 0;
      const Real tau1 = m_fluid.tau1At(point);
      const Real dTau1 = -tau1 * tau1 *
        (4 * dMu / (rho * h * h) - 8 * mu * dh / (rho * h * h * h) + 2 * dSpeed / h -
          2 * speed * dh / (h * h));
      const Real denominator = rho / dt + rho / tau1;
      const Real dBase = rho * dTau1 / (tau1 * tau1 * denominator * denominator);
      const Real ts = m_fluid.vmsTauAt(point), dTs = cfg.vmsScale * dBase;
      const Real tp = m_fluid.tauPAt(point), dTp = cfg.pspgScale * dBase;
      const Real tc = m_fluid.tauCAt(point);
      const Real dTc = cfg.gradDivScale * rho * .25 *
        (2 * h * dh / tau1 - h * h * dTau1 / (tau1 * tau1));
      const M identity = M::Identity(3, 3);
      const M momentum = 2 * mu * D - pressure * identity + tc * U.trace() * identity;
      const M dMomentum =
        2 * dMu * D + 2 * mu * dD + (dTc * U.trace() + tc * dU.trace()) * identity;
      const Real pressureValue = U.trace() + cfg.pressurePenalty * pressure;
      const V pressureGradient = tp * rm;
      ShapeCoefficients result;
      result.value = dF + dilation * F;
      result.gradient = dMomentum + dilation * momentum - momentum * G.transpose();
      result.supgResidual = rho * ((dTs + dilation * ts) * rm + ts * drm);
      result.supgBase = rho * ts * rm;
      result.supgTransport = da - G * a;
      result.pressure = dU.trace() + dilation * pressureValue;
      result.pressureGradient =
        dTp * rm + tp * drm + dilation * pressureGradient - G * pressureGradient;
      samples.push_back({reference, result});
      return result;
    };
    const VectorFunction vectorCoefficient(
      3, [&](const Point& point) { return coefficients(point).value; });
    const MatrixFunction matrixCoefficient(
      3, 3, [&](const Point& point) { return coefficients(point).gradient; });
    const RealFunction pressureCoefficient(
      [&](const Point& point) { return coefficients(point).pressure; });
    const VectorFunction gradientCoefficient(
      3, [&](const Point& point) { return coefficients(point).pressureGradient; });
    const VectorFunction supgResidual(
      3, [&](const Point& p) { return coefficients(p).supgResidual; });
    const VectorFunction supgBase(
      3, [&](const Point& p) { return coefficients(p).supgBase; });
    const VectorFunction supgTransport(
      3, [&](const Point& p) { return coefficients(p).supgTransport; });
    const auto normal = BoundaryNormal(m_fluid.m_mesh);
    const VectorFunction robinDerivative(
      3, [&](const Point& point) -> Math::SpatialVector<Real> {
        const Real inv =
          m_fluid.m_inverseInterfaceStretches.at(point.getPolytope().getIndex());
        return inv *
          (deltaFluidTraction(point) -
            m_robinCoefficient * direction(point) / m_fluid.m_run.timeStep);
      });
    const VectorFunction convectionDerivative(
      3, [&](const Point& point) -> Math::SpatialVector<Real> {
        const Math::SpatialVector<Real> n = normal.getValue(point),
                                        theta = direction(point);
        const Math::SpatialMatrix<Real> G = Jacobian(direction).getValue(point);
        const Real normalStretch = n.dot(G * n),
                   surfaceDilation = G.trace() - normalStretch;
        const Math::SpatialVector<Real> dn = -G.transpose() * n + normalStretch * n;
        const Math::SpatialVector<Real> a =
          m_fluid.m_uOld(point) - m_fluid.m_meshVelocity(point);
        const Real rate = (-theta / m_fluid.m_run.timeStep).dot(n) + a.dot(dn) +
          surfaceDilation * a.dot(n);
        return (-.5 * m_fluid.m_cfg.rho * rate) * m_fluid.m_u.getSolution()(point);
      });
    Variational::LinearForm<Fluid::VectorFESType, ::Vec> flowDerivative(m_fluid.m_v);
    Variational::LinearForm<Fluid::ScalarFESType, ::Vec> pressureDerivative(m_fluid.m_q);
    auto robinForm = BoundaryIntegral(robinDerivative, m_fluid.m_v);
    robinForm.setOrder(2);
    auto convectionForm = BoundaryIntegral(convectionDerivative, m_fluid.m_v);
    convectionForm.setOrder(3);
    flowDerivative = FlowIntegral(vectorCoefficient, m_fluid.m_v) +
      FlowIntegral(matrixCoefficient, Jacobian(m_fluid.m_v)) +
      FlowIntegral(supgResidual,
        Mult(Jacobian(m_fluid.m_v), m_fluid.m_uOld - m_fluid.m_meshVelocity)) +
      FlowIntegral(supgBase, Mult(Jacobian(m_fluid.m_v), supgTransport)) +
      robinForm.over(m_fluid.m_cfg.wall) + convectionForm.over(m_fluid.m_cfg.wall);
    pressureDerivative = FlowIntegral(pressureCoefficient, m_fluid.m_q) +
      FlowIntegral(gradientCoefficient, Grad(m_fluid.m_q));
    Variational::LinearForm<VectorSpace, ::Vec> solidLoad(m_solid.m_test);
    solidLoad = BoundaryIntegral(deltaG, m_solid.m_test).over(Attribute{4});
    const VectorFunction zero(
      3, [](const Point&) { return Math::SpatialVector<Real>{0, 0, 0}; });
    auto fixed = DirichletBC(m_solid.m_increment, zero);
    fixed.on(Attribute{1}, Attribute{2}, Attribute{3}, Attribute{5}, Attribute{6});
    fixed.assemble();

    context.multiply = [&](::Vec input, ::Vec output) {
      ++sensitivityActions;
      const double solidSensitivityStart = MPI_Wtime();
      ::Vec dg = nullptr, dt = nullptr, jg = nullptr, jt = nullptr;
      ierr = VecNestGetSubVec(input, 0, &dg);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecNestGetSubVec(input, 1, &dt);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecNestGetSubVec(output, 0, &jg);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecNestGetSubVec(output, 1, &jt);
      assert(ierr == PETSC_SUCCESS);
      deltaG.setData(dg);
      deltaT.setData(dt);
      solidLoad.assemble();
      std::visit(
        [&](const auto& dofs) {
          for (const auto& [index, value] : dofs)
          {
            const PetscInt i = static_cast<PetscInt>(index);
            const PetscScalar z = 0;
            ierr = VecSetValues(solidLoad.getVector(), 1, &i, &z, INSERT_VALUES);
            assert(ierr == PETSC_SUCCESS);
          }
        },
        fixed.getDOFs());
      ierr = VecAssemblyBegin(solidLoad.getVector());
      assert(ierr == PETSC_SUCCESS);
      ierr = VecAssemblyEnd(solidLoad.getVector());
      assert(ierr == PETSC_SUCCESS);
      ierr = VecSet(solidRhs, 0);
      assert(ierr == PETSC_SUCCESS);
      copyBlock(solidLoad.getVector(), solidRhs, m_solid.m_mechanics.getTestOffsets()[0]);
      ierr = VecSet(solidDelta, 0);
      assert(ierr == PETSC_SUCCESS);
      ierr = KSPSolve(m_solid.m_linearSolver.getHandle(), solidRhs, solidDelta);
      assert(ierr == PETSC_SUCCESS);
      ::KSPConvergedReason solidReason;
      ierr = KSPGetConvergedReason(m_solid.m_linearSolver.getHandle(), &solidReason);
      assert(ierr == PETSC_SUCCESS);
      if (solidReason <= 0)
        Alert::Exception() << "Analytic solid sensitivity solve failed" << Alert::Raise;
      deltaD.setData(solidDelta, m_solid.m_mechanics.getTrialOffsets()[0]);
      solidSensitivitySeconds += MPI_Wtime() - solidSensitivityStart;
      const double aleSensitivityStart = MPI_Wtime();

      baseWall.setData(m_fluid.m_wallDisplacement.sync().getData());
      std::vector<Math::SpatialPoint> coordinates(m_referenceVertices.size());
      for (auto vertex = m_fluid.m_mesh.getVertex(); vertex; ++vertex)
        coordinates[vertex->getIndex()] = vertex->getCoordinates();
      transfer(m_solidToFluid, deltaD, m_fluid.m_wallDisplacement, 1e-3);
      restoreReference();
      m_motionTrial.getSolution() = Math::SpatialVector<Real>{0, 0, 0};
      m_motion.assemble();
      m_motion.solve(m_motionSolver);
      ::KSPConvergedReason meshReason;
      ierr = KSPGetConvergedReason(m_motionSolver.getHandle(), &meshReason);
      assert(ierr == PETSC_SUCCESS);
      if (meshReason <= 0)
      {
        Alert::Exception() << "Analytic ALE sensitivity solve failed, reason "
                           << int(meshReason) << Alert::Raise;
      }
      direction.setData(m_motionTrial.getSolution().getData());
      m_motionTrial.getSolution().setData(m_fluid.m_meshDisplacement.sync().getData());
      m_fluid.m_wallDisplacement.setData(baseWall.getData());
      for (Index i = 0; i < coordinates.size(); ++i)
        m_fluid.m_mesh.setVertexCoordinates(i, coordinates[i]);
      aleSensitivitySeconds += MPI_Wtime() - aleSensitivityStart;
      const double fluidDerivativeStart = MPI_Wtime();
      transfer(m_solidToFluid, deltaT, deltaFluidTraction, 1e3);
      transfer(m_solidToFluid, deltaFluidTraction, deltaT, 1e-3, SCATTER_REVERSE);
      // All coefficients depend on the current direction. Never retain them
      // across MatMult calls or residual/geometry updates.
      coefficientCache.clear();
      flowDerivative.assemble();
      pressureDerivative.assemble();
      ierr = VecSet(fluidRhs, 0);
      assert(ierr == PETSC_SUCCESS);
      copyBlock(flowDerivative.getVector(), fluidRhs, m_fluid.m_flow.getTestOffsets()[0]);
      copyBlock(
        pressureDerivative.getVector(), fluidRhs, m_fluid.m_flow.getTestOffsets()[1]);
      ierr = VecScale(fluidRhs, -1);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecSet(fluidDelta, 0);
      assert(ierr == PETSC_SUCCESS);
      fluidDerivativeSeconds += MPI_Wtime() - fluidDerivativeStart;
      const double fluidSensitivityStart = MPI_Wtime();
      ierr = KSPSolve(m_fluid.m_flowKSP.getHandle(), fluidRhs, fluidDelta);
      assert(ierr == PETSC_SUCCESS);
      fluidSensitivitySeconds += MPI_Wtime() - fluidSensitivityStart;
      ::KSPConvergedReason fluidReason;
      ierr = KSPGetConvergedReason(m_fluid.m_flowKSP.getHandle(), &fluidReason);
      assert(ierr == PETSC_SUCCESS);
      if (fluidReason <= 0)
        Alert::Exception() << "Analytic fluid sensitivity solve failed" << Alert::Raise;
      deltaVelocity.setData(fluidDelta, m_fluid.m_flow.getTrialOffsets()[0]);
      transfer(m_solidToFluid, deltaVelocity, deltaImportedU, 1e3, SCATTER_REVERSE);
      transfer(
        m_solidToFluid, direction, deltaRate, 1e3 / m_solid.m_timeStep, SCATTER_REVERSE);
      // DH_T=dT+alpha(du-dv); DH_g=DH_T+alpha du. J=I-DH.
      VectorField ht(m_solid.m_space), hg(m_solid.m_space);
      ht = deltaT + m_solid.m_robinCoefficient * (deltaImportedU - deltaRate);
      hg = ht + m_solid.m_robinCoefficient * deltaImportedU;
      ierr = VecWAXPY(jg, -1, hg.sync().getData(), dg);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecWAXPY(jt, -1, ht.sync().getData(), dt);
      assert(ierr == PETSC_SUCCESS);
    };

    context.observe = [&]() {
      m_residualDifference = m_solid.m_displacement - m_previousDisplacementIteration;
      const Real change = m_residualDifference.norm(NORM_2);
      const Real state = m_solid.m_displacement.norm(NORM_2);
      m_diagnostics.displacementChange = state > 0 ? change / state : 1e-3 * change;
      m_previousDisplacementIteration.setData(m_solid.m_displacement.sync().getData());
    };
    const bool checkTangent = m_parameters.checkTangent;
    bool checked = false;
    context.linearize = [&](::Vec input) {
      if (!evaluate(input, residual))
      {
        Alert::Exception() << "Invalid interface linearization state: " << context.failure
                           << Alert::Raise;
      }
      const auto assembleTangent = [&]() {
        const double setupStart = MPI_Wtime();
        m_solid.m_mechanics.assemble(AssemblyTarget::LHS);
        ::Mat ks = m_solid.m_mechanics.getLinearSystem().getOperator();
        ierr = KSPSetOperators(m_solid.m_linearSolver.getHandle(), ks, ks);
        assert(ierr == PETSC_SUCCESS);
        ierr = KSPSetReusePreconditioner(m_solid.m_linearSolver.getHandle(), PETSC_FALSE);
        assert(ierr == PETSC_SUCCESS);
        ierr = KSPSetReusePreconditioner(m_fluid.m_flowKSP.getHandle(), PETSC_FALSE);
        assert(ierr == PETSC_SUCCESS);
        tangentSetupSeconds += MPI_Wtime() - setupStart;
      };
      assembleTangent();
      if (checkTangent && !checked)
      {
        checked = true;
        ::Vec v = nullptr, jv = nullptr, xh = nullptr, fh = nullptr, base = nullptr;
        ierr = VecDuplicate(input, &v);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecDuplicate(input, &jv);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecDuplicate(input, &xh);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecDuplicate(input, &fh);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecDuplicate(input, &base);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecCopy(residual, v);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecCopy(residual, base);
        assert(ierr == PETSC_SUCCESS);
        PetscReal norm = 0;
        ierr = VecNorm(v, NORM_2, &norm);
        assert(ierr == PETSC_SUCCESS);
        if (norm > 0)
        {
          m_solid.m_displacement.sync();
          m_fluid.m_u.getSolution().sync();
          VectorField taylorBaseD(m_solid.m_displacement);
          Fluid::VectorGridFunctionType taylorBaseU(m_fluid.m_u.getSolution());
          ::Vec baseY = nullptr, yCandidate = nullptr, equationError = nullptr;
          VecDuplicate(fluidDelta, &baseY);
          VecDuplicate(fluidDelta, &yCandidate);
          VecDuplicate(fluidDelta, &equationError);
          VecCopy(m_fluid.m_flow.getLinearSystem().getSolution(), baseY);
          Fluid::VectorGridFunctionType taylorBaseMesh(m_fluid.m_meshDisplacement);
          context.multiply(v, jv);
          for (const Real h : {.1, .05, .025, .0125})
          {
            ierr = VecWAXPY(xh, h, v, input);
            assert(ierr == PETSC_SUCCESS);
            if (!evaluate(xh, fh))
            {
              Alert::Exception()
                << "Invalid interface Taylor point: " << context.failure << Alert::Raise;
            }
            ierr = VecAXPY(fh, -1, base);
            assert(ierr == PETSC_SUCCESS);
            ierr = VecAXPY(fh, -h, jv);
            assert(ierr == PETSC_SUCCESS);
            ierr = VecNorm(fh, NORM_2, &norm);
            assert(ierr == PETSC_SUCCESS);
            VectorField errorD(m_solid.m_space);
            errorD = m_solid.m_displacement - taylorBaseD - h * deltaD;
            Fluid::VectorGridFunctionType errorU(m_fluid.m_vh);
            errorU = m_fluid.m_u.getSolution() - taylorBaseU - h * deltaVelocity;
            const Real remainderD = errorD.norm(NORM_2), remainderU = errorU.norm(NORM_2);
            Fluid::VectorGridFunctionType errorMesh(m_fluid.m_vh);
            errorMesh = m_fluid.m_meshDisplacement - taylorBaseMesh - h * direction;
            const Real remainderMesh = errorMesh.norm(NORM_2);
            VecWAXPY(yCandidate, h, fluidDelta, baseY);
            MatMult(
              m_fluid.m_flow.getLinearSystem().getOperator(), yCandidate, equationError);
            VecAXPY(equationError, -1, m_fluid.m_flow.getLinearSystem().getVector());
            PetscReal equationRemainder = 0;
            VecNorm(equationError, NORM_2, &equationRemainder);
            if (context.root)
            {
              std::cout << "TAYLOR_FLUID_EQUATION h=" << h
                        << " residual=" << equationRemainder << " mesh=" << remainderMesh
                        << std::endl;
            }
            if (context.root)
            {
              std::cout << "TAYLOR_PARTS h=" << h << " solid=" << remainderD
                        << " fluid=" << remainderU << std::endl;
            }
            if (context.root)
            {
              std::cout << "TAYLOR h=" << h << " remainder=" << std::setprecision(16)
                        << norm << std::endl;
            }
          }
          VecDestroy(&baseY);
          VecDestroy(&yCandidate);
          VecDestroy(&equationError);
        }
        for (::Vec* v : {&v, &jv, &xh, &fh, &base})
        {
          ierr = VecDestroy(v);
          assert(ierr == PETSC_SUCCESS);
        }
        if (!evaluate(input, residual))
        {
          Alert::Exception() << "Cannot restore interface Taylor state: "
                             << context.failure << Alert::Raise;
        }
        assembleTangent();
      }
    };

    ::SNES outer = nullptr;
    ::Mat jacobian = nullptr;
    ierr = SNESCreate(communicator, &outer);
    assert(ierr == PETSC_SUCCESS);
    ierr = SNESSetOptionsPrefix(outer, "atrium_fsi_interface_");
    assert(ierr == PETSC_SUCCESS);
    ierr = SNESSetFunction(outer, residual, InterfaceNewtonContext::function, &context);
    assert(ierr == PETSC_SUCCESS);
    PetscInt local = 0, global = 0;
    ierr = VecGetLocalSize(x, &local);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecGetSize(x, &global);
    assert(ierr == PETSC_SUCCESS);
    ierr =
      MatCreateShell(communicator, local, local, global, global, &context, &jacobian);
    assert(ierr == PETSC_SUCCESS);
    ierr = MatShellSetOperation(
      jacobian, MATOP_MULT, (void (*)(void))InterfaceNewtonContext::mult);
    assert(ierr == PETSC_SUCCESS);
    ierr = SNESSetJacobian(
      outer, jacobian, jacobian, InterfaceNewtonContext::jacobian, &context);
    assert(ierr == PETSC_SUCCESS);
    ierr = SNESSetType(outer, SNESNEWTONLS);
    assert(ierr == PETSC_SUCCESS);
    ierr = SNESSetTolerances(
      outer, residualAtol, 0, 0, MaximumInterfaceIterations, MaximumInterfaceEvaluations);
    assert(ierr == PETSC_SUCCESS);
    ::KSP interfaceKsp = nullptr;
    ::PC interfacePc = nullptr;
    ierr = SNESGetKSP(outer, &interfaceKsp);
    assert(ierr == PETSC_SUCCESS);
    ierr = KSPSetType(interfaceKsp, KSPGMRES);
    assert(ierr == PETSC_SUCCESS);
    ierr = KSPSetTolerances(interfaceKsp, LinearRelativeTolerance,
      LinearAbsoluteTolerance, LinearDivergenceTolerance, MaximumLinearIterations);
    assert(ierr == PETSC_SUCCESS);
    ierr = KSPGetPC(interfaceKsp, &interfacePc);
    assert(ierr == PETSC_SUCCESS);
    ierr = PCSetType(interfacePc, PCNONE);
    assert(ierr == PETSC_SUCCESS);
    ierr = SNESMonitorSet(outer, InterfaceNewtonContext::monitor, &context, nullptr);
    assert(ierr == PETSC_SUCCESS);
    ierr = SNESSetFromOptions(outer);
    assert(ierr == PETSC_SUCCESS);
    ierr = SNESSolve(outer, nullptr, x);
    if (ierr != PETSC_SUCCESS)
    {
      Alert::Exception() << "Analytic interface Newton: " << context.failure
                         << Alert::Raise;
    }
    ::SNESConvergedReason reason;
    PetscInt iterations = 0;
    ierr = SNESGetConvergedReason(outer, &reason);
    assert(ierr == PETSC_SUCCESS);
    ierr = SNESGetIterationNumber(outer, &iterations);
    assert(ierr == PETSC_SUCCESS);
    if (reason <= 0)
    {
      Alert::Exception() << "Analytic interface Newton failed, reason " << int(reason)
                         << Alert::Raise;
    }
    if (!evaluate(x, residual))
    {
      Alert::Exception() << "Invalid accepted interface state: " << context.failure
                         << Alert::Raise;
    }
    PetscReal norm = 0;
    ierr = VecNorm(residual, NORM_2, &norm);
    assert(ierr == PETSC_SUCCESS);
    m_diagnostics.iterations = iterations;
    m_diagnostics.converged = true;
    PetscReal update = 0;
    ierr = VecNorm(residual, NORM_INFINITY, &update);
    assert(ierr == PETSC_SUCCESS);
    m_diagnostics.tractionDataUpdate = 1e3 * update;
    m_importedTraction.setData(m_trialTraction.sync().getData());
    transfer(m_solidToFluid, m_importedTraction, m_fluid.m_wallTraction, 1e3);
    if (context.root)
    {
      Log("FSI/Newton") << "converged: iterations=" << iterations << " residual=" << norm
                        << " evaluations=" << context.evaluations
                        << " tangent_products=" << context.products << Alert::Raise;
    }
    if (context.root)
    {
      Log("FSI/Newton timing") << "cache_hits=" << stateCacheHits
                               << " actual_evaluations=" << context.evaluations
                               << " rejected_trials=" << context.rejectedTrials
                               << " sensitivity_actions=" << sensitivityActions
                               << " residual_seconds=" << residualSeconds
                               << " tangent_setup_seconds=" << tangentSetupSeconds
                               << " solid_sensitivity_seconds=" << solidSensitivitySeconds
                               << " ale_sensitivity_seconds=" << aleSensitivitySeconds
                               << " fluid_derivative_seconds=" << fluidDerivativeSeconds
                               << " fluid_sensitivity_seconds=" << fluidSensitivitySeconds
                               << Alert::Raise;
    }
    ierr = SNESDestroy(&outer);
    assert(ierr == PETSC_SUCCESS);
    ierr = MatDestroy(&jacobian);
    assert(ierr == PETSC_SUCCESS);
    for (::Vec* v : {&x, &residual, &solidRhs, &solidDelta, &fluidRhs, &fluidDelta,
           &evaluatedInput, &evaluatedResidual})
    {
      ierr = VecDestroy(v);
      assert(ierr == PETSC_SUCCESS);
    }
    for (int i = 0; i < 2; ++i)
    {
      ierr = VecDestroy(&xParts[i]);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecDestroy(&rParts[i]);
      assert(ierr == PETSC_SUCCESS);
    }
  }
}
