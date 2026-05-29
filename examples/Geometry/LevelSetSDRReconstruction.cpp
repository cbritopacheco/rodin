/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/Math.h>
#include <Rodin/Solver/HouseholderQR.h>
#include <Rodin/Solver/NewtonSolver.h>
#include <Rodin/Variational.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Rodin;
using namespace Rodin::Geometry;

namespace
{
  using Vec2 = Math::SpatialVector<Real>;

  Vec2 vec2(Real x = 0, Real y = 0)
  {
    Vec2 out(2);
    out(0) = x;
    out(1) = y;
    return out;
  }

  Real dot(const Vec2& a, const Vec2& b)
  {
    return a.dot(b);
  }

  Real cross(const Vec2& a, const Vec2& b)
  {
    return a(0) * b(1) - a(1) * b(0);
  }

  Real norm(const Vec2& a)
  {
    return a.norm();
  }

  Vec2 fromPoint(const Math::SpatialPoint& p)
  {
    return p;
  }

  Math::SpatialPoint toPoint(const Vec2& p)
  {
    return p;
  }

  enum class PhaseMomentMap
  {
    Tanh,
    Sign,
    Clamped
  };

  Real applyPhaseMomentMap(Real phi, Real epsilon, PhaseMomentMap map)
  {
    switch (map)
    {
      case PhaseMomentMap::Tanh:
        return std::tanh(phi / epsilon);
      case PhaseMomentMap::Sign:
        return phi < 0 ? -1 : 1;
      case PhaseMomentMap::Clamped:
        return std::clamp(phi / epsilon, Real(-1), Real(1));
    }
    return std::tanh(phi / epsilon);
  }

  struct CircleLevelSet
  {
    Real cx = 0.51;
    Real cy = 0.48;
    Real radius = 0.31;

    Real phi(const Vec2& p) const
    {
      const Real dx = p(0) - cx;
      const Real dy = p(1) - cy;
      return std::sqrt(dx * dx + dy * dy) - radius;
    }

    Vec2 grad(const Vec2& p) const
    {
      const Real dx = p(0) - cx;
      const Real dy = p(1) - cy;
      const Real r = std::max(std::sqrt(dx * dx + dy * dy), Real(1e-14));
      return vec2(dx / r, dy / r);
    }

    Math::SpatialMatrix<Real> hess(const Vec2& p) const
    {
      const Real dx = p(0) - cx;
      const Real dy = p(1) - cy;
      const Real r = std::max(std::sqrt(dx * dx + dy * dy), Real(1e-14));
      Math::SpatialMatrix<Real> H(2, 2);
      H(0, 0) = 1 / r - dx * dx / (r * r * r);
      H(0, 1) = -dx * dy / (r * r * r);
      H(1, 0) = H(0, 1);
      H(1, 1) = 1 / r - dy * dy / (r * r * r);
      return H;
    }
  };

  constexpr std::array<std::array<Real, 3>, 3> TriangleBarycentricQuadrature = {{
    {{ Real(2) / 3, Real(1) / 6, Real(1) / 6 }},
    {{ Real(1) / 6, Real(2) / 3, Real(1) / 6 }},
    {{ Real(1) / 6, Real(1) / 6, Real(2) / 3 }}
  }};

  struct CellData
  {
    Index index = 0;
    std::array<Index, 3> vertices = {{ 0, 0, 0 }};
    std::array<Vec2, 3> x;
    Real area = 0;
    Real detScale = 0;
    Real moment = 0;
  };

  Vec2 interpolate(
      const std::array<Vec2, 3>& values,
      const std::array<Real, 3>& bary)
  {
    return bary[0] * values[0] + bary[1] * values[1] + bary[2] * values[2];
  }

  std::vector<CellData> collectCells(
      const MeshBase& mesh,
      const CircleLevelSet& levelSet,
      Real epsilon,
      PhaseMomentMap momentMap)
  {
    std::vector<CellData> cells;
    cells.reserve(mesh.getCellCount());

    for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
    {
      const auto& cellPolytope = *cellIt;
      const auto& vertices = cellPolytope.getVertices();
      if (vertices.size() != 3)
        throw std::runtime_error("LevelSetSDRReconstruction expects triangular cells.");

      CellData data;
      data.index = cellPolytope.getIndex();
      for (size_t i = 0; i < 3; ++i)
      {
        data.vertices[i] = vertices[i];
        data.x[i] = fromPoint(mesh.getVertexCoordinates(vertices[i]));
      }

      Real det = cross(data.x[1] - data.x[0], data.x[2] - data.x[0]);
      if (det < 0)
      {
        std::swap(data.vertices[1], data.vertices[2]);
        std::swap(data.x[1], data.x[2]);
        det = -det;
      }
      if (det <= 0)
        throw std::runtime_error("LevelSetSDRReconstruction encountered a degenerate cell.");

      data.detScale = det;
      data.area = Real(0.5) * det;

      Real moment = 0;
      for (const auto& bary : TriangleBarycentricQuadrature)
      {
        const Vec2 xq = interpolate(data.x, bary);
        moment += applyPhaseMomentMap(levelSet.phi(xq), epsilon, momentMap);
      }
      data.moment = moment / TriangleBarycentricQuadrature.size();
      cells.push_back(std::move(data));
    }
    return cells;
  }

  Real facetLength(const MeshBase& mesh, Index facet)
  {
    const auto face = mesh.getFace(facet);
    const auto& vertices = face->getVertices();
    const Vec2 a = fromPoint(mesh.getVertexCoordinates(vertices[0]));
    const Vec2 b = fromPoint(mesh.getVertexCoordinates(vertices[1]));
    return norm(b - a);
  }

  Math::SpatialMatrix<Real> movedJacobian(
      const CellData& cell,
      const std::vector<int>& activeIndex,
      const Math::Vector<Real>& u)
  {
    auto displacement = [&](Index vertex)
    {
      const int active = activeIndex[vertex];
      if (active < 0)
        return vec2();
      return vec2(u(2 * active), u(2 * active + 1));
    };

    const Vec2 y0 = cell.x[0] + displacement(cell.vertices[0]);
    const Vec2 y1 = cell.x[1] + displacement(cell.vertices[1]);
    const Vec2 y2 = cell.x[2] + displacement(cell.vertices[2]);

    Math::SpatialMatrix<Real> A(2, 2);
    A(0, 0) = y1(0) - y0(0);
    A(1, 0) = y1(1) - y0(1);
    A(0, 1) = y2(0) - y0(0);
    A(1, 1) = y2(1) - y0(1);
    return A;
  }

  struct Parameters
  {
    Real rhoS = 1;
    Real gamma = 1e-1;
    Real beta = 1e-4;
    Real jMin = 1e-8;
    Real delta = 0;
    Real hRef = 0;
    Real jacobianStep = 1e-6;
    Boolean gaussianWeight = true;
  };

  struct ReconstructionAssembly
  {
    Math::Vector<Real> residual;
    Real sdrNormSquared = 0;
    Real sdrEnergy = 0;
    Real shapeEnergy = 0;
    Real detEnergy = 0;
    Real totalEnergy = 0;
    Real minDetRatio = std::numeric_limits<Real>::infinity();
    Real maxShapeQuality = 0;
  };

  class DisplacementProblem
  {
    public:
      DisplacementProblem(
          const std::vector<CellData>& cells,
          const ClassifiedSignedDistance& distance,
          const CircleLevelSet& levelSet,
          const std::vector<int>& activeIndex,
          Parameters parameters)
        : m_cells(cells),
          m_distance(distance),
          m_levelSet(levelSet),
          m_activeIndex(activeIndex),
          m_parameters(parameters)
      {
        m_activeCount = 0;
        for (const int active : activeIndex)
        {
          if (active >= 0)
            m_activeCount = std::max(m_activeCount, static_cast<size_t>(active + 1));
        }
        m_domainMeasure = 0;
        m_bandMeasure = 0;
        for (const CellData& cell : m_cells)
        {
          m_domainMeasure += cell.area;
          for (const auto& bary : TriangleBarycentricQuadrature)
          {
            const Vec2 xq = interpolate(cell.x, bary);
            if (std::abs(m_distance(toPoint(xq))) <= m_parameters.delta)
            {
              m_bandMeasure +=
                cell.area / static_cast<Real>(TriangleBarycentricQuadrature.size());
            }
          }
        }
        if (m_domainMeasure <= 0)
          throw std::runtime_error("LevelSetSDRReconstruction has zero domain measure.");
        if (m_bandMeasure <= 0)
          throw std::runtime_error("LevelSetSDRReconstruction has empty SDR narrow band.");
        if (m_parameters.hRef <= 0)
        {
          m_parameters.hRef =
            std::sqrt(m_domainMeasure / static_cast<Real>(m_cells.size()));
        }
      }

      size_t getActiveCount() const
      {
        return m_activeCount;
      }

      Vec2 displacement(Index vertex, const Math::Vector<Real>& u) const
      {
        const int active = m_activeIndex[vertex];
        if (active < 0)
          return vec2();
        return vec2(u(2 * active), u(2 * active + 1));
      }

      ReconstructionAssembly assemble(const Math::Vector<Real>& u) const
      {
        ReconstructionAssembly out;
        out.residual = Math::Vector<Real>::Zero(2 * m_activeCount);

        for (const CellData& cell : m_cells)
        {
          addSDR(cell, u, out);
          addShapeAndBarrier(cell, u, out);
        }
        out.sdrEnergy = Real(0.5) * out.sdrNormSquared;
        out.totalEnergy =
          m_parameters.rhoS * out.sdrEnergy
          + m_parameters.gamma * out.shapeEnergy
          + m_parameters.beta * out.detEnergy;
        return out;
      }

      const Parameters& getParameters() const
      {
        return m_parameters;
      }

      Real getDomainMeasure() const
      {
        return m_domainMeasure;
      }

      Real getBandMeasure() const
      {
        return m_bandMeasure;
      }

    private:
      void addToResidual(ReconstructionAssembly& out, Index vertex, const Vec2& value) const
      {
        const int active = m_activeIndex[vertex];
        if (active < 0)
          return;
        out.residual(2 * active) += value(0);
        out.residual(2 * active + 1) += value(1);
      }

      void addSDR(
          const CellData& cell,
          const Math::Vector<Real>& u,
          ReconstructionAssembly& out) const
      {
        std::array<Vec2, 3> displacementAtVertex;
        for (size_t i = 0; i < 3; ++i)
          displacementAtVertex[i] = displacement(cell.vertices[i], u);

        for (const auto& bary : TriangleBarycentricQuadrature)
        {
          const Vec2 xq = interpolate(cell.x, bary);
          const Real s = m_distance(toPoint(xq));
          if (std::abs(s) > m_parameters.delta)
            continue;

          const Vec2 uq = interpolate(displacementAtVertex, bary);
          const Vec2 yq = xq + uq;
          const Real residual = m_levelSet.phi(yq) - s;
          const Vec2 gradPhi = m_levelSet.grad(yq);
          const Real weight =
            m_parameters.gaussianWeight
            ? std::exp(-s * s / (2 * m_parameters.delta * m_parameters.delta))
            : Real(1);
          const Real physicalWeight =
            cell.area / static_cast<Real>(TriangleBarycentricQuadrature.size());
          const Real normalizer =
            Real(1) / (m_bandMeasure * m_parameters.hRef * m_parameters.hRef);

          out.sdrNormSquared +=
            weight * residual * residual * physicalWeight * normalizer;
          for (size_t i = 0; i < 3; ++i)
          {
            addToResidual(
                out,
                cell.vertices[i],
                (m_parameters.rhoS * weight * residual * bary[i] * physicalWeight * normalizer) * gradPhi);
          }
        }
      }

      void addShapeAndBarrier(
          const CellData& cell,
          const Math::Vector<Real>& u,
          ReconstructionAssembly& out) const
      {
        const Math::SpatialMatrix<Real> A = movedJacobian(cell, m_activeIndex, u);
        const Real J = A.determinant();
        const Real j = J / cell.detScale;
        if (J <= 0 || j <= m_parameters.jMin)
        {
          throw std::runtime_error(
              "Invalid determinant in LevelSetSDRReconstruction at cell "
              + std::to_string(cell.index) + ": det(A)=" + std::to_string(J)
              + ", det(A)/J_scale=" + std::to_string(j));
        }

        const Math::SpatialMatrix<Real> invT = A.inverse().transpose();
        const Real n = A.squaredNorm();
        const Real qShape = Real(0.5) * n / J;
        out.minDetRatio = std::min(out.minDetRatio, j);
        out.maxShapeQuality = std::max(out.maxShapeQuality, qShape);
        const Real averageWeight = cell.area / m_domainMeasure;
        out.shapeEnergy += averageWeight * (qShape - Real(1));
        out.detEnergy += -averageWeight * std::log(j - m_parameters.jMin);

        Math::SpatialMatrix<Real> PShape =
          A - Real(0.5) * n * invT;
        PShape *= Real(1) / J;

        Math::SpatialMatrix<Real> PDet = invT;
        PDet *= -(j / (j - m_parameters.jMin));

        const Math::SpatialMatrix<Real> P =
          averageWeight * (m_parameters.gamma * PShape + m_parameters.beta * PDet);

        const Vec2 col0 = vec2(P(0, 0), P(1, 0));
        const Vec2 col1 = vec2(P(0, 1), P(1, 1));
        addToResidual(
            out,
            cell.vertices[0],
            -(col0 + col1));
        addToResidual(out, cell.vertices[1], col0);
        addToResidual(out, cell.vertices[2], col1);
      }

      const std::vector<CellData>& m_cells;
      const ClassifiedSignedDistance& m_distance;
      const CircleLevelSet& m_levelSet;
      const std::vector<int>& m_activeIndex;
      Parameters m_parameters;
      size_t m_activeCount = 0;
      Real m_domainMeasure = 0;
      Real m_bandMeasure = 0;
  };

  std::vector<int> buildActiveVertexMap(
      const MeshBase& mesh,
      const ClassifiedSignedDistance& sLF,
      Real activeBand)
  {
    std::vector<int> activeIndex(mesh.getVertexCount(), -1);
    int active = 0;
    for (Index v = 0; v < mesh.getVertexCount(); ++v)
    {
      const Real s = sLF(mesh.getVertexCoordinates(v));
      if (std::abs(s) <= activeBand)
        activeIndex[v] = active++;
    }
    return activeIndex;
  }

  Math::Matrix<Real> finiteDifferenceJacobian(
      const DisplacementProblem& problem,
      const Math::Vector<Real>& u,
      const Math::Vector<Real>& residual)
  {
    Math::Matrix<Real> jacobian(residual.size(), residual.size());
    for (Index i = 0; i < static_cast<Index>(residual.size()); ++i)
    {
      Math::Vector<Real> plus = u;
      Math::Vector<Real> minus = u;
      const Real step =
        problem.getParameters().jacobianStep
        * std::max(problem.getParameters().hRef, std::abs(u(i)));
      plus(i) += step;
      minus(i) -= step;
      jacobian.col(i) =
        (problem.assemble(plus).residual - problem.assemble(minus).residual)
        / (2 * step);
    }
    return jacobian;
  }

  Math::Vector<Real> deterministicDirection(size_t size, Real phase)
  {
    Math::Vector<Real> out(size);
    for (Index i = 0; i < static_cast<Index>(size); ++i)
    {
      const Real k = static_cast<Real>(i + 1);
      out(i) = std::sin((phase + Real(0.37)) * k)
             + Real(0.5) * std::cos((phase + Real(0.11)) * (k + Real(1)));
    }
    const Real nrm = out.norm();
    if (nrm > 0)
      out /= nrm;
    return out;
  }

  Math::Vector<Real> deterministicState(size_t size, Real hRef, Real scale, Real phase)
  {
    return (scale * hRef) * deterministicDirection(size, phase);
  }

  Real relativeVectorError(
      const Math::Vector<Real>& a,
      const Math::Vector<Real>& b)
  {
    return (a - b).norm()
      / std::max({ Real(1e-14), a.norm(), b.norm() });
  }

  struct FDCheckSummary
  {
    Real maxRelativeError = 0;
  };

  struct NewtonConvergenceSummary
  {
    Real finalObservedOrder = 0;
    Boolean quadratic = false;
  };

  FDCheckSummary runFDConsistencyChecks(
      const DisplacementProblem& problem,
      const Math::Vector<Real>& finalState)
  {
    const Real hRef = problem.getParameters().hRef;
    const std::array<Real, 3> epsFactors{{ 1e-3, 1e-4, 1e-5 }};

    struct Case
    {
      std::string name;
      Math::Vector<Real> state;
      Real phase = 0;
    };

    std::vector<Case> cases;
    cases.push_back({ "zero", Math::Vector<Real>::Zero(finalState.size()), Real(0.17) });
    cases.push_back({
        "small-perturbation",
        deterministicState(finalState.size(), hRef, Real(2e-2), Real(0.53)),
        Real(0.71) });
    cases.push_back({ "newton-solution", finalState, Real(1.31) });

    FDCheckSummary summary;
    std::cout << "\nResidual-Jacobian FD checks\n";
    for (const Case& c : cases)
    {
      const ReconstructionAssembly assembly = problem.assemble(c.state);
      const Math::Matrix<Real> jacobian =
        finiteDifferenceJacobian(problem, c.state, assembly.residual);
      const Math::Vector<Real> direction =
        deterministicDirection(finalState.size(), c.phase);
      const Math::Vector<Real> jv = jacobian * direction;

      for (const Real epsFactor : epsFactors)
      {
        const Real eps = epsFactor * hRef;
        const Math::Vector<Real> plus = c.state + eps * direction;
        const Math::Vector<Real> minus = c.state - eps * direction;
        const Math::Vector<Real> fd =
          (problem.assemble(plus).residual - problem.assemble(minus).residual)
          / (2 * eps);
        const Real rel = relativeVectorError(fd, jv);
        summary.maxRelativeError = std::max(summary.maxRelativeError, rel);
        std::cout
          << "  " << c.name
          << ", eps/h=" << epsFactor
          << ": rel ||FD-Jv||=" << rel << '\n';
      }
    }
    std::cout << "  max relative FD error: " << summary.maxRelativeError << '\n';
    if (summary.maxRelativeError > 1e-5)
    {
      throw std::runtime_error(
          "Residual-Jacobian FD consistency check failed: max relative error "
          + std::to_string(summary.maxRelativeError));
    }
    return summary;
  }

  NewtonConvergenceSummary printNewtonConvergenceDiagnostics(
      const std::vector<Real>& residuals)
  {
    NewtonConvergenceSummary summary;
    if (residuals.size() < 3)
      return summary;

    std::cout << "\nNewton local convergence checks\n";
    for (size_t i = 1; i + 1 < residuals.size(); ++i)
    {
      if (residuals[i - 1] <= 0 || residuals[i] <= 0 || residuals[i + 1] <= 0)
        continue;
      const Real observedOrder =
        std::log(residuals[i + 1] / residuals[i])
        / std::log(residuals[i] / residuals[i - 1]);
      const Real quadraticRatio =
        residuals[i + 1] / (residuals[i] * residuals[i]);
      summary.finalObservedOrder = observedOrder;
      std::cout
        << "  it " << i << " -> " << (i + 1)
        << ": order~" << observedOrder
        << ", r_next/r^2=" << quadraticRatio << '\n';
    }
    summary.quadratic = summary.finalObservedOrder >= Real(1.9);
    std::cout
      << "  quadratic convergence verified: "
      << (summary.quadratic ? "yes" : "no")
      << " (final order~" << summary.finalObservedOrder << ")\n";
    if (!summary.quadratic)
      throw std::runtime_error("Newton did not reach the quadratic convergence regime.");
    return summary;
  }

  using DenseLinearSystem =
    Math::LinearSystem<Math::Matrix<Real>, Math::Vector<Real>>;

  class ReconstructionNewtonProblem final
    : public Variational::ProblemBase<DenseLinearSystem>
  {
    public:
      using Parent = Variational::ProblemBase<DenseLinearSystem>;
      using ProblemBodyType = typename Parent::ProblemBodyType;

      ReconstructionNewtonProblem(
          DisplacementProblem& problem,
          Math::Vector<Real>& solution)
        : m_problem(problem),
          m_solution(solution)
      {}

      Parent& operator=(const ProblemBodyType&) override
      {
        return *this;
      }

      void solve(Solver::LinearSolverBase<DenseLinearSystem>& solver) override
      {
        assemble();
        solver.solve(m_system);
      }

      ReconstructionNewtonProblem& assemble() override
      {
        m_lastAssembly = m_problem.get().assemble(m_solution.get());
        m_system.getOperator() =
          finiteDifferenceJacobian(
              m_problem.get(), m_solution.get(), m_lastAssembly.residual);
        m_system.getVector() = -m_lastAssembly.residual;
        m_system.getSolution().resize(m_system.getVector().size());
        return *this;
      }

      DenseLinearSystem& getLinearSystem() override
      {
        return m_system;
      }

      const DenseLinearSystem& getLinearSystem() const override
      {
        return m_system;
      }

      const ReconstructionAssembly& getLastAssembly() const
      {
        return m_lastAssembly;
      }

      ReconstructionNewtonProblem* copy() const noexcept override
      {
        return new ReconstructionNewtonProblem(*this);
      }

    private:
      std::reference_wrapper<DisplacementProblem> m_problem;
      std::reference_wrapper<Math::Vector<Real>> m_solution;
      DenseLinearSystem m_system;
      ReconstructionAssembly m_lastAssembly;
  };

  Real interfacePhiNorm(
      const MeshBase& mesh,
      const ClassifiedSignedDistance& sLF,
      const CircleLevelSet& levelSet,
      const DisplacementProblem& problem,
      const Math::Vector<Real>& u)
  {
    Real weightedSum = 0;
    Real totalWeight = 0;
    for (const Index facet : sLF.getInterfaceFacets())
    {
      const auto face = mesh.getFace(facet);
      const auto& vertices = face->getVertices();
      const Vec2 a = fromPoint(mesh.getVertexCoordinates(vertices[0]));
      const Vec2 b = fromPoint(mesh.getVertexCoordinates(vertices[1]));
      const Vec2 midpoint = Real(0.5) * (a + b);
      const Vec2 umid =
        Real(0.5) * (problem.displacement(vertices[0], u) + problem.displacement(vertices[1], u));
      const Real value = levelSet.phi(midpoint + umid);
      const Real length = norm(b - a);
      weightedSum += value * value * length;
      totalWeight += length;
    }
    return std::sqrt(weightedSum / std::max(totalWeight, Real(1e-30)));
  }

  void writeVertexCSV(
      const std::string& filename,
      const MeshBase& mesh,
      const ClassifiedSignedDistance& sLF,
      const CircleLevelSet& levelSet,
      const DisplacementProblem& problem,
      const std::vector<int>& activeIndex,
      const Math::Vector<Real>& u)
  {
    std::ofstream out(filename);
    out << std::setprecision(16);
    out << "vertex,x,y,active,u_x,u_y,s_lf,phi,phi_moved\n";
    for (Index v = 0; v < mesh.getVertexCount(); ++v)
    {
      const Vec2 x = fromPoint(mesh.getVertexCoordinates(v));
      const Vec2 uv = problem.displacement(v, u);
      out << v << ','
          << x(0) << ',' << x(1) << ','
          << (activeIndex[v] >= 0 ? 1 : 0) << ','
          << uv(0) << ',' << uv(1) << ','
          << sLF(toPoint(x)) << ','
          << levelSet.phi(x) << ','
          << levelSet.phi(x + uv) << '\n';
    }
  }

  void writeCellCSV(
      const std::string& filename,
      const std::vector<CellData>& cells,
      const std::vector<int>& labels,
      const ClassifiedSignedDistance& sLF,
      const CircleLevelSet& levelSet,
      const std::vector<int>& activeIndex,
      const Math::Vector<Real>& u)
  {
    std::ofstream out(filename);
    out << std::setprecision(16);
    out << "cell,label,moment,centroid_x,centroid_y,s_lf,phi,phi_moved,det_ratio,q_shape\n";
    for (const CellData& cell : cells)
    {
      const std::array<Real, 3> bary{{ Real(1) / 3, Real(1) / 3, Real(1) / 3 }};
      const Vec2 centroid = interpolate(cell.x, bary);
      std::array<Vec2, 3> displacementAtVertex;
      for (size_t i = 0; i < 3; ++i)
      {
        const int active = activeIndex[cell.vertices[i]];
        displacementAtVertex[i] =
          active < 0 ? vec2() : vec2(u(2 * active), u(2 * active + 1));
      }
      const Vec2 uCentroid = interpolate(displacementAtVertex, bary);
      const Math::SpatialMatrix<Real> A = movedJacobian(cell, activeIndex, u);
      const Real J = A.determinant();
      const Real qShape = Real(0.5) * A.squaredNorm() / J;
      out << cell.index << ','
          << labels[cell.index] << ','
          << cell.moment << ','
          << centroid(0) << ',' << centroid(1) << ','
          << sLF(toPoint(centroid)) << ','
          << levelSet.phi(centroid) << ','
          << levelSet.phi(centroid + uCentroid) << ','
          << J / cell.detScale << ','
          << qShape << '\n';
    }
  }

  void writeXDMF(
      const std::string& stem,
      const LocalMesh& mesh,
      const LocalMesh& moved,
      const std::vector<CellData>& cells,
      const std::vector<int>& labels,
      const ClassifiedSignedDistance& sLF,
      const CircleLevelSet& levelSet,
      const DisplacementProblem& problem,
      const std::vector<int>& activeIndex,
      const Math::Vector<Real>& u)
  {
    Variational::P0 classifiedCellFES(mesh);
    Variational::GridFunction cellLabel(classifiedCellFES);
    Variational::GridFunction phaseMoment(classifiedCellFES);
    for (const CellData& cell : cells)
    {
      cellLabel[cell.index] = labels[cell.index];
      phaseMoment[cell.index] = cell.moment;
    }

    Variational::P1 classifiedNodeFES(mesh);
    Variational::GridFunction phi(classifiedNodeFES);
    Variational::GridFunction signedDistance(classifiedNodeFES);
    phi = [&](const Point& p)
    {
      return levelSet.phi(fromPoint(p.getCoordinates()));
    };
    signedDistance = [&](const Point& p)
    {
      return sLF(p);
    };

    Variational::P0 fittedCellFES(moved);
    Variational::GridFunction detRatio(fittedCellFES);
    Variational::GridFunction shapeQuality(fittedCellFES);
    for (const CellData& cell : cells)
    {
      const Math::SpatialMatrix<Real> A = movedJacobian(cell, activeIndex, u);
      const Real J = A.determinant();
      detRatio[cell.index] = J / cell.detScale;
      shapeQuality[cell.index] = Real(0.5) * A.squaredNorm() / J;
    }

    Variational::P1 fittedNodeFES(moved);
    Variational::GridFunction phiMoved(fittedNodeFES);
    phiMoved = [&](const Point& p)
    {
      return levelSet.phi(fromPoint(p.getCoordinates()));
    };

    Variational::P1 vectorFES(moved, moved.getSpaceDimension());
    Variational::GridFunction displacement(vectorFES);
    const Index vertexCount = moved.getVertexCount();
    for (Index vertex = 0; vertex < vertexCount; ++vertex)
    {
      const Vec2 uv = problem.displacement(vertex, u);
      displacement[vertex] = uv(0);
      displacement[vertex + vertexCount] = uv(1);
    }

    IO::XDMF xdmf(stem);
    auto classified = xdmf.grid("classified");
    classified.setMesh(mesh);
    classified.add("cell_label", cellLabel, IO::XDMF::Center::Cell);
    classified.add("phase_moment", phaseMoment, IO::XDMF::Center::Cell);
    classified.add("phi", phi, IO::XDMF::Center::Node);
    classified.add("s_lf", signedDistance, IO::XDMF::Center::Node);

    auto fitted = xdmf.grid("fitted");
    fitted.setMesh(moved);
    fitted.add("phi_moved", phiMoved, IO::XDMF::Center::Node);
    fitted.add("displacement", displacement, IO::XDMF::Center::Node);
    fitted.add("det_ratio", detRatio, IO::XDMF::Center::Cell);
    fitted.add("shape_quality", shapeQuality, IO::XDMF::Center::Cell);

    xdmf.write().close();
  }
}

int main(int, char**)
{
  constexpr size_t n = 17;
  const Real h = Real(1) / static_cast<Real>(n - 1);
  const Real epsilon = 1.25 * h;
  const Real lambdaC = 0.008;

  LocalMesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { n, n });
  mesh.scale(h);
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 2);
  mesh.getConnectivity().compute(2, 2);

  const CircleLevelSet levelSet;
  const auto cells =
    collectCells(mesh, levelSet, epsilon, PhaseMomentMap::Tanh);

  std::vector<Real> volumes(cells.size());
  std::vector<Real> moments(cells.size());
  for (const CellData& cell : cells)
  {
    volumes[cell.index] = cell.area;
    moments[cell.index] = cell.moment;
  }

  std::vector<MinSTCut::Edge> graphEdges;
  for (auto faceIt = mesh.getFace(); faceIt; ++faceIt)
  {
    const Index facet = faceIt->getIndex();
    const auto& incident = mesh.getConnectivity().getIncidence({ 1, 2 }, facet);
    if (incident.size() == 2)
    {
      graphEdges.push_back({
          incident[0],
          incident[1],
          lambdaC * facetLength(mesh, facet),
          facet });
    }
  }

  const MinSTCut cut;
  const MinSTCut::Result classified =
    cut.classify(volumes, moments, graphEdges);

  const auto insideComponents =
    mesh.ccl(2,
        [&](const Polytope&, const Polytope&)
        {
          return true;
        },
        [&](const Polytope& cell)
        {
          return classified.labels[cell.getIndex()] == MinSTCut::Inside;
        });

  std::vector<Index> interfaceFacets;
  interfaceFacets.reserve(classified.cutEdges.size());
  for (const MinSTCut::Edge& edge : classified.cutEdges)
  {
    if (edge.index != MinSTCut::InvalidIndex)
      interfaceFacets.push_back(edge.index);
  }

  for (Index cell = 0; cell < classified.labels.size(); ++cell)
  {
    mesh.setAttribute(
        { mesh.getDimension(), cell },
        classified.labels[cell] == MinSTCut::Inside ? 1 : 2);
  }
  for (const Index facet : interfaceFacets)
    mesh.setAttribute({ mesh.getDimension() - 1, facet }, 10);

  mesh.save("LevelSetSDRReconstruction_LF.mesh", IO::FileFormat::MFEM);

  const ClassifiedSignedDistance sLF(mesh, classified.labels, interfaceFacets);

  Parameters parameters;
  parameters.delta = 1.75 * h;
  parameters.hRef = h;
  const Real activeBand = 2.75 * h;
  const auto activeIndex = buildActiveVertexMap(mesh, sLF, activeBand);

  DisplacementProblem problem(
      cells, sLF, levelSet, activeIndex, parameters);

  Math::Vector<Real> u =
    Math::Vector<Real>::Zero(2 * problem.getActiveCount());
  if (u.size() == 0)
    throw std::runtime_error("LevelSetSDRReconstruction found no active displacement dofs.");

  const ReconstructionAssembly initial = problem.assemble(u);
  const Real initialInterfacePhi =
    interfacePhiNorm(mesh, sLF, levelSet, problem, u);
  const Real initialSDR = std::sqrt(initial.sdrNormSquared);

  constexpr size_t maxNewtonIterations = 8;
  constexpr Real absoluteTolerance = 1e-9;
  constexpr Real relativeTolerance = 1e-7;

  ReconstructionNewtonProblem newtonProblem(problem, u);
  Solver::HouseholderQR linearSolver(newtonProblem);
  Solver::NewtonSolver newton(linearSolver);
  std::vector<Real> newtonResiduals;
  newton
    .setMaxIterations(maxNewtonIterations)
    .setDampingFactor(1.0)
    .setStepTolerance(0.0)
    .setAbsoluteTolerance(absoluteTolerance)
    .setRelativeTolerance(relativeTolerance)
    .setMonitor([&](const auto& report)
    {
      const ReconstructionAssembly& assembly = newtonProblem.getLastAssembly();
      std::cout
        << "Newton " << report.iterations
        << ": ||R||=" << report.final_residual
        << ", SDR_RMS=" << std::sqrt(assembly.sdrNormSquared)
        << ", E=" << assembly.totalEnergy
        << ", min j=" << assembly.minDetRatio
        << ", max Q_shape=" << assembly.maxShapeQuality
        << '\n';
      newtonResiduals.push_back(report.final_residual);
    });

  newton.solve(u);
  const auto newtonReport = newton.getReport();

  const ReconstructionAssembly final = problem.assemble(u);
  const Real finalInterfacePhi =
    interfacePhiNorm(mesh, sLF, levelSet, problem, u);
  const Real finalSDR = std::sqrt(final.sdrNormSquared);
  const NewtonConvergenceSummary newtonCheck =
    printNewtonConvergenceDiagnostics(newtonResiduals);
  const FDCheckSummary fdCheck = runFDConsistencyChecks(problem, u);

  LocalMesh moved(mesh);
  for (Index vertex = 0; vertex < moved.getVertexCount(); ++vertex)
  {
    const Vec2 x = fromPoint(mesh.getVertexCoordinates(vertex));
    const Vec2 uv = problem.displacement(vertex, u);
    moved.setVertexCoordinates(vertex, toPoint(x + uv));
  }
  moved.save("LevelSetSDRReconstruction_HF.mesh", IO::FileFormat::MFEM);
  writeXDMF(
      "LevelSetSDRReconstruction",
      mesh,
      moved,
      cells,
      classified.labels,
      sLF,
      levelSet,
      problem,
      activeIndex,
      u);

  writeVertexCSV(
      "LevelSetSDRReconstruction_vertices.csv",
      mesh, sLF, levelSet, problem, activeIndex, u);
  writeCellCSV(
      "LevelSetSDRReconstruction_cells.csv",
      cells, classified.labels, sLF, levelSet, activeIndex, u);

  std::cout << "\nDiagnostics\n";
  std::cout << "  inside cells: " << classified.insideCells.size() << '\n';
  std::cout << "  outside cells: " << classified.outsideCells.size() << '\n';
  std::cout << "  inside components: " << insideComponents.getCount() << '\n';
  std::cout << "  interface facets: " << interfaceFacets.size() << '\n';
  std::cout << "  active displacement dofs: " << u.size() << '\n';
  std::cout << "  h_ref: " << problem.getParameters().hRef << '\n';
  std::cout << "  |Omega_h|: " << problem.getDomainMeasure() << '\n';
  std::cout << "  |Omega_delta|: " << problem.getBandMeasure() << '\n';
  std::cout << "  initial ||phi|| on Gamma_h^LF: " << initialInterfacePhi << '\n';
  std::cout << "  final ||phi(X+u)|| on Gamma_h^LF: " << finalInterfacePhi << '\n';
  std::cout << "  initial normalized SDR residual RMS: " << initialSDR << '\n';
  std::cout << "  final normalized SDR residual RMS: " << finalSDR << '\n';
  std::cout << "  initial weighted energies:"
            << " SDR=" << parameters.rhoS * initial.sdrEnergy
            << ", shape=" << parameters.gamma * initial.shapeEnergy
            << ", det=" << parameters.beta * initial.detEnergy
            << ", total=" << initial.totalEnergy << '\n';
  std::cout << "  final weighted energies:"
            << " SDR=" << parameters.rhoS * final.sdrEnergy
            << ", shape=" << parameters.gamma * final.shapeEnergy
            << ", det=" << parameters.beta * final.detEnergy
            << ", total=" << final.totalEnergy << '\n';
  std::cout << "  min j_K^u: " << final.minDetRatio << '\n';
  std::cout << "  max Q_shape: " << final.maxShapeQuality << '\n';
  std::cout << "  final Newton observed order: "
            << newtonCheck.finalObservedOrder << '\n';
  std::cout << "  max residual-Jacobian FD relative error: "
            << fdCheck.maxRelativeError << '\n';
  std::cout << "  Newton iteration count: " << newtonReport.iterations << '\n';
  std::cout << "  Newton converged: " << (newtonReport.converged ? "yes" : "no") << '\n';
  std::cout << "  wrote LevelSetSDRReconstruction_LF.mesh\n";
  std::cout << "  wrote LevelSetSDRReconstruction_HF.mesh\n";
  std::cout << "  wrote LevelSetSDRReconstruction.xdmf\n";
  std::cout << "  wrote LevelSetSDRReconstruction_vertices.csv\n";
  std::cout << "  wrote LevelSetSDRReconstruction_cells.csv\n";

  return 0;
}
