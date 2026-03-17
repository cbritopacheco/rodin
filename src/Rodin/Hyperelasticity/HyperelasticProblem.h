/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file HyperelasticProblem.h
 * @brief Nonlinear hyperelastic problem with Newton-Raphson solver.
 *
 * Assembles the internal-force residual and tangent stiffness matrix for a
 * displacement-based finite-strain formulation and solves via Newton
 * iterations with a user-chosen linear solver back-end.
 *
 * ## Weak form (total Lagrangian)
 * Find @f$ \mathbf{u} \in V @f$ such that for all @f$ \mathbf{v} \in V @f$:
 * @f[
 *   \int_{\Omega_0} \mathbf{P}(\mathbf{F}) : \nabla \mathbf{v}\, dX
 *   = \int_{\Omega_0} \mathbf{f}_0 \cdot \mathbf{v}\, dX
 *   + \int_{\Gamma_N} \mathbf{t}_0 \cdot \mathbf{v}\, dS
 * @f]
 * where @f$ \mathbf{F} = \mathbf{I} + \nabla \mathbf{u} @f$.
 */
#ifndef RODIN_HYPERELASTICITY_HYPERELASTICPROBLEM_H
#define RODIN_HYPERELASTICITY_HYPERELASTICPROBLEM_H

#include <vector>
#include <functional>
#include <cmath>
#include <cassert>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

#include "Rodin/Types.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"

#include "Rodin/Alert/Info.h"
#include "Rodin/Alert/Raise.h"
#include "Rodin/Alert/Warning.h"
#include "Rodin/Alert/Success.h"

#include "Rodin/Geometry/Mesh.h"

#include "Rodin/Variational/P1/P1.h"
#include "Rodin/Variational/FiniteElementSpace.h"

#include "ForwardDecls.h"
#include "ConstitutiveLaw.h"

namespace Rodin::Hyperelasticity
{
  /**
   * @brief Nonlinear hyperelastic problem (displacement-based, total
   *        Lagrangian).
   *
   * Assembles element-level residual and tangent contributions, applies
   * Dirichlet boundary conditions, and solves the resulting nonlinear
   * system with Newton-Raphson iterations using Eigen's SparseLU.
   *
   * @tparam ConstitutiveLaw A type derived from
   *         ConstitutiveLawBase\<ConstitutiveLaw\>.
   *
   * ## Usage
   * @code{.cpp}
   * Mesh mesh;
   * mesh = mesh.UniformGrid(Polytope::Type::Triangle, {8, 8});
   * mesh.getConnectivity().compute(1, 2);
   *
   * size_t d = mesh.getSpaceDimension();
   * P1 fes(mesh, d);
   *
   * NeoHookean material(mu, lambda);
   * HyperelasticProblem problem(fes, material);
   * problem.addDirichletBC(GammaD, [](auto) { return Eigen::Vector2d::Zero(); });
   * problem.addNeumannBC(GammaN, [](auto) { return Eigen::Vector2d(0, -1); });
   * problem.solve();
   * @endcode
   */
  template <class ConstitutiveLaw>
  class HyperelasticProblem final
  {
    public:
      /// Mesh type
      using MeshType = Geometry::Mesh<Context::Local>;

      /// Finite-element space type (vector P1)
      using FESType = Variational::P1<Math::Vector<Real>, MeshType>;

      /// Boundary-value callback: maps a spatial coordinate to a vector value
      using BCFunction = std::function<Eigen::VectorXd(const Eigen::VectorXd&)>;

      /**
       * @brief Constructs a hyperelastic problem.
       * @param fes  Vector-valued P1 finite-element space
       * @param law  Hyperelastic constitutive law
       */
      HyperelasticProblem(const FESType& fes, const ConstitutiveLaw& law)
        : m_fes(fes),
          m_law(law),
          m_mesh(fes.getMesh()),
          m_maxIt(50),
          m_atol(1e-10),
          m_rtol(1e-8)
      {}

      /** @brief Sets maximum Newton iterations. */
      HyperelasticProblem& setMaxIterations(size_t n)
      { m_maxIt = n; return *this; }

      /** @brief Sets absolute convergence tolerance on the residual norm. */
      HyperelasticProblem& setAbsoluteTolerance(Real t)
      { m_atol = t; return *this; }

      /** @brief Sets relative convergence tolerance on the residual norm. */
      HyperelasticProblem& setRelativeTolerance(Real t)
      { m_rtol = t; return *this; }

      /**
       * @brief Adds a Dirichlet boundary condition.
       * @param attr Boundary attribute
       * @param fn   Callback returning the prescribed displacement at a point
       */
      HyperelasticProblem& addDirichletBC(Geometry::Attribute attr, BCFunction fn)
      {
        m_dirichlet.push_back({attr, std::move(fn)});
        return *this;
      }

      /**
       * @brief Adds a Neumann boundary condition (reference traction).
       * @param attr Boundary attribute
       * @param fn   Callback returning the traction at a point
       */
      HyperelasticProblem& addNeumannBC(Geometry::Attribute attr, BCFunction fn)
      {
        m_neumann.push_back({attr, std::move(fn)});
        return *this;
      }

      /**
       * @brief Sets a body force in the reference configuration.
       * @param fn Callback returning the body force at a point
       */
      HyperelasticProblem& setBodyForce(BCFunction fn)
      {
        m_bodyForce = std::move(fn);
        return *this;
      }

      /**
       * @brief Solves the nonlinear problem with Newton-Raphson iterations.
       *
       * Uses Eigen SparseLU as the linear solver back-end for each Newton
       * step.
       */
      void solve()
      {
        const size_t ndof = m_fes.get().getSize();
        m_solution.setZero(ndof);

        // Apply initial Dirichlet values
        applyDirichletValues(m_solution);

        // Identify constrained DOFs
        std::vector<bool> isConstrained(ndof, false);
        markConstrainedDOFs(isConstrained);

        Eigen::VectorXd R(ndof);
        Eigen::SparseMatrix<Real> K(ndof, ndof);

        Real r0 = 0.0;
        for (size_t it = 0; it < m_maxIt; ++it)
        {
          // Assemble residual and tangent
          assembleSystem(m_solution, R, K, isConstrained);

          const Real rNorm = R.norm();
          if (it == 0) r0 = rNorm;

          Alert::Info() << "  Newton it " << it
                        << "  |R| = " << rNorm << Alert::Raise;

          if (rNorm <= m_atol || (r0 > 0.0 && rNorm <= m_rtol * r0))
          {
            Alert::Success() << "Newton converged after " << it
                             << " iterations (|R| = " << rNorm << ")."
                             << Alert::Raise;
            return;
          }

          // Solve K * du = -R
          Eigen::SparseLU<Eigen::SparseMatrix<Real>> solver;
          solver.compute(K);
          if (solver.info() != Eigen::Success)
          {
            Alert::Warning() << "SparseLU factorisation failed at Newton it "
                             << it << "." << Alert::Raise;
            return;
          }
          Eigen::VectorXd du = solver.solve(-R);

          m_solution += du;
        }

        Alert::Warning() << "Newton did not converge within "
                         << m_maxIt << " iterations." << Alert::Raise;
      }

      /** @brief Returns the computed displacement DOF vector. */
      const Eigen::VectorXd& getSolution() const { return m_solution; }

      /** @brief Returns the computed displacement DOF vector (mutable). */
      Eigen::VectorXd& getSolution() { return m_solution; }

    private:
      // ---------------------------------------------------------------
      //  Assembly helpers
      // ---------------------------------------------------------------

      /**
       * @brief Assembles the global residual R and tangent K.
       *
       * Constrained DOFs get their row/column zeroed with a unit diagonal
       * and a zero residual entry.
       */
      void assembleSystem(
          const Eigen::VectorXd& u,
          Eigen::VectorXd& R,
          Eigen::SparseMatrix<Real>& K,
          const std::vector<bool>& isConstrained) const
      {
        const auto& mesh = m_mesh.get();
        const auto& fes  = m_fes.get();
        const size_t d    = mesh.getSpaceDimension();
        const size_t nv   = mesh.getVertexCount();
        const size_t ndof = fes.getSize();
        const size_t dim  = mesh.getDimension();

        R.setZero(ndof);

        // Triplets for sparse tangent
        using Triplet = Eigen::Triplet<Real>;
        std::vector<Triplet> triplets;
        triplets.reserve(mesh.getCellCount() * (d * (dim + 1)) * (d * (dim + 1)));

        // Loop over cells
        for (auto cell = mesh.getCell(); !cell.end(); ++cell)
        {
          const Index cellIdx = cell->getIndex();
          const auto& verts = mesh.getConnectivity().getPolytope(dim, cellIdx);
          const size_t npe = static_cast<size_t>(verts.size()); // nodes per element

          // Gather vertex coordinates (reference config) and displacements
          Eigen::MatrixXd X(d, npe);   // reference coordinates
          Eigen::MatrixXd U(d, npe);   // displacement
          for (size_t a = 0; a < npe; ++a)
          {
            const Index vi = verts(static_cast<Eigen::Index>(a));
            auto coords = mesh.getVertexCoordinates(vi);
            for (size_t i = 0; i < d; ++i)
            {
              X(i, a) = coords(i);
              U(i, a) = u(vi + i * nv);   // block DOF ordering
            }
          }

          // Compute element volume and shape-function gradients
          // For P1 on a simplex: gradN_ref are constant per element
          Eigen::MatrixXd dN_dxi;  // (dim x npe) reference gradients
          computeRefGradients(dim, npe, dN_dxi);

          // Jacobian of reference-to-physical mapping: dx/dxi
          Eigen::MatrixXd Jmap = X * dN_dxi.transpose(); // (d x dim)
          const Real detJ = simplexDetJ(Jmap, d, dim);
          assert(detJ > 0.0);
          const Real vol = std::abs(detJ) / factorial(dim);

          // Physical gradients: dN/dX = dN/dxi * inv(dxi/dX)
          // For simplex in same-dim space: dX/dxi = Jmap, so dN/dX = inv(Jmap)^T * dN/dxi^T ... 
          // Actually: grad_X N_a = Jmap^{-T} * grad_xi N_a
          Eigen::MatrixXd JmapInvT;
          if (d == dim)
            JmapInvT = Jmap.inverse().transpose();
          else
            JmapInvT = (Jmap.transpose() * Jmap).inverse() * Jmap.transpose(); // pseudo-inverse transpose

          Eigen::MatrixXd dN_dX = JmapInvT * dN_dxi; // (d x npe)

          // Deformation gradient: F = I + grad_X u = I + U * dN_dX^T
          // grad_X u_{ij} = sum_a U(i,a) * dN_dX(j,a)
          Eigen::MatrixXd gradU = U * dN_dX.transpose(); // (d x d)
          Eigen::MatrixXd F = Eigen::MatrixXd::Identity(d, d) + gradU;

          // Evaluate constitutive law
          Eigen::MatrixXd P(d, d);
          m_law.get().stress(F, P);

          Eigen::MatrixXd C(d * d, d * d);
          m_law.get().tangent(F, C);

          // Build local DOF index map: localDOF(a, i) → global DOF
          // Block ordering: vertex v, component c → v + c * nv
          std::vector<Index> dofMap(npe * d);
          for (size_t a = 0; a < npe; ++a)
            for (size_t i = 0; i < d; ++i)
              dofMap[a * d + i] = verts(static_cast<Eigen::Index>(a)) + i * nv;

          // Assemble residual:  R_I += int P_{iJ} dN_a/dX_J dV
          for (size_t a = 0; a < npe; ++a)
            for (size_t i = 0; i < d; ++i)
            {
              Real val = 0.0;
              for (size_t J = 0; J < d; ++J)
                val += P(i, J) * dN_dX(J, a);
              R(dofMap[a * d + i]) += val * vol;
            }

          // Assemble tangent:
          //   K_{(a,i),(b,k)} += int A_{iJkL} dN_a/dX_J dN_b/dX_L dV
          for (size_t a = 0; a < npe; ++a)
            for (size_t i = 0; i < d; ++i)
              for (size_t b = 0; b < npe; ++b)
                for (size_t k = 0; k < d; ++k)
                {
                  Real val = 0.0;
                  for (size_t J = 0; J < d; ++J)
                    for (size_t L = 0; L < d; ++L)
                      val += C(i * d + J, k * d + L) * dN_dX(J, a) * dN_dX(L, b);
                  triplets.emplace_back(
                    static_cast<int>(dofMap[a * d + i]),
                    static_cast<int>(dofMap[b * d + k]),
                    val * vol);
                }
        }

        // Add Neumann traction to residual (subtract, since R = F_int - F_ext)
        assembleNeumann(u, R);

        // Subtract body-force contribution from residual
        assembleBodyForce(R);

        // Build sparse tangent
        K.resize(ndof, ndof);
        K.setZero();
        K.setFromTriplets(triplets.begin(), triplets.end());

        // Apply Dirichlet BC: zero constrained rows/cols, set diagonal to 1
        for (size_t i = 0; i < ndof; ++i)
        {
          if (!isConstrained[i]) continue;
          R(i) = 0.0;
          for (Eigen::SparseMatrix<Real>::InnerIterator it(K, i); it; ++it)
            it.valueRef() = 0.0;
          K.coeffRef(i, i) = 1.0;
        }
        K.makeCompressed();

        // Also zero columns for constrained DOFs
        Eigen::SparseMatrix<Real> Kt = K.transpose();
        for (size_t i = 0; i < ndof; ++i)
        {
          if (!isConstrained[i]) continue;
          for (Eigen::SparseMatrix<Real>::InnerIterator it(Kt, i); it; ++it)
          {
            if (static_cast<size_t>(it.row()) != i)
              it.valueRef() = 0.0;
          }
        }
        K = Kt.transpose();
        K.makeCompressed();
      }

      /**
       * @brief Assembles Neumann traction into the residual.
       *
       * Subtracts the external traction contribution from R:
       * @f$ R_I -= \int_{\Gamma_N} t_i N_a dS @f$
       */
      void assembleNeumann(
          const Eigen::VectorXd& /*u*/,
          Eigen::VectorXd& R) const
      {
        if (m_neumann.empty()) return;

        const auto& mesh = m_mesh.get();
        const size_t d   = mesh.getSpaceDimension();
        const size_t nv  = mesh.getVertexCount();
        const size_t faceDim = mesh.getDimension() - 1;

        for (auto face = mesh.getBoundary(); !face.end(); ++face)
        {
          const auto attr = face->getAttribute();
          if (!attr.has_value()) continue;

          for (const auto& [bcAttr, fn] : m_neumann)
          {
            if (attr.value() != bcAttr) continue;

            const Index faceIdx = face->getIndex();
            const auto& verts = mesh.getConnectivity().getPolytope(faceDim, faceIdx);
            const size_t nfv = static_cast<size_t>(verts.size());

            // Compute face centroid for traction evaluation
            Eigen::VectorXd centroid = Eigen::VectorXd::Zero(d);
            for (size_t a = 0; a < nfv; ++a)
            {
              auto coords = mesh.getVertexCoordinates(verts(static_cast<Eigen::Index>(a)));
              for (size_t i = 0; i < d; ++i)
                centroid(i) += coords(i);
            }
            centroid /= static_cast<Real>(nfv);

            // Compute face measure (length in 2D, area in 3D)
            Real measure = computeFaceMeasure(mesh, verts, nfv, d);

            // Evaluate traction
            Eigen::VectorXd t = fn(centroid);

            // Distribute equally to face nodes (P1 lumping)
            const Real w = measure / static_cast<Real>(nfv);
            for (size_t a = 0; a < nfv; ++a)
            {
              const Index vi = verts(static_cast<Eigen::Index>(a));
              for (size_t i = 0; i < d; ++i)
                R(vi + i * nv) -= t(i) * w;
            }
          }
        }
      }

      /**
       * @brief Subtracts body-force contribution from the residual.
       */
      void assembleBodyForce(Eigen::VectorXd& R) const
      {
        if (!m_bodyForce) return;

        const auto& mesh = m_mesh.get();
        const size_t d   = mesh.getSpaceDimension();
        const size_t nv  = mesh.getVertexCount();
        const size_t dim = mesh.getDimension();

        for (auto cell = mesh.getCell(); !cell.end(); ++cell)
        {
          const Index cellIdx = cell->getIndex();
          const auto& verts = mesh.getConnectivity().getPolytope(dim, cellIdx);
          const size_t npe = static_cast<size_t>(verts.size());

          Eigen::MatrixXd X(d, npe);
          Eigen::VectorXd centroid = Eigen::VectorXd::Zero(d);
          for (size_t a = 0; a < npe; ++a)
          {
            auto coords = mesh.getVertexCoordinates(verts(static_cast<Eigen::Index>(a)));
            for (size_t i = 0; i < d; ++i)
            {
              X(i, a) = coords(i);
              centroid(i) += coords(i);
            }
          }
          centroid /= static_cast<Real>(npe);

          Eigen::MatrixXd dN_dxi;
          computeRefGradients(dim, npe, dN_dxi);
          Eigen::MatrixXd Jmap = X * dN_dxi.transpose();
          const Real detJ = simplexDetJ(Jmap, d, dim);
          const Real vol = std::abs(detJ) / factorial(dim);

          Eigen::VectorXd f = m_bodyForce(centroid);
          const Real w = vol / static_cast<Real>(npe);
          for (size_t a = 0; a < npe; ++a)
          {
            const Index vi = verts(static_cast<Eigen::Index>(a));
            for (size_t i = 0; i < d; ++i)
              R(vi + i * nv) -= f(i) * w;
          }
        }
      }

      // ---------------------------------------------------------------
      //  Dirichlet helpers
      // ---------------------------------------------------------------

      void markConstrainedDOFs(std::vector<bool>& isConstrained) const
      {
        const auto& mesh = m_mesh.get();
        const size_t d   = mesh.getSpaceDimension();
        const size_t nv  = mesh.getVertexCount();
        const size_t faceDim = mesh.getDimension() - 1;

        for (auto face = mesh.getBoundary(); !face.end(); ++face)
        {
          const auto attr = face->getAttribute();
          if (!attr.has_value()) continue;

          for (const auto& [bcAttr, fn] : m_dirichlet)
          {
            if (attr.value() != bcAttr) continue;

            const Index faceIdx = face->getIndex();
            const auto& verts = mesh.getConnectivity().getPolytope(faceDim, faceIdx);
            for (Eigen::Index a = 0; a < verts.size(); ++a)
            {
              const Index vi = verts(a);
              for (size_t c = 0; c < d; ++c)
                isConstrained[vi + c * nv] = true;
            }
          }
        }
      }

      void applyDirichletValues(Eigen::VectorXd& u) const
      {
        const auto& mesh = m_mesh.get();
        const size_t d   = mesh.getSpaceDimension();
        const size_t nv  = mesh.getVertexCount();
        const size_t faceDim = mesh.getDimension() - 1;

        for (auto face = mesh.getBoundary(); !face.end(); ++face)
        {
          const auto attr = face->getAttribute();
          if (!attr.has_value()) continue;

          for (const auto& [bcAttr, fn] : m_dirichlet)
          {
            if (attr.value() != bcAttr) continue;

            const Index faceIdx = face->getIndex();
            const auto& verts = mesh.getConnectivity().getPolytope(faceDim, faceIdx);
            for (Eigen::Index a = 0; a < verts.size(); ++a)
            {
              const Index vi = verts(a);
              auto coords = mesh.getVertexCoordinates(vi);
              Eigen::VectorXd pt(d);
              for (size_t c = 0; c < d; ++c)
                pt(c) = coords(c);

              Eigen::VectorXd val = fn(pt);
              for (size_t c = 0; c < d; ++c)
                u(vi + c * nv) = val(c);
            }
          }
        }
      }

      // ---------------------------------------------------------------
      //  Geometry helpers
      // ---------------------------------------------------------------

      /**
       * @brief Reference-element shape-function gradients for P1 simplices.
       *
       * For a simplex with npe vertices in dim-dimensional reference space
       * the gradients are constant:
       *   N_0 = 1 - xi_1 - ... - xi_dim
       *   N_a = xi_a   for a = 1..dim
       */
      static void computeRefGradients(size_t dim, size_t npe,
                                      Eigen::MatrixXd& dN_dxi)
      {
        assert(npe == dim + 1); // P1 simplex
        (void)npe;
        dN_dxi.setZero(dim, dim + 1);
        for (size_t i = 0; i < dim; ++i)
        {
          dN_dxi(i, 0) = -1.0;           // dN_0/dxi_i
          dN_dxi(i, i + 1) = 1.0;        // dN_{i+1}/dxi_i
        }
      }

      /**
       * @brief Determinant of the Jacobian for a simplex mapping.
       *
       * When d == dim, this is the ordinary determinant.
       * When d > dim (e.g. surface mesh), uses the Gram determinant.
       */
      static Real simplexDetJ(const Eigen::MatrixXd& J,
                               size_t d, size_t dim)
      {
        if (d == dim)
          return J.determinant();
        // Gram determinant for embedded simplices
        return std::sqrt((J.transpose() * J).determinant());
      }

      /** @brief Factorial helper. */
      static Real factorial(size_t n)
      {
        Real f = 1.0;
        for (size_t k = 2; k <= n; ++k)
          f *= static_cast<Real>(k);
        return f;
      }

      /**
       * @brief Computes the measure (length / area) of a boundary face.
       */
      static Real computeFaceMeasure(
          const MeshType& mesh,
          const Eigen::Ref<const Eigen::VectorXi>& verts,
          size_t nfv, size_t d)
      {
        if (nfv == 2)
        {
          // Edge in 2D
          auto p0 = mesh.getVertexCoordinates(verts(0));
          auto p1 = mesh.getVertexCoordinates(verts(1));
          Real len = 0.0;
          for (size_t i = 0; i < d; ++i)
          {
            Real dx = p1(i) - p0(i);
            len += dx * dx;
          }
          return std::sqrt(len);
        }
        if (nfv == 3)
        {
          // Triangle face in 3D
          auto p0 = mesh.getVertexCoordinates(verts(0));
          auto p1 = mesh.getVertexCoordinates(verts(1));
          auto p2 = mesh.getVertexCoordinates(verts(2));
          Eigen::Vector3d e1, e2;
          e1.setZero(); e2.setZero();
          for (size_t i = 0; i < d && i < 3; ++i)
          {
            e1(i) = p1(i) - p0(i);
            e2(i) = p2(i) - p0(i);
          }
          return 0.5 * e1.cross(e2).norm();
        }
        return 0.0;
      }

      // ---------------------------------------------------------------
      //  Data members
      // ---------------------------------------------------------------

      std::reference_wrapper<const FESType>          m_fes;
      std::reference_wrapper<const ConstitutiveLaw>  m_law;
      std::reference_wrapper<const MeshType>         m_mesh;

      size_t m_maxIt;
      Real   m_atol;
      Real   m_rtol;

      struct BC { Geometry::Attribute attr; BCFunction fn; };
      std::vector<BC> m_dirichlet;
      std::vector<BC> m_neumann;
      BCFunction m_bodyForce;

      Eigen::VectorXd m_solution;
  };
}

#endif
