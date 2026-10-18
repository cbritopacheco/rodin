/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_WORKLOADS_MASSWORKLOAD_H
#define RODIN_WORKLOADS_MASSWORKLOAD_H

/**
 * @file
 * @brief Backend-independent mass/reaction workloads and constant-field oracles.
 */
#include <benchmark/benchmark.h>
#include "Rodin/Variational/H1.h"
#include "Rodin/Variational/P0g.h"
#include "convergence/Convergence.h"
#include "convergence/CurvedGeometry.h"
#include "QuadratureReference.h"

namespace Rodin::Tests::Benchmarks
{
  enum class MassFamily
  {
    H1,
    P1,
    P0,
    P0g
  };

  /**
   * @brief Defines mass and reaction operators and projection right-hand sides.
   *
   * The forms are @f$a_\rho(u,v)=\int_\Omega\rho u\cdot\overline v\,dx@f$
   * and @f$L(v)=a_\rho(c,v)@f$. Mass uses @f$\rho=1@f$; reaction uses
   * @f$\rho=s\alpha(1+\sum_{j=0}^{d-1}x_j)@f$, with @f$\alpha=1@f$
   * over the reals and @f$\alpha=1+i/2@f$ over the complex numbers.
   * For the quadratic image of the unit box, the independent moments are
   * @f[
   * V_d(a)=\int_\Omega1\,dx=
   * \begin{cases}1+a,&d=1,\\1,&d\geq2,\end{cases}
   * \qquad
   * W_d(a)=\int_\Omega(1+\sum_jx_j)\,dx=
   * \begin{cases}(1+a)+(1+a)^2/2,&d=1,\\
   * 1+d/2+a/3,&d\geq2.\end{cases}
   * @f]
   * Hence @f$Mc=b@f$, while @f$c^*Mc=|c|^2V_d(a)@f$ for mass and
   * @f$c^*Mc=s\alpha|c|^2W_d(a)@f$ for reaction. Affine cases have @f$a=0@f$.
   *
   * ## Architecture
   *
   * Geometry is installed before distribution. Spaces use public constructors
   * and constants use DOF functionals, not coefficient-layout assumptions.
   * Rank-local pointwise geometry checks supplement the original-entry-loop
   * matrix oracle, which shares mapped quadrature with production assembly.
   * A referenced coefficient scale permits untimed state updates. Three mesh
   * sizes separate operator from load timing. Fixed order eight defines the
   * implemented quadrature operator; exact full matrix entries on all curved
   * families are not inferred from the constant-field moments.
   *
   * @tparam K Polynomial degree for the H1 family.
   * @tparam Scalar Real or complex coefficient field.
   * @tparam Vector Whether the field has one component per spatial dimension.
   * @tparam Family Finite-element space family.
   */
  template <size_t K, class Scalar, bool Vector, MassFamily Family>
  class MassWorkload
  {
    public:
      /// Dimensionless amplitude of the exact quadratic map.
      static constexpr Real GeometryAmplitude = Real(1) / 10;

      /**
       * @brief Constructs the physical benchmark domain before distribution.
       *
       * Curved cases use @f$\Phi(\xi)=\xi+a\xi_0^2e_{d-1}@f$ with
       * @f$a=0.1@f$. Its measure is one for @f$d\geq2@f$ and
       * @f$1+a@f$ in one dimension. No topology or ownership is changed.
       *
       * @param geometry Reference-cell family.
       * @param pointsPerAxis Grid points per coordinate axis.
       * @param curved Whether to install exact quadratic geometry.
       * @returns The local mesh, ready for use or MPI partitioning.
       */
      static auto makeMesh(
        Geometry::Polytope::Type geometry, size_t pointsPerAxis, bool curved)
      {
        auto mesh = Convergence::UniformGrid(geometry).makeMesh(pointsPerAxis);
        if (curved)
        {
          Convergence::CurvedGeometry map(mesh,
            Convergence::CurvedGeometry<decltype(mesh)>::Map::Quadratic,
            GeometryAmplitude);
          map.template install<2>();
        }
        return mesh;
      }

      /**
       * @brief Checks mapped quadrature against independent pointwise geometry.
       *
       * Every held cell, including MPI halos, is checked without communication.
       * Physical coordinates and Jacobian entries must agree exactly with
       * direct map evaluation at the same reference point. This is a numerical
       * equivalence check, not a geometric identity or ownership decision.
       * Curved cells must also retain factor degree two after distribution.
       *
       * @tparam Mesh Local or distributed mesh type.
       * @param mesh Mesh carrying the geometry being measured.
       * @param geometry Reference-cell family, including on empty ranks.
       * @param curved Whether exact quadratic geometry was requested.
       * @returns Whether every held cell passes the independent checks.
       */
      template <class Mesh>
      static bool checkGeometry(
        const Mesh& mesh, Geometry::Polytope::Type geometry, bool curved)
      {
        const auto& formula = QF::PolytopeQuadratureFormula::get(8, geometry);
        for (auto cell = mesh.getCell(); cell; ++cell)
        {
          if (curved && cell->getTransformation().getFactorOrder() != 2)
            return false;
          const Geometry::PolytopeQuadrature quadrature(*cell, formula);
          for (size_t qp = 0; qp < formula.getSize(); ++qp)
          {
            const Geometry::Point independent(*cell, formula.getPoint(qp));
            const auto& mapped = quadrature.getPoint(qp);
            const auto& x = mapped.getPhysicalCoordinates();
            const auto& reference = independent.getPhysicalCoordinates();
            const auto& jacobian = mapped.getJacobian();
            const auto& referenceJacobian = independent.getJacobian();
            for (size_t j = 0; j < mesh.getSpaceDimension(); ++j)
            {
              if (x(j) != reference(j))
                return false;
              for (size_t i = 0; i < cell->getDimension(); ++i)
              {
                if (jacobian(j, i) != referenceJacobian(j, i))
                  return false;
              }
            }
          }
        }
        return true;
      }

      /**
       * @brief Constructs the requested scalar or vector space.
       * @tparam Mesh Local or distributed mesh type.
       * @param mesh Physical mesh defining the space.
       * @param dim Spatial dimension and vector component count.
       * @returns The requested finite-element space.
       */
      template <class Mesh>
      static auto makeSpace(const Mesh& mesh, size_t dim)
      {
        using Range = std::conditional_t<Vector, Math::SpatialVector<Scalar>, Scalar>;
        using FES =
          std::conditional_t<Family == MassFamily::H1, Variational::H1<K, Range, Mesh>,
            std::conditional_t<Family == MassFamily::P1, Variational::P1<Range, Mesh>,
              std::conditional_t<Family == MassFamily::P0, Variational::P0<Range, Mesh>,
                Variational::P0g<Range, Mesh>>>>;
        if constexpr (Family == MassFamily::H1)
        {
          if constexpr (Vector)
            return FES(std::integral_constant<size_t, K>{}, mesh, dim);
          else
            return FES(std::integral_constant<size_t, K>{}, mesh);
        }
        else if constexpr (Vector)
          return FES(mesh, dim);
        else
          return FES(mesh);
      }

      /**
       * @brief Defines the forms and their independent constant-field action.
       *
       * The coefficient borrows @p scale for untimed replacement checks.
       * A reference configuration changes the matrix-entry algorithm only;
       * independent pointwise geometry is checked separately.
       *
       * @tparam Reference Whether to use the original matrix-entry loop.
       * @tparam U Trial-function type.
       * @tparam V Test-function type.
       * @tparam A Bilinear-form type.
       * @tparam B Linear-form type.
       * @param u Trial function whose solution receives the constant field.
       * @param v Test function.
       * @param matrix Operator being configured.
       * @param load Projection load being configured.
       * @param dim Domain dimension.
       * @param reaction Whether to use the variable reaction coefficient.
       * @param scale Borrowed reaction-coefficient multiplier.
       * @param curved Whether the domain is the quadratic image of the unit box.
       * @returns The analytic action at unit scale.
       */
      template <bool Reference = false, class U, class V, class A, class B>
      static Scalar configure(U& u, V& v, A& matrix, B& load, size_t dim, bool reaction,
        const Real& scale, bool curved)
      {
        using namespace Variational;
        const Scalar phase = []() {
          if constexpr (std::is_same_v<Scalar, Complex>)
            return Scalar(1, 2);
          else
            return Scalar(1);
        }();
        Scalar alpha = 1;
        if constexpr (std::is_same_v<Scalar, Complex>)
          if (reaction)
            alpha = Scalar(1, 0.5);
        auto field = [&]() {
          if constexpr (Vector)
            return VectorFunction(dim, [dim, phase](const Geometry::Point&) {
              Math::SpatialVector<Scalar> value(static_cast<std::uint8_t>(dim));
              for (size_t j = 0; j < dim; ++j)
                value(j) = Scalar(j + 1) * phase;
              return value;
            });
          else if constexpr (std::is_same_v<Scalar, Complex>)
            return ComplexFunction(phase);
          else
            return RealFunction(phase);
        }();
        auto coefficient = [&]() {
          const auto evaluate = [dim, reaction, alpha, &scale](const Geometry::Point& p) {
            Real rho = 1;
            if (reaction)
              for (size_t j = 0; j < dim; ++j)
                rho += p(j);
            return scale * alpha * rho;
          };
          if constexpr (std::is_same_v<Scalar, Complex>)
            return ComplexFunction(evaluate);
          else
            return RealFunction(evaluate);
        }();
        const auto assign = [&matrix](auto integral) {
          integral.setOrder(8);
          if constexpr (Reference)
            matrix = ReferenceIntegral(integral);
          else
            matrix = integral;
        };
        if (reaction)
          assign(Integral(u * coefficient, v));
        else
          assign(Integral(u, v));
        auto rhs = Integral(field * coefficient, v);
        rhs.setOrder(8);
        load = rhs;
        u.getSolution() = field;
        const Real square = Vector ? Real(dim * (dim + 1) * (2 * dim + 1)) / 6 : 1;
        const Real phaseSquare = std::abs(phase) * std::abs(phase);
        const Real amplitude = curved ? GeometryAmplitude : 0;
        const Real length = 1 + amplitude;
        const Real measure = dim == 1 ? length : 1;
        const Real weightedMeasure =
          dim == 1 ? length + length * length / 2 : 1 + Real(dim) / 2 + amplitude / 3;
        return alpha * square * phaseSquare * (reaction ? weightedMeasure : measure);
      }

      /**
       * @brief Registers affine and curved operator/load studies at three sizes.
       * @tparam Callback Benchmark implementation callable.
       * @param callback Implementation receiving geometry and workload selections.
       * @param manual Whether aligned manual timing is required for MPI.
       */
      template <class Callback>
      static void registerCases(Callback callback, bool manual)
      {
        using Geometry::Polytope;
        const std::string family = Family == MassFamily::H1 ? "H1/P" + std::to_string(K)
          : Family == MassFamily::P1                        ? "P1"
          : Family == MassFamily::P0                        ? "P0"
                                                            : "P0g";
        for (auto geometry :
          {Polytope::Type::Segment, Polytope::Type::Triangle,
            Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
            Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
        {
          for (bool reaction : {false, true})
          {
            for (bool rhs : {false, true})
            {
              for (bool curved : {false, true})
              {
                const std::string name = std::string(curved ? "CurvedQ2/" : "") +
                  (reaction ? "Reaction/" : "Mass/") + (rhs ? "Load/" : "Operator/") +
                  family + "/" + (std::is_same_v<Scalar, Complex> ? "Complex" : "Real") +
                  (Vector ? "Vector/" : "Scalar/") +
                  std::string(Convergence::UniformGrid::getGeometryName(geometry));
                auto* entry = benchmark::RegisterBenchmark(name.c_str(),
                  [callback, geometry, reaction, rhs, curved](
                    auto& state) { callback(state, geometry, reaction, rhs, curved); })
                                ->Arg(3)
                                ->Arg(5)
                                ->Arg(9);
                if (manual)
                  entry->Iterations(3)->UseManualTime();
                else
                  entry->UseRealTime();
              }
            }
          }
        }
      }
  };
}
#endif
