/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_TESTS_BENCHMARKS_MATRIXMASSWORKLOAD_H
#define RODIN_TESTS_BENCHMARKS_MATRIXMASSWORKLOAD_H

#include "MassWorkload.h"

namespace Rodin::Tests::Benchmarks
{
  /**
   * @brief Defines matrix-field reaction and constant-field projection forms.
   *
   * The forms are @f$a_L(U,V)=\int_\Omega(\rho U):\overline V\,dx@f$
   * and @f$a_R(U,V)=\int_\Omega(U\rho):\overline V\,dx@f$, with
   * @f$\rho=s\alpha(1+\sum_jx_j)@f$. Over the reals @f$\alpha=1@f$;
   * over the complex numbers @f$\alpha=1+i/2@f$. A representable constant
   * @f$C_{rc}=\gamma(1+r n_c+c)@f$, with @f$\gamma=1@f$ or
   * @f$\gamma=1+2i@f$. Let @f$\mathbf c@f$ denote its DOF vector.
   * Then @f$\mathbf b=A\mathbf c@f$ and
   * @f$\mathbf c^*A\mathbf c=s\alpha\|C\|_F^2\int_\Omega(1+\sum_jx_j)\,dx@f$.
   *
   * Geometry construction and validation are inherited from the scalar
   * workload; the analytic image-of-unit-box moments are used unchanged.
   * Only field construction and the forms change. Neither DOF
   * numbering nor nodal interpolation is assumed: constants are assigned
   * through the public grid-function interpolation contract.
   *
   * @tparam K H1 degree.
   * @tparam Scalar Real or complex field scalar.
   * @tparam Rows Matrix row count.
   * @tparam Columns Matrix column count.
   * @tparam Left Whether the scalar coefficient precedes the trial field.
   */
  template <size_t K, class Scalar, size_t Rows, size_t Columns, bool Left>
  class MatrixMassWorkload : public MassWorkload<K, Scalar, false, MassFamily::H1>
  {
    public:
      using Base = MassWorkload<K, Scalar, false, MassFamily::H1>;

      /**
       * @brief Constructs the matrix space through its public constructor.
       * @tparam Mesh Mesh type, without an assembly-backend policy.
       * @param mesh Physical mesh.
       * @param dim Domain dimension, unused because the field shape is independent.
       * @returns Matrix-valued H1 space.
       */
      template <class Mesh>
      static auto makeSpace(const Mesh& mesh, [[maybe_unused]] size_t dim)
      {
        return Variational::H1<K, Math::SpatialMatrix<Scalar>, Mesh>(
          std::integral_constant<size_t, K>{}, mesh, Rows, Columns);
      }

      /**
       * @brief Configures matrix-field forms and the independent constant action.
       * @tparam Reference Whether the original entry loop supplies the operator.
       * @tparam U Trial-function type.
       * @tparam V Test-function type.
       * @tparam A Bilinear-form type.
       * @tparam B Linear-form type.
       * @param u Trial function receiving the interpolated constant matrix.
       * @param v Test function.
       * @param matrix Public operator form.
       * @param load Public projection-load form.
       * @param dim Domain dimension.
       * @param reaction Whether a variable scalar coefficient is requested.
       * @param scale Borrowed coefficient multiplier.
       * @param curved Whether quadratic geometry is installed.
       * @returns Analytic constant-field action at unit scale.
       */
      template <bool Reference = false, class U, class V, class A, class B>
      static Scalar configure(U& u, V& v, A& matrix, B& load, size_t dim, bool reaction,
        const Real& scale, bool curved)
      {
        using namespace Variational;
        const Scalar phase = [] {
          if constexpr (std::is_same_v<Scalar, Complex>)
            return Scalar(1, 2);
          else
            return Scalar(1);
        }();
        Scalar alpha = 1;
        if constexpr (std::is_same_v<Scalar, Complex>)
        {
          if (reaction)
            alpha = Scalar(1, 0.5);
        }
        Math::SpatialMatrix<Scalar> constant(Rows, Columns);
        for (size_t r = 0; r < Rows; ++r)
        {
          for (size_t c = 0; c < Columns; ++c)
            constant(r, c) = phase * Scalar(1 + r * Columns + c);
        }
        const auto field = MatrixFunction(constant);
        const auto evaluate = [dim, reaction, alpha, &scale](const Geometry::Point& p) {
          Real rho = 1;
          if (reaction)
          {
            for (size_t j = 0; j < dim; ++j)
              rho += p(j);
          }
          return scale * alpha * rho;
        };
        const auto coefficient = [&]() {
          if constexpr (std::is_same_v<Scalar, Complex>)
            return ComplexFunction(evaluate);
          else
            return RealFunction(evaluate);
        }();
        auto integral = [&]() {
          if constexpr (Left)
            return Integral(coefficient * u, v);
          else
            return Integral(u * coefficient, v);
        }();
        integral.setOrder(8);
        if constexpr (Reference)
          matrix = ReferenceIntegral(integral);
        else
          matrix = integral;
        auto rhs = Integral(field * coefficient, v);
        rhs.setOrder(8);
        load = rhs;
        u.getSolution() = field;
        constexpr size_t Components = Rows * Columns;
        const Real square =
          Real(Components * (Components + 1) * (2 * Components + 1)) / 6;
        const Real phaseSquare = std::is_same_v<Scalar, Complex> ? Real(5) : Real(1);
        const Real amplitude = curved ? Base::GeometryAmplitude : 0;
        const Real length = 1 + amplitude;
        const Real measure = dim == 1 ? length : 1;
        const Real weightedMeasure =
          dim == 1 ? length + length * length / 2 : 1 + Real(dim) / 2 + amplitude / 3;
        return alpha * square * phaseSquare * (reaction ? weightedMeasure : measure);
      }
  };
}

#endif
