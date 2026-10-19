/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_H1_H1ELEMENT_HPP
#define RODIN_VARIATIONAL_H1_H1ELEMENT_HPP

/**
 * @file
 * @brief Template implementation for H1Element class.
 *
 * This file contains the template member function implementations for
 * H1Element, including:
 * - Basis function evaluation using modal Dubiner basis
 * - First-order derivative computation via chain rule
 * - Vandermonde matrix inversion for nodal-to-modal conversion
 *
 * The implementation uses:
 * - **Segments**: Classical Lagrange interpolation with GLL nodes
 * - **Triangles**: Dubiner modal basis with Fekete nodes
 * - **Quadrilaterals**: Tensor product of 1D GLL Lagrange bases
 * - **Tetrahedra**: 3D Dubiner modal basis with Fekete nodes
 * - **Wedges**: Product of triangle Fekete basis and 1D GLL basis
 */

#include "Rodin/Math/Common.h"

#include "Rodin/QF/QuadratureFormula.h"
#include "Rodin/Variational/ReferenceTabulationCache.h"

#include "H1Element.h"

#include "LagrangeBasis.h"
#include "Dubiner.h"
#include "Fekete.h"
#include "WarpBlend.h"
#include "GLL.h"
#include "LegendrePolynomial.h"

namespace Rodin::Variational
{
  /// @brief Indexing of the tensor-product modes on the reference pyramid.
  template <size_t K>
  struct PyramidIndex
  {
    /// @brief Number of entries.
    static constexpr size_t Count = (K + 1) * (K + 2) * (2 * K + 3) / 6;

    /// @brief Pair of tensor-product indices of a pyramid mode.
    struct IJ
    {
      /// @brief Index of the cached polytope.
      size_t i;

      /// @brief Second tensor-product index.
      size_t j;
    };

    /**
     * @brief Gets the index offset of a pyramid layer.
     * @returns The index offset of a pyramid layer.
     * @param layer Lattice layer index.
     */
    static constexpr size_t getLayerOffset(size_t layer)
    {
      size_t out = 0;
      for (size_t k = 0; k < layer; ++k)
      {
        const size_t n = K - k + 1;
        out += n * n;
      }
      return out;
    }

    /**
     * @brief Gets the linear index of a tensor-product mode.
     * @param i Index of the requested entry.
     * @param j Index of the second coordinate.
     * @returns The linear index of a tensor-product mode.
     * @param k Index of the requested lattice entry.
     */
    static constexpr size_t getIndex(size_t i, size_t j, size_t k)
    {
      const size_t n = K - k + 1;
      return getLayerOffset(k) + j * n + i;
    }

    /**
     * @brief Decodes a linear index into its tensor-product indices.
     * @param idx Index of the requested entry.
     * @param i Index of the requested entry.
     * @param j Index of the second coordinate.
     * @param k Storage for the decoded lattice index.
     */
    static constexpr void decode(size_t idx, size_t& i, size_t& j, size_t& k)
    {
      size_t rem = idx;
      for (k = 0; k <= K; ++k)
      {
        const size_t n = K - k + 1;
        const size_t layer = n * n;
        if (rem < layer)
        {
          j = rem / n;
          i = rem % n;
          return;
        }
        rem -= layer;
      }

      i = j = k = 0;
    }

    /**
     * @brief Gets the number of modes on a triangular lattice.
     * @returns The number of modes on a triangular lattice.
     * @param alpha Basis or lattice index in the reference element.
     */
    static constexpr IJ getTriangleLattice(size_t alpha)
    {
      size_t pos = 0;
      for (size_t j = 0; j <= K; ++j)
      {
        for (size_t i = 0; i <= K - j; ++i, ++pos)
        {
          if (pos == alpha)
            return IJ{i, j};
        }
      }
      return IJ{0, 0};
    }

    /**
     * @brief Gets the index of a mode on a pyramid side.
     * @param local Index in the local numbering.
     * @returns The index of a mode on a pyramid side.
     * @param alpha Basis or lattice index in the reference element.
     */
    static constexpr size_t getSideIndex(size_t local, size_t alpha)
    {
      const auto ij = getTriangleLattice(alpha);
      const size_t i = ij.i;
      const size_t k = ij.j;
      const size_t n = K - k;

      switch (local)
      {
        case 1:
          return getIndex(i, 0, k);
        case 2:
          return getIndex(n, i, k);
        case 3:
          return getIndex(n - i, n, k);
        case 4:
          return getIndex(0, n - i, k);
        default:
          return 0;
      }
    }
  };

  /// @brief Bernstein basis on the reference pyramid.
  template <size_t K>
  class BernsteinPyramid
  {
    public:
      /**
       * @brief Gets a binomial coefficient.
       * @returns A binomial coefficient.
       * @param n Nonnegative upper argument of the binomial coefficient.
       * @param k Lower argument of the binomial coefficient.
       */
      static Real getBinomial(size_t n, size_t k)
      {
        if (k > n)
          return 0;
        if (k > n - k)
          k = n - k;

        Real out = 1;
        for (size_t i = 1; i <= k; ++i)
        {
          out *= static_cast<Real>(n + 1 - i);
          out /= static_cast<Real>(i);
        }
        return out;
      }

      /**
       * @brief Gets the basis function of a local degree of freedom.
       * @param i Local basis-function index.
       * @returns Value of the selected local basis function at the evaluation point.
       * @param n Basis or lattice index in the reference element.
       * @param x Reference-coordinate component at which the basis is evaluated.
       */
      static Real getBasis(size_t n, size_t i, Real x)
      {
        if (i > n)
          return 0;

        Real out = getBinomial(n, i);
        for (size_t p = 0; p < i; ++p)
          out *= x;
        for (size_t p = 0; p < n - i; ++p)
          out *= (Real(1) - x);
        return out;
      }

      /**
       * @brief Gets the derivative of the basis function.
       * @param i Index of the requested entry.
       * @returns Derivative evaluated at the supplied point.
       * @param n Basis or lattice index in the reference element.
       * @param x Reference-coordinate component at which the basis is evaluated.
       */
      static Real getDerivative(size_t n, size_t i, Real x)
      {
        if (n == 0)
          return 0;

        Real left = 0;
        if (i > 0)
          left = getBasis(n - 1, i - 1, x);

        Real right = 0;
        if (i < n)
          right = getBasis(n - 1, i, x);

        return static_cast<Real>(n) * (left - right);
      }
  };

  /// @brief Modal basis on the reference pyramid.
  template <size_t K>
  class PyramidModal
  {
    public:
      /**
       * @brief Gets the basis function of a local degree of freedom.
       * @param r Point at which the operation is evaluated.
       * @returns Value of the selected local basis function at the evaluation point.
       * @param mode Basis-mode index.
       */
      static Real getBasis(size_t mode, const Math::SpatialPoint& r)
      {
        size_t i, j, k;
        PyramidIndex<K>::decode(mode, i, j, k);

        const Real z = r.z();
        const Real q = Real(1) - z;
        if (q == 0.0)
          return k == K ? Real(1) : Real(0);

        const size_t n = K - k;
        const Real a = r.x() / q;
        const Real b = r.y() / q;

        return BernsteinPyramid<K>::getBasis(n, i, a) *
          BernsteinPyramid<K>::getBasis(n, j, b) * BernsteinPyramid<K>::getBasis(K, k, z);
      }

      /**
       * @brief Gets the derivative of the basis function.
       * @param r Point at which the operation is evaluated.
       * @returns Derivative evaluated at the supplied point.
       * @param mode Basis-mode index.
       * @param deriv Coordinate direction of differentiation.
       */
      static Real getDerivative(size_t mode, size_t deriv, const Math::SpatialPoint& r)
      {
        size_t i, j, k;
        PyramidIndex<K>::decode(mode, i, j, k);

        const Real z = r.z();
        const Real q = Real(1) - z;
        if (q == 0.0)
          return 0;

        const size_t n = K - k;
        const Real a = r.x() / q;
        const Real b = r.y() / q;

        const Real Bz = BernsteinPyramid<K>::getBasis(K, k, z);
        if (deriv == 0)
          return BernsteinPyramid<K>::getDerivative(n, i, a) *
            BernsteinPyramid<K>::getBasis(n, j, b) * Bz / q;
        const Real Ba = BernsteinPyramid<K>::getBasis(n, i, a);
        const Real dBb = BernsteinPyramid<K>::getDerivative(n, j, b);
        if (deriv == 1)
          return Ba * dBb * Bz / q;
        if (deriv == 2)
        {
          const Real Bb = BernsteinPyramid<K>::getBasis(n, j, b);
          const Real dBa = BernsteinPyramid<K>::getDerivative(n, i, a);
          const Real dBz = BernsteinPyramid<K>::getDerivative(K, k, z);
          return (dBa * (a / q) * Bb + Ba * dBb * (b / q)) * Bz + Ba * Bb * dBz;
        }
        return 0;
      }

      /**
       * @brief Fixed modal rows containing a value followed by three derivatives.
       *
       * Row @f$m@f$ contains @f$(\psi_m,\partial_0\psi_m,
       * \partial_1\psi_m,\partial_2\psi_m)@f$. Its compile-time layout
       * retains the original point-local modal workspace.
       */
      using Tabulation = std::array<std::array<Real, 4>, PyramidIndex<K>::Count>;

      /// @brief Selects the value column of caller-owned modal rows.
      struct Values
      {
        Tabulation&
          data; ///< Caller-owned modal rows; derivative columns are untouched.
      };

      /// @brief Selects the three derivative columns of caller-owned modal rows.
      struct Gradients
      {
        Tabulation& data; ///< Caller-owned modal rows; the value column is untouched.
      };

      /**
       * @brief Tabulates all four modal columns with shared point-local work.
       * @param output Caller-owned fixed-layout modal rows.
       * @param rc Three reference coordinates, with optional trailing padding.
       * @throws Alert::Exception If coordinates are invalid.
       */
      static void tabulate(Tabulation& output, const Math::SpatialPoint& rc)
      {
        tabulate<Tabulation&>(output, rc);
      }

      /**
       * @brief Tabulates modal values at one point.
       * @param values Mode-indexed output.
       * @param rc Reference coordinates.
       * @throws Alert::Exception If coordinates are invalid.
       */
      static void tabulate(Values values, const Math::SpatialPoint& rc)
      {
        tabulate<Values>(values, rc);
      }

      /**
       * @brief Tabulates modal reference gradients at one point.
       * @param gradients Mode-major output.
       * @param rc Reference coordinates.
       * @throws Alert::Exception If coordinates are invalid.
       */
      static void tabulate(Gradients gradients, const Math::SpatialPoint& rc)
      {
        tabulate<Gradients>(gradients, rc);
      }

    private:
      /**
       * @brief Computes the requested modal quantities with compile-time selection.
       * @tparam Output Selected modal columns or the complete tabulation reference.
       * @param output Caller-owned modal quantities.
       * @param rc Three reference coordinates.
       * @throws Alert::Exception If reference coordinates are missing.
       */
      template <class Output>
      static void tabulate(Output output, const Math::SpatialPoint& rc)
      {
        constexpr bool HasValues = !std::is_same_v<Output, Gradients>;
        constexpr bool HasGradients = !std::is_same_v<Output, Values>;
        auto& modal = [&]() -> Tabulation& {
          if constexpr (std::is_same_v<Output, Tabulation&>)
            return output;
          else
            return output.data;
        }();
        if (rc.size() < 3)
        {
          Alert::NamespacedException("Rodin::Variational::PyramidModal::tabulate")
            << "Modal tabulation requires three reference coordinates." << Alert::Raise;
          return;
        }
        constexpr size_t Count = PyramidIndex<K>::Count;
        for (size_t mode = 0; mode < Count; ++mode)
        {
          if constexpr (HasValues)
            modal[mode][0] = getBasis(mode, rc);
          if constexpr (HasGradients)
          {
            for (size_t direction = 0; direction < 3; ++direction)
              modal[mode][direction + 1] = getDerivative(mode, direction, rc);
          }
        }
      }
  };

  /// @brief Vandermonde matrix of the modal basis on the reference pyramid.
  template <size_t K>
  class VandermondePyramid
  {
    public:
      /**
       * @brief Gets the underlying matrix.
       * @returns The underlying matrix.
       */
      static const Math::Matrix<Real>& getMatrix()
      {
        static const Math::Matrix<Real> s_vandermonde = [] {
          constexpr size_t N = PyramidIndex<K>::Count;
          const auto& nodes =
            H1Element<K, Real>::getNodes(Geometry::Polytope::Type::Pyramid);
          Math::Matrix<Real> vandermonde(N, N);

          for (size_t node = 0; node < N; ++node)
          {
            for (size_t mode = 0; mode < N; ++mode)
              vandermonde(node, mode) = PyramidModal<K>::getBasis(mode, nodes[node]);
          }
          return vandermonde;
        }();

        return s_vandermonde;
      }

      /**
       * @brief Gets the inverse of the matrix.
       * @returns The inverse of the matrix.
       */
      static const Math::Matrix<Real>& getInverse()
      {
        static const Math::Matrix<Real> s_inv = [] {
          const auto& V = getMatrix();
#if EIGEN_VERSION_AT_LEAST(3, 4, 90)
          Eigen::BDCSVD<Math::Matrix<Real>, Eigen::ComputeThinU | Eigen::ComputeThinV>
            svd(V);
#else
          Eigen::BDCSVD<Math::Matrix<Real>> svd(
            V, Eigen::ComputeThinU | Eigen::ComputeThinV);
#endif
          const Math::Matrix<Real> I = Math::Matrix<Real>::Identity(V.rows(), V.cols());
          return svd.solve(I).eval();
        }();

        return s_inv;
      }
  };

  template <size_t K, class Scalar>
  template <Geometry::Polytope::Type ReferenceGeometry>
  class H1Element<K, Scalar>::Evaluation
  {
      static_assert(
        ReferenceGeometry == G::Tetrahedron || ReferenceGeometry == G::Pyramid);
      using Modal = std::conditional_t<ReferenceGeometry == G::Tetrahedron,
        DubinerTetrahedron<K>, PyramidModal<K>>;
      using Vandermonde = std::conditional_t<ReferenceGeometry == G::Tetrahedron,
        VandermondeTetrahedron<K>, VandermondePyramid<K>>;

    public:
      /**
       * @brief Evaluates shared modal derivatives and obtains the cached inverse once.
       * @param point Reference coordinates with at least three components.
       * @throws Alert::Exception If coordinates are missing.
       */
      explicit Evaluation(const Math::SpatialPoint& point)
        : m_inverse(Vandermonde::getInverse())
      {
        if (point.size() < 3)
          Alert::MemberFunctionException(*this, __func__)
            << "Reference gradient evaluation requires three coordinate components."
            << Alert::Raise;
        // The value column is unused here, but initialized so copies own fully
        // initialized rows without zeroing derivatives that tabulation overwrites.
        for (auto& row : m_modal)
          row[0] = 0;
        Modal::tabulate(typename Modal::Gradients{m_modal}, point);
      }

      /**
       * @brief Returns the number of scalar nodal basis functions.
       * @returns Number of scalar nodal basis functions in this evaluation.
       */
      static constexpr size_t getCount()
      {
        return std::tuple_size_v<typename Modal::Tabulation>;
      }

      /**
       * @brief Converts one modal gradient into the requested nodal gradient.
       *
       * Modes retain their original summation order and real accumulator before
       * promotion. No array of all nodal gradients is formed.
       * @param local Local nodal basis index.
       * @returns The three reference derivatives of the selected basis function.
       * @throws Alert::Exception If the basis index is outside the evaluation.
       */
      Math::SpatialVector<Scalar> getGradient(size_t local) const
      {
        if (local >= getCount())
          Alert::MemberFunctionException(*this, __func__)
            << "Basis index " << local << " is outside the reference evaluation."
            << Alert::Raise;
        Real dx = 0;
        Real dy = 0;
        Real dz = 0;
        for (size_t mode = 0; mode < getCount(); ++mode)
        {
          const Real weight = m_inverse(mode, local);
          dx += weight * m_modal[mode][1];
          dy += weight * m_modal[mode][2];
          dz += weight * m_modal[mode][3];
        }
        return Math::SpatialVector<Scalar>{Scalar(dx), Scalar(dy), Scalar(dz)};
      }

    private:
      typename Modal::Tabulation m_modal;
      const Math::Matrix<Real>& m_inverse;
  };

  /**
   * @brief Prepares shared modal derivatives for nodal gradient evaluation.
   * @tparam ReferenceGeometry Tetrahedron or pyramid reference geometry.
   * @param point Reference coordinates; optional trailing padding is ignored.
   * @returns Owned modal derivatives with shared inverse Vandermonde access.
   * @throws Alert::Exception If the geometry or coordinates are incompatible.
   */
  template <size_t K, class Scalar>
  template <Geometry::Polytope::Type ReferenceGeometry>
  auto H1Element<K, Scalar>::tabulate(
    const Math::SpatialPoint& point) const -> Evaluation<ReferenceGeometry>
  {
    if (this->getGeometry() != ReferenceGeometry || point.size() < 3)
      Alert::MemberFunctionException(*this, __func__)
        << "Reference gradient evaluation requires a matching geometry and three "
           "coordinates."
        << Alert::Raise;
    return Evaluation<ReferenceGeometry>(point);
  }

  template <size_t K, class Scalar>
  Optional<std::reference_wrapper<const typename H1Element<K, Scalar>::Tabulation>>
  H1Element<K, Scalar>::findTabulation(size_t identity) const
  {
    return ReferenceTabulationCache<H1Element<K, Scalar>, Scalar>::find(*this, identity);
  }

  template <size_t K, class Scalar>
  const typename H1Element<K, Scalar>::Tabulation& H1Element<K, Scalar>::getTabulation(
    const QF::QuadratureFormulaBase& qf) const
  {
    return ReferenceTabulationCache<H1Element<K, Scalar>, Scalar>::get(*this, qf);
  }

  template <size_t K, class Scalar>
  void H1Element<K, Scalar>::tabulate(
    Values values, const Math::SpatialPoint& point) const
  {
    tabulate<Values>(values, point);
  }

  template <size_t K, class Scalar>
  void H1Element<K, Scalar>::tabulate(
    Gradients gradients, const Math::SpatialPoint& point) const
  {
    tabulate<Gradients>(gradients, point);
  }

  template <size_t K, class Scalar>
  void H1Element<K, Scalar>::tabulate(
    Values values, Gradients gradients, const Math::SpatialPoint& point) const
  {
    tabulate<std::pair<Values, Gradients>>({values, gradients}, point);
  }

  template <size_t K, class Scalar>
  template <class Output>
  void H1Element<K, Scalar>::tabulate(
    Output output, const Math::SpatialPoint& point) const
  {
    constexpr bool HasValues = !std::is_same_v<Output, Gradients>;
    constexpr bool HasGradients = !std::is_same_v<Output, Values>;
    const auto values = [&]() -> Values {
      if constexpr (std::is_same_v<Output, Values>)
        return output;
      else if constexpr (HasValues)
        return output.first;
      else
        return {{}};
    }();
    const auto gradients = [&]() -> Gradients {
      if constexpr (std::is_same_v<Output, Gradients>)
        return output;
      else if constexpr (HasGradients)
        return output.second;
      else
        return {{}, 0};
    }();
    const size_t dimension =
      Geometry::Polytope::Traits(this->getGeometry()).getDimension();
    const size_t count = this->getCount();
    if (static_cast<size_t>(point.size()) < dimension)
      Alert::MemberFunctionException(*this, __func__)
        << "Reference tabulation requires " << dimension << " coordinate components, got "
        << point.size() << "." << Alert::Raise;
    if ((HasValues && values.data.size() != count) ||
      (HasGradients &&
        (gradients.dimension != dimension || gradients.data.size() != count * dimension)))
      Alert::MemberFunctionException(*this, __func__)
        << "Reference tabulation output does not match the basis count and dimension."
        << Alert::Raise;

    const auto contract = [&]<class Modal, class Vandermonde, size_t Count>() {
      typename Modal::Tabulation modal;
      if constexpr (HasValues && HasGradients)
        Modal::tabulate(modal, point);
      else if constexpr (HasValues)
        Modal::tabulate(typename Modal::Values{modal}, point);
      else
        Modal::tabulate(typename Modal::Gradients{modal}, point);
      const auto& inverse = Vandermonde::getInverse();
      for (size_t local = 0; local < Count; ++local)
      {
        if constexpr (HasValues)
        {
          using Accumulator =
            std::conditional_t<std::is_same_v<Modal, PyramidModal<K>>, Scalar, Real>;
          Accumulator value = 0;
          for (size_t mode = 0; mode < Count; ++mode)
            value += inverse(mode, local) * modal[mode][0];
          values.data[local] = Scalar(value);
        }
        if constexpr (HasGradients)
        {
          Math::SpatialVector<Scalar> gradient(3);
          for (size_t direction = 0; direction < 3; ++direction)
          {
            Real derivative = 0;
            for (size_t mode = 0; mode < Count; ++mode)
              derivative += inverse(mode, local) * modal[mode][direction + 1];
            gradient(direction) = Scalar(derivative);
          }
          for (size_t direction = 0; direction < 3; ++direction)
            gradients.data[3 * local + direction] = gradient(direction);
        }
      }
    };
    if (this->getGeometry() == Geometry::Polytope::Type::Tetrahedron)
    {
      contract.template operator()<DubinerTetrahedron<K>, VandermondeTetrahedron<K>,
        (K + 1) * (K + 2) * (K + 3) / 6>();
      return;
    }
    if (this->getGeometry() == Geometry::Polytope::Type::Pyramid)
    {
      contract.template
      operator()<PyramidModal<K>, VandermondePyramid<K>, PyramidIndex<K>::Count>();
      return;
    }
    for (size_t local = 0; local < count; ++local)
    {
      const auto basis = this->getBasis(local);
      if constexpr (HasValues)
        values.data[local] = basis(point);
      if constexpr (HasGradients)
      {
        Math::SpatialVector<Scalar> gradient(dimension);
        for (size_t direction = 0; direction < dimension; ++direction)
          gradient(direction) = basis.template getDerivative<1>(direction)(point);
        for (size_t direction = 0; direction < dimension; ++direction)
          gradients.data[local * dimension + direction] = gradient(direction);
      }
    }
  }

  template <size_t K, class Scalar>
  void H1Element<K, Scalar>::tabulate(
    Tabulation& table, const QF::QuadratureFormulaBase& formula) const
  {
    table.resize(formula.getSize(), this->getCount(),
      Geometry::Polytope::Traits(this->getGeometry()).getDimension());
    if (table.getSize() == 0)
      return;
    if (this->getGeometry() != Geometry::Polytope::Type::Tetrahedron &&
      this->getGeometry() != Geometry::Polytope::Type::Pyramid)
    {
      std::vector<std::reference_wrapper<const BasisFunction>> bases;
      bases.reserve(table.getCount());
      for (size_t local = 0; local < table.getCount(); ++local)
        bases.emplace_back(this->getBasis(local));
      for (size_t qp = 0; qp < table.getSize(); ++qp)
      {
        const auto& point = formula.getPoint(qp);
        if (static_cast<size_t>(point.size()) < table.getDimension())
          Alert::MemberFunctionException(*this, __func__)
            << "Reference tabulation requires " << table.getDimension()
            << " coordinate components, got " << point.size() << "." << Alert::Raise;
        for (size_t local = 0; local < table.getCount(); ++local)
        {
          const auto& basis = bases[local].get();
          table.getValues(qp).data[local] = basis(point);
          for (size_t direction = 0; direction < table.getDimension(); ++direction)
          {
            table.getGradients(qp).data[local * table.getDimension() + direction] =
              basis.template getDerivative<1>(direction)(point);
          }
        }
      }
      return;
    }
    const auto tabulate =
      [&]<class Modal, class Vandermonde, class Accumulator, size_t Count>() {
        assert(table.getCount() == Count && table.getDimension() == 3);
        const auto& inverse = Vandermonde::getInverse();
        typename Modal::Tabulation modal;
        for (size_t qp = 0; qp < table.getSize(); ++qp)
        {
          Modal::tabulate(modal, formula.getPoint(qp));
          for (size_t local = 0; local < Count; ++local)
          {
            Accumulator value = 0;
            for (size_t mode = 0; mode < Count; ++mode)
              value += inverse(mode, local) * modal[mode][0];
            table.getValues(qp).data[local] = Scalar(value);
            for (size_t direction = 0; direction < 3; ++direction)
            {
              Real derivative = 0;
              for (size_t mode = 0; mode < Count; ++mode)
                derivative += inverse(mode, local) * modal[mode][direction + 1];
              table.getGradients(qp).data[local * table.getDimension() + direction] =
                Scalar(derivative);
            }
          }
        }
      };
    if (this->getGeometry() == Geometry::Polytope::Type::Tetrahedron)
      tabulate.template operator()<DubinerTetrahedron<K>, VandermondeTetrahedron<K>, Real,
        FeketeTetrahedron<K>::Count>();
    else
      tabulate.template operator()<PyramidModal<K>, VandermondePyramid<K>, Scalar,
        PyramidIndex<K>::Count>();
  }

  template <size_t K, class Scalar>
  constexpr const Math::SpatialPoint& H1Element<K, Scalar>::getNode(size_t i) const
  {
    return this->getNodes(this->getGeometry())[i];
  }

  template <size_t K, class Scalar>
  const typename H1Element<K, Scalar>::LinearForm& H1Element<K, Scalar>::getLinearForm(
    size_t i) const
  {
    const Geometry::Polytope::Type g = this->getGeometry();

    // These objects depend only on the compile-time order and reference
    // geometry. Sharing their immutable storage avoids per-thread replicas.
    const auto makeLinearForms = [this, g] {
      std::vector<LinearForm> forms;
      const size_t count = getCount();
      forms.reserve(count);
      for (size_t j = 0; j < count; ++j)
        forms.emplace_back(j, g);
      return forms;
    };

    switch (g)
    {
      case Geometry::Polytope::Type::Point:
      {
        static const std::vector<LinearForm> s_lfs = makeLinearForms();
        return s_lfs[i];
      }
      case Geometry::Polytope::Type::Segment:
      {
        static const std::vector<LinearForm> s_lfs = makeLinearForms();
        return s_lfs[i];
      }
      case Geometry::Polytope::Type::Triangle:
      {
        static const std::vector<LinearForm> s_lfs = makeLinearForms();
        return s_lfs[i];
      }
      case Geometry::Polytope::Type::Quadrilateral:
      {
        static const std::vector<LinearForm> s_lfs = makeLinearForms();
        return s_lfs[i];
      }
      case Geometry::Polytope::Type::Tetrahedron:
      {
        static const std::vector<LinearForm> s_lfs = makeLinearForms();
        return s_lfs[i];
      }
      case Geometry::Polytope::Type::Pyramid:
      {
        static const std::vector<LinearForm> s_lfs = makeLinearForms();
        return s_lfs[i];
      }
      case Geometry::Polytope::Type::Wedge:
      {
        static const std::vector<LinearForm> s_lfs = makeLinearForms();
        return s_lfs[i];
      }
      case Geometry::Polytope::Type::Hexahedron:
      {
        static const std::vector<LinearForm> s_lfs = makeLinearForms();
        return s_lfs[i];
      }
    }

    // Fallback (should never happen)
    static const LinearForm s_null(0, Geometry::Polytope::Type::Point);
    assert(false);
    return s_null;
  }

  template <size_t K, class Scalar>
  const typename H1Element<K, Scalar>::BasisFunction& H1Element<K, Scalar>::getBasis(
    size_t i) const
  {
    const Geometry::Polytope::Type g = this->getGeometry();

    // Basis descriptors are immutable after construction and are safe to share
    // across concurrent element evaluations.
    const auto makeBasis = [this, g] {
      std::vector<BasisFunction> basis;
      const size_t count = getCount();
      basis.reserve(count);
      for (size_t j = 0; j < count; ++j)
        basis.emplace_back(j, g);
      return basis;
    };

    switch (g)
    {
      case Geometry::Polytope::Type::Point:
      {
        static const std::vector<BasisFunction> s_bs = makeBasis();
        return s_bs[i];
      }
      case Geometry::Polytope::Type::Segment:
      {
        static const std::vector<BasisFunction> s_bs = makeBasis();
        return s_bs[i];
      }
      case Geometry::Polytope::Type::Triangle:
      {
        static const std::vector<BasisFunction> s_bs = makeBasis();
        return s_bs[i];
      }
      case Geometry::Polytope::Type::Quadrilateral:
      {
        static const std::vector<BasisFunction> s_bs = makeBasis();
        return s_bs[i];
      }
      case Geometry::Polytope::Type::Tetrahedron:
      {
        static const std::vector<BasisFunction> s_bs = makeBasis();
        return s_bs[i];
      }
      case Geometry::Polytope::Type::Pyramid:
      {
        static const std::vector<BasisFunction> s_bs = makeBasis();
        return s_bs[i];
      }
      case Geometry::Polytope::Type::Wedge:
      {
        static const std::vector<BasisFunction> s_bs = makeBasis();
        return s_bs[i];
      }
      case Geometry::Polytope::Type::Hexahedron:
      {
        static const std::vector<BasisFunction> s_bs = makeBasis();
        return s_bs[i];
      }
    }

    // Fallback (should never happen)
    static const BasisFunction s_null(0, Geometry::Polytope::Type::Point);
    assert(false);
    return s_null;
  }

  template <size_t K, class Scalar>
  Scalar H1Element<K, Scalar>::BasisFunction::operator()(
    const Math::SpatialPoint& r) const
  {
    switch (m_g)
    {
      case Geometry::Polytope::Type::Point:
      {
        return Scalar(1);
      }
      case Geometry::Polytope::Type::Segment:
      {
        // Canonical [0,1] segment with GLL01 nodes
        return LagrangeBasisSegment<K>::getBasis(m_local, r.value());
      }
      case Geometry::Polytope::Type::Triangle:
      {
        // Modal Dubiner basis + Vandermonde on Fekete nodes
        const auto& inverse = VandermondeTriangle<K>::getInverse();

        Real rc, sc;
        DubinerTriangle<K>::getCollapsed(rc, sc, r.x(), r.y());

        Scalar result = Scalar(0);
        size_t modeIdx = 0;

        Rodin::Utility::ForIndex<K + 1>([&](auto pIdx) {
          constexpr size_t P = pIdx.value;
          Rodin::Utility::ForIndex<K + 1 - P>([&](auto qIdx) {
            constexpr size_t Q = qIdx.value;
            Real psi;
            DubinerTriangle<K>::template getBasis<P, Q>(psi, rc, sc);
            result += inverse(modeIdx, m_local) * psi;
            ++modeIdx;
          });
        });

        return result;
      }
      case Geometry::Polytope::Type::Quadrilateral:
      {
        // Tensor product GLL01 × GLL01 (Lagrange)
        const size_t jIdx = m_local / (K + 1);
        const size_t iIdx = m_local % (K + 1);
        return LagrangeBasisQuadrilateral<K>::getBasis(iIdx, jIdx, r.x(), r.y());
      }
      case Geometry::Polytope::Type::Tetrahedron:
      {
        // Modal Dubiner basis + Vandermonde on tetra Fekete nodes
        const auto& inverse = VandermondeTetrahedron<K>::getInverse();

        Real ac, bc, cc;
        DubinerTetrahedron<K>::getCollapsed(ac, bc, cc, r.x(), r.y(), r.z());

        Real result = 0;
        size_t modeIdx = 0;

        Rodin::Utility::ForIndex<K + 1>([&](auto pIdx) {
          constexpr size_t P = pIdx.value;
          Rodin::Utility::ForIndex<K + 1 - P>([&](auto qIdx) {
            constexpr size_t Q = qIdx.value;
            Rodin::Utility::ForIndex<K + 1 - P - Q>([&](auto rIdx) {
              constexpr size_t R = rIdx.value;
              Real psi;
              DubinerTetrahedron<K>::template getBasis<P, Q, R>(psi, ac, bc, cc);
              result += inverse(modeIdx, m_local) * psi;
              ++modeIdx;
            });
          });
        });

        return result;
      }
      case Geometry::Polytope::Type::Pyramid:
      {
        const auto& inverse = VandermondePyramid<K>::getInverse();
        constexpr size_t count = PyramidIndex<K>::Count;

        Scalar result = Scalar(0);
        for (size_t mode = 0; mode < count; ++mode)
          result += inverse(mode, m_local) * PyramidModal<K>::getBasis(mode, r);

        return result;
      }
      case Geometry::Polytope::Type::Wedge:
      {
        // Number of triangle DOFs (Fekete triangle)
        constexpr size_t ntri =
          FeketeTriangle<K>::Count; // same as FeketeTriangle<K>::Count

        const size_t k = m_local / ntri; // z-index, 0..K
        const size_t alpha = m_local % ntri; // triangle node index

        // --- triangle factor: nodal basis on Fekete nodes ---
        const auto& Vinv = VandermondeTriangle<K>::getInverse();

        Real rc, sc;
        DubinerTriangle<K>::getCollapsed(rc, sc, r.x(), r.y());

        Scalar triVal = Scalar(0);
        size_t modeIdx = 0;

        Rodin::Utility::ForIndex<K + 1>([&](auto pIdx) {
          constexpr size_t P = pIdx.value;
          Rodin::Utility::ForIndex<K + 1 - P>([&](auto qIdx) {
            constexpr size_t Q = qIdx.value;
            Real psi;
            DubinerTriangle<K>::template getBasis<P, Q>(psi, rc, sc);
            triVal += Vinv(modeIdx, alpha) * psi;
            ++modeIdx;
          });
        });

        // --- segment factor in z (GLL01, Lagrange) ---
        const Real z = r(2);
        const Real segVal = LagrangeBasisSegment<K>::getBasis(k, z);

        return triVal * segVal;
      }
      case Geometry::Polytope::Type::Hexahedron:
      {
        // Tensor product GLL01 × GLL01 × GLL01
        // Node ordering matches getNodes(): i + (K+1)*(j + (K+1)*k)
        const size_t n1 = K + 1;
        const size_t k = m_local / (n1 * n1);
        const size_t r2 = m_local % (n1 * n1);
        const size_t j = r2 / n1;
        const size_t i = r2 % n1;

        const Real x = r.x();
        const Real y = r.y();
        const Real z = r.z();

        const Real lx = LagrangeBasisSegment<K>::getBasis(i, x);
        const Real ly = LagrangeBasisSegment<K>::getBasis(j, y);
        const Real lz = LagrangeBasisSegment<K>::getBasis(k, z);

        return static_cast<Scalar>(lx * ly * lz);
      }
    }

    return Math::nan<Scalar>();
  }

  template <size_t K, class Scalar>
  template <size_t Order>
  Scalar H1Element<K, Scalar>::BasisFunction::DerivativeFunction<Order>::operator()(
    const Math::SpatialPoint& r) const
  {
    if constexpr (!std::is_same_v<Scalar, Real>)
    {
      // Reference polynomials and their coordinates are real even when the
      // coefficients of the finite element field are complex.
      return typename H1Element<K, Real>::BasisFunction(m_local, m_g)
        .template getDerivative<Order>(m_i)(r);
    }
    else if constexpr (Order == 0)
    {
      return BasisFunction(m_local, m_g)(r);
    }
    else if constexpr (Order == 1)
    {
      switch (m_g)
      {
        case Geometry::Polytope::Type::Point:
        {
          return Scalar(0);
        }

        case Geometry::Polytope::Type::Segment:
        {
          // \partialφ_i / \partialx on [0,1] with GLL01 nodes
          assert(m_i == 0);
          return LagrangeBasisSegment<K>::getDerivative(m_local, r.value());
        }

        case Geometry::Polytope::Type::Triangle:
        {
          // Dubiner modal gradients in reference coordinates.
          const auto& Vinv = VandermondeTriangle<K>::getInverse();
          Real a, b;
          DubinerTriangle<K>::getCollapsed(a, b, r.x(), r.y());
          const Real s = 1.0 - r.y();

          Scalar result = Scalar(0);
          size_t modeIdx = 0;

          Rodin::Utility::ForIndex<K + 1>([&](auto pIdx) {
            constexpr size_t P = pIdx.value;
            Rodin::Utility::ForIndex<K + 1 - P>([&](auto qIdx) {
              constexpr size_t Q = qIdx.value;

              Real dpsi_dx = 0.0, dpsi_dy = 0.0;
              DubinerTriangle<K>::template getReferenceGradientFromCollapsed<P, Q>(
                dpsi_dx, dpsi_dy, a, b, s);

              if (m_i == 0) // \partial/\partialx
                result += Vinv(modeIdx, m_local) * dpsi_dx;
              else if (m_i == 1) // \partial/\partialy
                result += Vinv(modeIdx, m_local) * dpsi_dy;

              ++modeIdx;
            });
          });

          return result;
        }

        case Geometry::Polytope::Type::Quadrilateral:
        {
          // Tensor product Lagrange on GLL01 × GLL01
          const size_t jIdx = m_local / (K + 1);
          const size_t iIdx = m_local % (K + 1);

          if (m_i == 0) // \partial/\partialx
          {
            return LagrangeBasisQuadrilateral<K>::getDerivative(
              iIdx, jIdx, 0, r.x(), r.y());
          }
          else if (m_i == 1) // \partial/\partialy
          {
            return LagrangeBasisQuadrilateral<K>::getDerivative(
              iIdx, jIdx, 1, r.x(), r.y());
          }
          return Scalar(0);
        }

        case Geometry::Polytope::Type::Tetrahedron:
        {
          // Dubiner modal gradients in reference coordinates.
          const auto& Vinv = VandermondeTetrahedron<K>::getInverse();
          Real a, b, c;
          DubinerTetrahedron<K>::getCollapsed(a, b, c, r.x(), r.y(), r.z());
          const Real s1 = 1.0 - r.y() - r.z();
          const Real s2 = 1.0 - r.z();

          Scalar result = Scalar(0);
          size_t modeIdx = 0;

          Rodin::Utility::ForIndex<K + 1>([&](auto pIdx) {
            constexpr size_t P = pIdx.value;
            Rodin::Utility::ForIndex<K + 1 - P>([&](auto qIdx) {
              constexpr size_t Q = qIdx.value;
              Rodin::Utility::ForIndex<K + 1 - P - Q>([&](auto rIdx) {
                constexpr size_t R = rIdx.value;

                Real dpsi_dx = 0.0;
                Real dpsi_dy = 0.0;
                Real dpsi_dz = 0.0;
                DubinerTetrahedron<K>::template getReferenceGradientFromCollapsed<P, Q,
                  R>(dpsi_dx, dpsi_dy, dpsi_dz, a, b, c, s1, s2);

                if (m_i == 0) // \partial/\partialx
                  result += Vinv(modeIdx, m_local) * dpsi_dx;
                else if (m_i == 1) // \partial/\partialy
                  result += Vinv(modeIdx, m_local) * dpsi_dy;
                else if (m_i == 2) // \partial/\partialz
                  result += Vinv(modeIdx, m_local) * dpsi_dz;

                ++modeIdx;
              });
            });
          });

          return result;
        }
        case Geometry::Polytope::Type::Pyramid:
        {
          const auto& inverse = VandermondePyramid<K>::getInverse();
          constexpr size_t count = PyramidIndex<K>::Count;

          Scalar result = Scalar(0);
          for (size_t mode = 0; mode < count; ++mode)
          {
            result +=
              inverse(mode, m_local) * PyramidModal<K>::getDerivative(mode, m_i, r);
          }

          return result;
        }

        case Geometry::Polytope::Type::Wedge:
        {
          constexpr size_t ntri = FeketeTriangle<K>::Count;

          const size_t k = m_local / ntri;
          const size_t alpha = m_local % ntri;

          const Real z = r(2);

          if (m_i < 2) // \partial/\partialx or \partial/\partialy
          {
            // --- triangle gradient in reference coordinates ---
            const auto& Vinv = VandermondeTriangle<K>::getInverse();
            Real a, b;
            DubinerTriangle<K>::getCollapsed(a, b, r.x(), r.y());
            const Real s = 1.0 - r.y();

            Scalar triDeriv = Scalar(0);
            size_t modeIdx = 0;

            Rodin::Utility::ForIndex<K + 1>([&](auto pIdx) {
              constexpr size_t P = pIdx.value;
              Rodin::Utility::ForIndex<K + 1 - P>([&](auto qIdx) {
                constexpr size_t Q = qIdx.value;

                Real dpsi_dx = 0.0, dpsi_dy = 0.0;
                DubinerTriangle<K>::template getReferenceGradientFromCollapsed<P, Q>(
                  dpsi_dx, dpsi_dy, a, b, s);

                if (m_i == 0) // \partial/\partialx
                  triDeriv += Vinv(modeIdx, alpha) * dpsi_dx;
                else if (m_i == 1) // \partial/\partialy
                  triDeriv += Vinv(modeIdx, alpha) * dpsi_dy;

                ++modeIdx;
              });
            });

            // --- segment value in z ---
            const Real segVal = LagrangeBasisSegment<K>::getBasis(k, z);
            return triDeriv * segVal;
          }
          else // m_i == 2 → \partial/\partialz
          {
            // --- triangle value (same as in BasisFunction wedge case) ---
            const auto& Vinv = VandermondeTriangle<K>::getInverse();

            Real rc, sc;
            DubinerTriangle<K>::getCollapsed(rc, sc, r.x(), r.y());

            Scalar triVal = Scalar(0);
            size_t modeIdx = 0;

            Rodin::Utility::ForIndex<K + 1>([&](auto pIdx) {
              constexpr size_t P = pIdx.value;
              Rodin::Utility::ForIndex<K + 1 - P>([&](auto qIdx) {
                constexpr size_t Q = qIdx.value;
                Real psi;
                DubinerTriangle<K>::template getBasis<P, Q>(psi, rc, sc);
                triVal += Vinv(modeIdx, alpha) * psi;
                ++modeIdx;
              });
            });

            // --- 1D derivative in z ---
            const Real dseg = LagrangeBasisSegment<K>::getDerivative(k, z);
            return triVal * dseg;
          }
        }
        case Geometry::Polytope::Type::Hexahedron:
        {
          // Tensor product derivative on [0,1]^3 with GLL nodes.
          // Node ordering: i + (K+1)*(j + (K+1)*k)
          const size_t n1 = K + 1;
          const size_t k = m_local / (n1 * n1);
          const size_t r2 = m_local % (n1 * n1);
          const size_t j = r2 / n1;
          const size_t i = r2 % n1;

          const Real x = r.x();
          const Real y = r.y();
          const Real z = r.z();

          const Real lx = LagrangeBasisSegment<K>::getBasis(i, x);
          const Real ly = LagrangeBasisSegment<K>::getBasis(j, y);
          const Real lz = LagrangeBasisSegment<K>::getBasis(k, z);

          const Real dlx = LagrangeBasisSegment<K>::getDerivative(i, x);
          const Real dly = LagrangeBasisSegment<K>::getDerivative(j, y);
          const Real dlz = LagrangeBasisSegment<K>::getDerivative(k, z);

          Real val = 0;
          if (m_i == 0) // \partial/\partialx
            val = dlx * ly * lz;
          else if (m_i == 1) // \partial/\partialy
            val = lx * dly * lz;
          else if (m_i == 2) // \partial/\partialz
            val = lx * ly * dlz;

          return static_cast<Scalar>(val);
        }
      }

      return Scalar(0);
    }
    else
    {
      // Higher-order derivatives not implemented
      return Scalar(0);
    }
  }
}

#endif
