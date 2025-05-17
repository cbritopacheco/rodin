/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_PK_PKELEMENT_H
#define RODIN_VARIATIONAL_PK_PKELEMENT_H

/**
 * @ingroup RodinDirectives
 * @brief Indicates the maximum vector dimension a PkElement
 */
#define RODIN_PK_MAX_VECTOR_DIMENSION 16

#include <Eigen/Dense>
#include <vector>
#include <array>
#include <cassert>
#include <cmath>

#include "ForwardDecls.h"
#include "Rodin/Types.h"
#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/Connectivity.h"
#include "Rodin/Geometry/GeometryIndexed.h"
#include "Rodin/Variational/FiniteElement.h"
#include "Rodin/Variational/FiniteElementSpace.h"

namespace Rodin::FormLanguage
{
  /**
   * @ingroup TraitsSpecializations
   */
  template <class Range>
  struct Traits<Variational::PkElement<Range>>
  {
    using RangeType = Range;
  };
}

namespace Rodin::Variational
{
  // ─────────────────────────────────────────────────────────────
  //  Internal helper functions (for nodal generation, monomials, etc.)
  // ─────────────────────────────────────────────────────────────
  namespace Internal
  {
    // Generate equispaced nodes on the reference polytope for a given degree.
    // (For brevity we implement only for Point, Segment, and Triangle.)
    inline Math::PointMatrix generateNodes(Geometry::Polytope::Type g, size_t degree)
    {
      if (g == Geometry::Polytope::Type::Point)
      {
        Math::PointMatrix nodes(1, 1);
        nodes(0, 0) = 0;
        return nodes;
      }
      else if (g == Geometry::Polytope::Type::Segment)
      {
        size_t numNodes = degree + 1;
        Math::PointMatrix nodes(1, numNodes);
        for (size_t i = 0; i < numNodes; i++)
          nodes(0, i) = static_cast<Real>(i) / degree;
        return nodes;
      }
      else if (g == Geometry::Polytope::Type::Triangle)
      {
        size_t numNodes = (degree + 1) * (degree + 2) / 2;
        Math::PointMatrix nodes(2, numNodes);
        size_t index = 0;
        for (size_t i = 0; i <= degree; i++)
        {
          for (size_t j = 0; j <= degree - i; j++)
          {
            nodes(0, index) = static_cast<Real>(i) / degree;
            nodes(1, index) = static_cast<Real>(j) / degree;
            index++;
          }
        }
        return nodes;
      }
      else
      {
        assert(false && "generateNodes not implemented for this polytope type.");
        return Math::PointMatrix();
      }
    }

    // Recursively generate multi-indices for monomials in dimension d
    // with total degree <= degree.
    inline void generateMultiIndices(size_t d, size_t degree,
                                     std::vector<std::vector<size_t>>& multiIndices,
                                     std::vector<size_t>& current, size_t pos, size_t remaining)
    {
      if (pos == d)
      {
        if (remaining == 0)
          multiIndices.push_back(current);
        return;
      }
      for (size_t i = 0; i <= remaining; i++)
      {
        current[pos] = i;
        generateMultiIndices(d, degree, multiIndices, current, pos + 1, remaining - i);
      }
    }

    inline std::vector<std::vector<size_t>> generateMultiIndices(size_t d, size_t degree)
    {
      std::vector<std::vector<size_t>> multiIndices;
      std::vector<size_t> current(d, 0);
      generateMultiIndices(d, degree, multiIndices, current, 0, degree);
      return multiIndices;
    }

    // Evaluate a monomial with given exponent vector at point r.
    inline Real evaluateMonomial(const Math::SpatialVector<Real>& r, const std::vector<size_t>& exponents)
    {
      Real prod = 1;
      for (size_t i = 0; i < exponents.size(); i++)
        prod *= std::pow(r[i], static_cast<Real>(exponents[i]));
      return prod;
    }

    // Evaluate all monomials (specified by multiIndices) at point r.
    inline std::vector<Real> evaluateMonomials(const Math::SpatialVector<Real>& r,
                                                 const std::vector<std::vector<size_t>>& multiIndices)
    {
      std::vector<Real> values;
      values.reserve(multiIndices.size());
      for (const auto& exponents : multiIndices)
        values.push_back(evaluateMonomial(r, exponents));
      return values;
    }

    // Build the Vandermonde matrix.
    // Rows: evaluation at each node; Columns: each monomial (ordered as in multiIndices).
    inline Math::Matrix<Real> buildVandermonde(const Math::PointMatrix& nodes,
                                                const std::vector<std::vector<size_t>>& multiIndices)
    {
      size_t numNodes = nodes.cols();
      size_t numMonomials = multiIndices.size();
      Math::Matrix<Real> V(numNodes, numMonomials);
      for (size_t j = 0; j < numNodes; j++)
      {
        Math::SpatialVector<Real> r(nodes.rows());
        for (size_t i = 0; i < nodes.rows(); i++)
          r[i] = nodes(i, j);
        auto mons = evaluateMonomials(r, multiIndices);
        for (size_t i = 0; i < numMonomials; i++)
          V(j, i) = mons[i];
      }
      return V;
    }

    // Invert a matrix using Eigen's built-in .inverse() method.
    inline Math::Matrix<Real> invertMatrix(const Math::Matrix<Real>& A)
    {
      return A.inverse();
    }
  } // namespace Internal

  // ─────────────────────────────────────────────────────────────
  //  PkElement<Real, Degree>: degree-k scalar Lagrange element (Real)
  // ─────────────────────────────────────────────────────────────
  template <size_t Degree>
  class PkElement<Real, Degree> final : public FiniteElementBase<PkElement<Real, Degree>>
  {
    using G = Geometry::Polytope::Type;
  public:
    using Parent = FiniteElementBase<PkElement<Real, Degree>>;
    using RangeType = Real;

    // Linear forms (evaluation at nodes)
    class LinearForm
    {
    public:
      constexpr LinearForm(size_t i, Geometry::Polytope::Type g)
        : m_i(i), m_g(g) {}
      template <class T>
      constexpr auto operator()(const T& v) const { return v(s_nodes[m_g].col(m_i)); }
    private:
      size_t m_i;
      Geometry::Polytope::Type m_g;
    };

    // Basis function stored in the monomial basis (its coefficients)
    class BasisFunction
    {
    public:
      BasisFunction(const std::vector<Real>& coeffs,
                    const std::vector<std::vector<size_t>>& multiIndices)
        : m_coeffs(coeffs), m_multiIndices(multiIndices) {}
      Real operator()(const Math::SpatialVector<Real>& r) const
      {
        auto mons = Internal::evaluateMonomials(r, m_multiIndices);
        Real value = 0;
        for (size_t i = 0; i < m_coeffs.size(); i++)
          value += m_coeffs[i] * mons[i];
        return value;
      }
    private:
      std::vector<Real> m_coeffs;
      std::vector<std::vector<size_t>> m_multiIndices;
    };

    // Gradient function computed by differentiating the monomial basis.
    class GradientFunction
    {
    public:
      GradientFunction(const std::vector<Real>& coeffs,
                       const std::vector<std::vector<size_t>>& multiIndices)
        : m_coeffs(coeffs), m_multiIndices(multiIndices) {}
      Math::SpatialVector<Real> operator()(const Math::SpatialVector<Real>& r) const
      {
        size_t dim = r.size();
        Math::SpatialVector<Real> grad = Math::SpatialVector<Real>::Zero(dim);
        // For each monomial term
        for (size_t i = 0; i < m_coeffs.size(); i++)
        {
          for (size_t d = 0; d < dim; d++)
          {
            if (m_multiIndices[i][d] > 0)
            {
              std::vector<size_t> exp = m_multiIndices[i];
              exp[d]--; // differentiate with respect to direction d
              grad[d] += m_coeffs[i] * m_multiIndices[i][d] * Internal::evaluateMonomial(r, exp);
            }
          }
        }
        return grad;
      }
    private:
      std::vector<Real> m_coeffs;
      std::vector<std::vector<size_t>> m_multiIndices;
    };

    constexpr PkElement() = default;
    constexpr PkElement(Geometry::Polytope::Type geometry)
      : Parent(geometry) {}
    constexpr PkElement(const PkElement& other) : Parent(other) {}
    constexpr PkElement(PkElement&& other) : Parent(std::move(other)) {}
    constexpr PkElement& operator=(const PkElement& other) { Parent::operator=(other); return *this; }
    constexpr PkElement& operator=(PkElement&& other) { Parent::operator=(std::move(other)); return *this; }

    constexpr size_t getCount() const { return s_nodes[this->getGeometry()].cols(); }
    constexpr const Math::PointMatrix& getNodes() const { return s_nodes[this->getGeometry()]; }
    const auto& getLinearForm(size_t i) const { assert(i < getCount()); return s_ls[this->getGeometry()][i]; }
    const auto& getBasis(size_t i) const { assert(i < getCount()); return s_basis[this->getGeometry()][i]; }
    const auto& getGradient(size_t i) const { assert(i < getCount()); return s_gradient[this->getGeometry()][i]; }
    constexpr size_t getOrder() const { return Degree; }
    
  private:
    // Static members are indexed by the reference polytope type.
    static const Geometry::GeometryIndexed<Math::PointMatrix> s_nodes;
    static const Geometry::GeometryIndexed<std::vector<LinearForm>> s_ls;
    static const Geometry::GeometryIndexed<std::vector<BasisFunction>> s_basis;
    static const Geometry::GeometryIndexed<std::vector<GradientFunction>> s_gradient;
  };

  // ─────────────────────────────────────────────────────────────
  //  Static member definitions for PkElement<Real, Degree>
  //
  // For each reference polytope (iterating over Geometry::Polytope::Types),
  // we generate the nodal coordinates, build the Vandermonde system, invert it
  // (using Eigen’s .inverse()), and then build the basis/gradient functions.
  template <size_t Degree>
  const Geometry::GeometryIndexed<Math::PointMatrix> PkElement<Real, Degree>::s_nodes = []() {
    Geometry::GeometryIndexed<Math::PointMatrix> nodes;
    for (auto g : Geometry::Polytope::Types)
    {
      nodes[g] = Internal::generateNodes(g, Degree);
    }
    return nodes;
  }();

  template <size_t Degree>
  const Geometry::GeometryIndexed<std::vector<typename PkElement<Real, Degree>::LinearForm>> 
  PkElement<Real, Degree>::s_ls = []() {
    Geometry::GeometryIndexed<std::vector<LinearForm>> ls;
    for (auto g : Geometry::Polytope::Types)
    {
      Math::PointMatrix N = s_nodes[g];
      size_t numNodes = N.cols();
      std::vector<LinearForm> forms;
      forms.reserve(numNodes);
      for (size_t i = 0; i < numNodes; i++)
        forms.push_back(LinearForm(i, g));
      ls[g] = forms;
    }
    return ls;
  }();

  template <size_t Degree>
  const Geometry::GeometryIndexed<std::vector<typename PkElement<Real, Degree>::BasisFunction>>
  PkElement<Real, Degree>::s_basis = []() {
    Geometry::GeometryIndexed<std::vector<PkElement<Real, Degree>::BasisFunction>> basis;
    for (auto g : Geometry::Polytope::Types)
    {
      Math::PointMatrix N = s_nodes[g];
      size_t numNodes = N.cols();
      size_t dim = Geometry::Polytope::getGeometryDimension(g);
      auto multiIndices = Internal::generateMultiIndices(dim, Degree);
      Math::Matrix<Real> V = Internal::buildVandermonde(N, multiIndices);
      Math::Matrix<Real> invV = Internal::invertMatrix(V);
      std::vector<PkElement<Real, Degree>::BasisFunction> funcs;
      funcs.reserve(numNodes);
      // Each column of invV gives the coefficients for one Lagrange basis function.
      for (size_t i = 0; i < numNodes; i++)
      {
        std::vector<Real> coeffs;
        coeffs.reserve(invV.rows());
        for (size_t j = 0; j < invV.rows(); j++)
          coeffs.push_back(invV(j, i));
        funcs.push_back(BasisFunction(coeffs, multiIndices));
      }
      basis[g] = funcs;
    }
    return basis;
  }();

  template <size_t Degree>
  const Geometry::GeometryIndexed<std::vector<typename PkElement<Real, Degree>::GradientFunction>>
  PkElement<Real, Degree>::s_gradient = []() {
    Geometry::GeometryIndexed<std::vector<PkElement<Real, Degree>::GradientFunction>> grads;
    for (auto g : Geometry::Polytope::Types)
    {
      Math::PointMatrix N = s_nodes[g];
      size_t numNodes = N.cols();
      size_t dim = Geometry::Polytope::getGeometryDimension(g);
      auto multiIndices = Internal::generateMultiIndices(dim, Degree);
      Math::Matrix<Real> V = Internal::buildVandermonde(N, multiIndices);
      Math::Matrix<Real> invV = Internal::invertMatrix(V);
      std::vector<PkElement<Real, Degree>::GradientFunction> gfs;
      gfs.reserve(numNodes);
      for (size_t i = 0; i < numNodes; i++)
      {
        std::vector<Real> coeffs;
        coeffs.reserve(invV.rows());
        for (size_t j = 0; j < invV.rows(); j++)
          coeffs.push_back(invV(j, i));
        gfs.push_back(GradientFunction(coeffs, multiIndices));
      }
      grads[g] = gfs;
    }
    return grads;
  }();

  // ─────────────────────────────────────────────────────────────
  //  PkElement<Complex, Degree>: degree-k scalar Lagrange element (Complex)
  // ─────────────────────────────────────────────────────────────
  template <size_t Degree>
  class PkElement<Complex, Degree> final : public FiniteElementBase<PkElement<Complex, Degree>>
  {
    using G = Geometry::Polytope::Type;
  public:
    using Parent = FiniteElementBase<PkElement<Complex, Degree>>;
    using RangeType = Complex;

    // Linear forms (evaluation at nodes) -- note the extra factor.
    class LinearForm
    {
    public:
      constexpr LinearForm(size_t i, Geometry::Polytope::Type g)
        : m_i(i), m_g(g) {}
      template <class T>
      constexpr auto operator()(const T& v) const
      {
        // Multiply by the complex constant (1,-1) as in our design.
        return 0.5 * std::conj(v(s_nodes[m_g].col(m_i))) * Complex(1, -1);
      }
    private:
      size_t m_i;
      Geometry::Polytope::Type m_g;
    };

    // Basis function stored in the monomial basis.
    class BasisFunction
    {
    public:
      BasisFunction(const std::vector<Real>& coeffs,
                    const std::vector<std::vector<size_t>>& multiIndices)
        : m_coeffs(coeffs), m_multiIndices(multiIndices) {}
      Complex operator()(const Math::SpatialVector<Real>& r) const
      {
        auto mons = Internal::evaluateMonomials(r, m_multiIndices);
        Real sum = 0;
        for (size_t i = 0; i < m_coeffs.size(); i++)
          sum += m_coeffs[i] * mons[i];
        return sum * Complex(1, -1);
      }
    private:
      std::vector<Real> m_coeffs;
      std::vector<std::vector<size_t>> m_multiIndices;
    };

    // Gradient function (complex-valued)
    class GradientFunction
    {
    public:
      GradientFunction(const std::vector<Real>& coeffs,
                       const std::vector<std::vector<size_t>>& multiIndices)
        : m_coeffs(coeffs), m_multiIndices(multiIndices) {}
      Math::SpatialVector<Complex> operator()(const Math::SpatialVector<Real>& r) const
      {
        size_t dim = r.size();
        Math::SpatialVector<Complex> grad = Math::SpatialVector<Complex>::Zero(dim);
        for (size_t i = 0; i < m_coeffs.size(); i++)
        {
          for (size_t d = 0; d < dim; d++)
          {
            if (m_multiIndices[i][d] > 0)
            {
              std::vector<size_t> exp = m_multiIndices[i];
              exp[d]--;
              grad[d] += m_coeffs[i] * m_multiIndices[i][d] *
                         Internal::evaluateMonomial(r, exp);
            }
          }
        }
        // Multiply the entire gradient by the complex constant.
        for (size_t d = 0; d < dim; d++)
          grad[d] *= Complex(1, -1);
        return grad;
      }
    private:
      std::vector<Real> m_coeffs;
      std::vector<std::vector<size_t>> m_multiIndices;
    };

    constexpr PkElement() = default;
    constexpr PkElement(Geometry::Polytope::Type geometry)
      : Parent(geometry) {}
    constexpr PkElement(const PkElement& other) : Parent(other) {}
    constexpr PkElement(PkElement&& other) : Parent(std::move(other)) {}
    constexpr PkElement& operator=(const PkElement& other) { Parent::operator=(other); return *this; }
    constexpr PkElement& operator=(PkElement&& other) { Parent::operator=(std::move(other)); return *this; }

    constexpr size_t getCount() const { return s_nodes[this->getGeometry()].cols(); }
    constexpr const Math::PointMatrix& getNodes() const { return s_nodes[this->getGeometry()]; }
    const auto& getLinearForm(size_t i) const { assert(i < getCount()); return s_ls[this->getGeometry()][i]; }
    const auto& getBasis(size_t i) const { assert(i < getCount()); return s_basis[this->getGeometry()][i]; }
    const auto& getGradient(size_t i) const { assert(i < getCount()); return s_gradient[this->getGeometry()][i]; }
    constexpr size_t getOrder() const { return Degree; }
    
  private:
    static const Geometry::GeometryIndexed<Math::PointMatrix> s_nodes;
    static const Geometry::GeometryIndexed<std::vector<LinearForm>> s_ls;
    static const Geometry::GeometryIndexed<std::vector<BasisFunction>> s_basis;
    static const Geometry::GeometryIndexed<std::vector<GradientFunction>> s_gradient;
  };

  template <size_t Degree>
  const Geometry::GeometryIndexed<Math::PointMatrix> PkElement<Complex, Degree>::s_nodes = PkElement<Real, Degree>::s_nodes;

  template <size_t Degree>
  const Geometry::GeometryIndexed<std::vector<typename PkElement<Complex, Degree>::LinearForm>>
  PkElement<Complex, Degree>::s_ls = []() {
    Geometry::GeometryIndexed<std::vector<LinearForm>> ls;
    for (auto g : Geometry::Polytope::Types)
    {
      Math::PointMatrix N = s_nodes[g];
      size_t numNodes = N.cols();
      std::vector<LinearForm> forms;
      forms.reserve(numNodes);
      for (size_t i = 0; i < numNodes; i++)
        forms.push_back(LinearForm(i, g));
      ls[g] = forms;
    }
    return ls;
  }();

  template <size_t Degree>
  const Geometry::GeometryIndexed<std::vector<typename PkElement<Complex, Degree>::BasisFunction>>
  PkElement<Complex, Degree>::s_basis = []() {
    Geometry::GeometryIndexed<std::vector<typename PkElement<Complex, Degree>::BasisFunction>> basis;
    for (auto g : Geometry::Polytope::Types)
    {
      Math::PointMatrix N = s_nodes[g];
      size_t numNodes = N.cols();
      size_t dim = Geometry::Polytope::getGeometryDimension(g);
      auto multiIndices = Internal::generateMultiIndices(dim, Degree);
      Math::Matrix<Real> V = Internal::buildVandermonde(N, multiIndices);
      Math::Matrix<Real> invV = Internal::invertMatrix(V);
      std::vector<typename PkElement<Complex, Degree>::BasisFunction> funcs;
      funcs.reserve(numNodes);
      for (size_t i = 0; i < numNodes; i++)
      {
        std::vector<Real> coeffs;
        coeffs.reserve(invV.rows());
        for (size_t j = 0; j < invV.rows(); j++)
          coeffs.push_back(invV(j, i));
        funcs.push_back(BasisFunction(coeffs, multiIndices));
      }
      basis[g] = funcs;
    }
    return basis;
  }();

  template <size_t Degree>
  const Geometry::GeometryIndexed<std::vector<typename PkElement<Complex, Degree>::GradientFunction>>
  PkElement<Complex, Degree>::s_gradient = []() {
    Geometry::GeometryIndexed<std::vector<typename PkElement<Complex, Degree>::GradientFunction>> grads;
    for (auto g : Geometry::Polytope::Types)
    {
      Math::PointMatrix N = s_nodes[g];
      size_t numNodes = N.cols();
      size_t dim = Geometry::Polytope::getGeometryDimension(g);
      auto multiIndices = Internal::generateMultiIndices(dim, Degree);
      Math::Matrix<Real> V = Internal::buildVandermonde(N, multiIndices);
      Math::Matrix<Real> invV = Internal::invertMatrix(V);
      std::vector<typename PkElement<Complex, Degree>::GradientFunction> gfs;
      gfs.reserve(numNodes);
      for (size_t i = 0; i < numNodes; i++)
      {
        std::vector<Real> coeffs;
        coeffs.reserve(invV.rows());
        for (size_t j = 0; j < invV.rows(); j++)
          coeffs.push_back(invV(j, i));
        gfs.push_back(GradientFunction(coeffs, multiIndices));
      }
      grads[g] = gfs;
    }
    return grads;
  }();

  // ─────────────────────────────────────────────────────────────
  //  PkElement<Math::Vector<Real>, Degree>: degree-k vector Lagrange element
  // ─────────────────────────────────────────────────────────────


  template <size_t Degree>
  class PkElement<Math::Vector<Real>, Degree> final : public FiniteElementBase<PkElement<Math::Vector<Real>, Degree>>
  {
    using G = Geometry::Polytope::Type;
  public:
    using Parent = FiniteElementBase<PkElement<Math::Vector<Real>, Degree>>;
    using RangeType = Math::Vector<Real>;

    // Linear form: selects the appropriate component.
    class LinearForm
    {
    public:
      constexpr LinearForm() : m_vdim(0), m_i(0), m_g(Geometry::Polytope::Type::Point) {}
      constexpr LinearForm(size_t vdim, size_t i, Geometry::Polytope::Type g)
        : m_vdim(vdim), m_i(i), m_g(g)
      {
        assert(m_vdim > 0);
      }
      template <class T>
      constexpr auto operator()(const T& v) const
      {
        return v(s_nodes[m_vdim][m_g].col(m_i)).coeff(static_cast<size_t>(m_i % m_vdim));
      }
    private:
      size_t m_vdim;
      size_t m_i;
      Geometry::Polytope::Type m_g;
    };

    // Basis function: its scalar part is as in the scalar case, placed in the appropriate component.
    class BasisFunction
    {
    public:
      constexpr BasisFunction() : m_vdim(0), m_i(0), m_g(Geometry::Polytope::Type::Point) {}
      constexpr BasisFunction(size_t vdim, size_t i, Geometry::Polytope::Type g)
        : m_vdim(vdim), m_i(i), m_g(g)
      {
        assert(m_vdim > 0);
      }
      // Evaluate and write into out.
      void operator()(Math::Vector<Real>& out, const Math::SpatialVector<Real>& r) const
      {
        out = Math::Vector<Real>::Zero(m_vdim);
        size_t comp = m_i % m_vdim;
        size_t node = m_i / m_vdim;
        // Use the scalar PkElement to get the value at the node.
        auto val = PkElement<Real, Degree>(m_g).getBasis(node)(r);
        out[comp] = val;
      }
    private:
      size_t m_vdim;
      size_t m_i;
      Geometry::Polytope::Type m_g;
    };

    // Jacobian function: derivative with respect to spatial coordinates.
    class JacobianFunction
    {
    public:
      constexpr JacobianFunction() : m_vdim(0), m_i(0), m_g(Geometry::Polytope::Type::Point) {}
      constexpr JacobianFunction(size_t vdim, size_t i, Geometry::Polytope::Type g)
        : m_vdim(vdim), m_i(i), m_g(g)
      {
        assert(m_vdim > 0);
      }
      void operator()(Math::SpatialMatrix<Real>& out, const Math::SpatialVector<Real>& r) const
      {
        out = Math::SpatialMatrix<Real>::Zero(m_vdim, Geometry::Polytope::getGeometryDimension(m_g));
        size_t comp = m_i % m_vdim;
        size_t node = m_i / m_vdim;
        // Use the scalar PkElement gradient for the node.
        auto grad = PkElement<Real, Degree>(m_g).getGradient(node)(r);
        out.row(comp) = grad.transpose();
      }
    private:
      size_t m_vdim;
      size_t m_i;
      Geometry::Polytope::Type m_g;
    };

    PkElement() = default;
    constexpr PkElement(size_t vdim, Geometry::Polytope::Type geometry)
      : Parent(geometry), m_vdim(vdim) {}
    constexpr PkElement(const PkElement& other)
      : Parent(other), m_vdim(other.m_vdim) {}
    constexpr PkElement(PkElement&& other)
      : Parent(std::move(other)), m_vdim(other.m_vdim) {}
    PkElement& operator=(const PkElement& other)
    {
      if (this != &other)
      {
        Parent::operator=(other);
        m_vdim = other.m_vdim;
      }
      return *this;
    }
    PkElement& operator=(PkElement&& other)
    {
      if (this != &other)
      {
        Parent::operator=(std::move(other));
        m_vdim = other.m_vdim;
      }
      return *this;
    }
    constexpr size_t getCount() const
    {
      return s_nodes[m_vdim][this->getGeometry()].cols();
    }
    constexpr const Math::PointMatrix& getNodes() const
    {
      return s_nodes[m_vdim][this->getGeometry()];
    }
    const auto& getLinearForm(size_t i) const
    {
      assert(i < getCount());
      return s_ls[m_vdim][this->getGeometry()][i];
    }
    const auto& getBasis(size_t i) const
    {
      assert(i < getCount());
      return s_basis[m_vdim][this->getGeometry()][i];
    }
    const auto& getJacobian(size_t i) const
    {
      assert(i < getCount());
      return s_jacobian[m_vdim][this->getGeometry()][i];
    }
    constexpr size_t getOrder() const { return Degree; }
    
  private:
    // For vector elements we store an array (indexed by the vector dimension)
    // of static nodal/basis/jacobian data.
    static const std::array<Geometry::GeometryIndexed<Math::PointMatrix>, RODIN_PK_MAX_VECTOR_DIMENSION> s_nodes;
    static const std::array<Geometry::GeometryIndexed<std::vector<LinearForm>>, RODIN_PK_MAX_VECTOR_DIMENSION> s_ls;
    static const std::array<Geometry::GeometryIndexed<std::vector<BasisFunction>>, RODIN_PK_MAX_VECTOR_DIMENSION> s_basis;
    static const std::array<Geometry::GeometryIndexed<std::vector<JacobianFunction>>, RODIN_PK_MAX_VECTOR_DIMENSION> s_jacobian;
    
    size_t m_vdim;
  };

  // ─────────────────────────────────────────────────────────────
  //  Static member definitions for PkElement<Math::Vector<Real>, Degree>
  // (These are built by “lifting” the scalar data via replication over components.)
  namespace Internal
  {
    template <size_t Degree>
    std::array<Geometry::GeometryIndexed<Math::PointMatrix>, RODIN_PK_MAX_VECTOR_DIMENSION>
    initVectorPkNodes()
    {
      std::array<Geometry::GeometryIndexed<Math::PointMatrix>, RODIN_PK_MAX_VECTOR_DIMENSION> res;
      // vdim = 0 is unused.
      for (size_t vdim = 1; vdim < RODIN_PK_MAX_VECTOR_DIMENSION; vdim++)
      {
        for (auto g : Geometry::Polytope::Types)
        {
          if (g == Geometry::Polytope::Type::Point)
          {
            res[vdim][g] = Math::PointMatrix::Zero(1, vdim);
          }
          else
          {
            auto scalarEl = PkElement<Real, Degree>(g);
            const Math::PointMatrix& scalarNodes = scalarEl.getNodes();
            const size_t d = Geometry::Polytope::getGeometryDimension(g);
            res[vdim][g].resize(d, scalarNodes.cols() * vdim);
            for (size_t i = 0; i < scalarNodes.cols(); i++)
            {
              for (size_t c = 0; c < vdim; c++)
                res[vdim][g].col(i * vdim + c) = scalarNodes.col(i);
            }
          }
        }
      }
      return res;
    }

    template <size_t Degree>
    std::array<Geometry::GeometryIndexed<std::vector<typename PkElement<Math::Vector<Real>, Degree>::LinearForm>>, RODIN_PK_MAX_VECTOR_DIMENSION>
    initVectorPkLinearForms()
    {
      std::array<Geometry::GeometryIndexed<std::vector<typename PkElement<Math::Vector<Real>, Degree>::LinearForm>>, RODIN_PK_MAX_VECTOR_DIMENSION> res;
      for (size_t vdim = 1; vdim < RODIN_PK_MAX_VECTOR_DIMENSION; vdim++)
      {
        for (auto g : Geometry::Polytope::Types)
        {
          const size_t numNodes = Internal::generateNodes(g, Degree).cols();
          res[vdim][g].resize(numNodes * vdim);
          for (size_t i = 0; i < numNodes * vdim; i++)
            res[vdim][g][i] = typename PkElement<Math::Vector<Real>, Degree>::LinearForm(vdim, i, g);
        }
      }
      return res;
    }

    template <size_t Degree>
    std::array<Geometry::GeometryIndexed<std::vector<typename PkElement<Math::Vector<Real>, Degree>::BasisFunction>>, RODIN_PK_MAX_VECTOR_DIMENSION>
    initVectorPkBasis()
    {
      std::array<Geometry::GeometryIndexed<std::vector<typename PkElement<Math::Vector<Real>, Degree>::BasisFunction>>, RODIN_PK_MAX_VECTOR_DIMENSION> res;
      for (size_t vdim = 1; vdim < RODIN_PK_MAX_VECTOR_DIMENSION; vdim++)
      {
        for (auto g : Geometry::Polytope::Types)
        {
          const size_t numNodes = Internal::generateNodes(g, Degree).cols();
          res[vdim][g].resize(numNodes * vdim);
          for (size_t i = 0; i < numNodes * vdim; i++)
            res[vdim][g][i] = typename PkElement<Math::Vector<Real>, Degree>::BasisFunction(vdim, i, g);
        }
      }
      return res;
    }

    template <size_t Degree>
    std::array<Geometry::GeometryIndexed<std::vector<typename PkElement<Math::Vector<Real>, Degree>::JacobianFunction>>, RODIN_PK_MAX_VECTOR_DIMENSION>
    initVectorPkJacobian()
    {
      std::array<Geometry::GeometryIndexed<std::vector<typename PkElement<Math::Vector<Real>, Degree>::JacobianFunction>>, RODIN_PK_MAX_VECTOR_DIMENSION> res;
      for (size_t vdim = 1; vdim < RODIN_PK_MAX_VECTOR_DIMENSION; vdim++)
      {
        for (auto g : Geometry::Polytope::Types)
        {
          const size_t numNodes = Internal::generateNodes(g, Degree).cols();
          res[vdim][g].resize(numNodes * vdim);
          for (size_t i = 0; i < numNodes * vdim; i++)
            res[vdim][g][i] = typename PkElement<Math::Vector<Real>, Degree>::JacobianFunction(vdim, i, g);
        }
      }
      return res;
    }
  } // namespace Internal

  template <size_t Degree>
  const std::array<Geometry::GeometryIndexed<Math::PointMatrix>, RODIN_PK_MAX_VECTOR_DIMENSION>
  PkElement<Math::Vector<Real>, Degree>::s_nodes = Internal::initVectorPkNodes<Degree>();

  template <size_t Degree>
  const std::array<Geometry::GeometryIndexed<std::vector<typename PkElement<Math::Vector<Real>, Degree>::LinearForm>>, RODIN_PK_MAX_VECTOR_DIMENSION>
  PkElement<Math::Vector<Real>, Degree>::s_ls = Internal::initVectorPkLinearForms<Degree>();

  template <size_t Degree>
  const std::array<Geometry::GeometryIndexed<std::vector<typename PkElement<Math::Vector<Real>, Degree>::BasisFunction>>, RODIN_PK_MAX_VECTOR_DIMENSION>
  PkElement<Math::Vector<Real>, Degree>::s_basis = Internal::initVectorPkBasis<Degree>();

  template <size_t Degree>
  const std::array<Geometry::GeometryIndexed<std::vector<typename PkElement<Math::Vector<Real>, Degree>::JacobianFunction>>, RODIN_PK_MAX_VECTOR_DIMENSION>
  PkElement<Math::Vector<Real>, Degree>::s_jacobian = Internal::initVectorPkJacobian<Degree>();

} // namespace Rodin::Variational

#endif

