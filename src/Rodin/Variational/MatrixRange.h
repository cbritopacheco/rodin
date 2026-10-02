/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Matrix finite elements and component-wise spatial evaluation.
 */
#ifndef RODIN_VARIATIONAL_MATRIXRANGE_H
#define RODIN_VARIATIONAL_MATRIXRANGE_H

#include <map>
#include <type_traits>
#include <vector>

#include "Rodin/Math/SpatialMatrix.h"
#include "Rodin/Math/SpatialTensor.h"
#include "FiniteElement.h"
#include "FiniteElementSpace.h"
#include "IntegrationPoint.h"
#include "ShapeFunction.h"

namespace Rodin::Variational::Detail
{
  template <class ScalarElement>
  class MatrixElement;
}

namespace Rodin::FormLanguage
{
  /// @brief Type traits for the matrix or tensor expression specialization.
  template <class ScalarElement>
  struct Traits<Variational::Detail::MatrixElement<ScalarElement>>
  {
      /// @brief Scalar type of matrix or tensor entries.
      using ScalarType = typename Traits<ScalarElement>::ScalarType;
      /// @brief Evaluated matrix, tensor, or scalar range type.
      using RangeType = Math::SpatialMatrix<ScalarType>;
  };
}

namespace Rodin::Variational::Detail
{
  /// @brief Rejects zero or oversized matrix axis extents.
  inline void checkMatrixDimensions(size_t rows, size_t cols)
  {
    if (rows == 0 || cols == 0 || rows > RODIN_MAXIMAL_SPACE_DIMENSION ||
      cols > RODIN_MAXIMAL_SPACE_DIMENSION)
      Alert::Exception() << "SpatialMatrix ranges require 1 to 3 rows and columns."
                         << Alert::Raise;
  }

  /** Componentwise nodal element. Components are ordered by row, then column. */
  template <class ScalarElement>
  class MatrixElement : public FiniteElementBase<MatrixElement<ScalarElement>>
  {
    public:
      /// @brief Scalar type of matrix or tensor entries.
      using ScalarType = typename FormLanguage::Traits<ScalarElement>::ScalarType;
      /// @brief Evaluated matrix, tensor, or scalar range type.
      using RangeType = Math::SpatialMatrix<ScalarType>;
      /// @brief CRTP or finite element base class.
      using Parent = FiniteElementBase<MatrixElement>;
      /// @brief Basis type of the underlying scalar element.
      using ScalarBasis = decltype(std::declval<ScalarElement>().getBasis(0));
      /// @brief Nodal functional type of the underlying scalar element.
      using ScalarLinearForm = decltype(std::declval<ScalarElement>().getLinearForm(0));

      /// @brief Replicates a scalar reference element for a rectangular matrix range.
      MatrixElement(Geometry::Polytope::Type geometry, size_t rows, size_t cols)
        : MatrixElement(ScalarElement(geometry), rows, cols)
      {}

      /// @brief Replicates a scalar reference element for a rectangular matrix range.
      MatrixElement(const ScalarElement& scalar, size_t rows, size_t cols)
        : Parent(scalar.getGeometry()),
          m_scalar(scalar),
          m_rows(rows),
          m_cols(cols)
      {
        checkMatrixDimensions(rows, cols);
      }

      /// @brief Matrix-range finite element or expression specialization.
      class BasisFunction
      {
        public:
          /// @brief Selects a matrix unit multiplied by a scalar basis function.
          BasisFunction(ScalarBasis basis, size_t rows, size_t cols, size_t component)
            : m_basis(std::move(basis)),
              m_rows(rows),
              m_cols(cols),
              m_component(component)
          {}

          /// @brief Evaluates the selected matrix basis or its component nodal functional.
          RangeType operator()(const Math::SpatialPoint& point) const
          {
            RangeType value(m_rows, m_cols);
            value.setZero();
            value(m_component / m_cols, m_component % m_cols) = m_basis(point);
            return value;
          }

          /// @brief Returns a reference-coordinate derivative of the basis.
          template <size_t Order>
          auto getDerivative(size_t direction) const
          {
            return [derivative = m_basis.template getDerivative<Order>(direction),
                     rows = m_rows, cols = m_cols, component = m_component](
                     const Math::SpatialPoint& point) -> RangeType {
              RangeType value(rows, cols);
              value.setZero();
              value(component / cols, component % cols) = derivative(point);
              return value;
            };
          }

        private:
          ScalarBasis m_basis;
          size_t m_rows, m_cols, m_component;
      };

      /// @brief Matrix-range finite element or expression specialization.
      class LinearForm
      {
        public:
          /// @brief Selects a matrix entry for the scalar nodal functional.
          LinearForm(ScalarLinearForm form, size_t rows, size_t cols, size_t component)
            : m_form(std::move(form)),
              m_rows(rows),
              m_cols(cols),
              m_component(component)
          {}

          /// @brief Evaluates the selected matrix basis or its component nodal functional.
          template <class Callable>
          ScalarType operator()(const Callable& function) const
          {
            return m_form([&](const Math::SpatialPoint& point) -> ScalarType {
              const auto value = function(point);
              if (value.rows() != m_rows || value.cols() != m_cols)
                Alert::Exception()
                  << "Matrix value does not match the finite element range."
                  << Alert::Raise;
              return value(m_component / m_cols, m_component % m_cols);
            });
          }

        private:
          ScalarLinearForm m_form;
          size_t m_rows, m_cols, m_component;
      };

      /// @brief Contracts each matrix entry with the scalar reference element.
      template <class Coefficient>
      void evaluate(
        RangeType& out, Coefficient&& coefficient, const Math::SpatialPoint& point) const
      {
        out.resize(m_rows, m_cols);
        const size_t components = m_rows * m_cols;
        for (size_t c = 0; c < components; ++c)
          m_scalar.evaluate(
            out(c / m_cols, c % m_cols),
            [&](size_t a) { return coefficient(a * components + c); }, point);
      }

      /// @brief Returns the number of local matrix basis functions.
      size_t getCount() const
      {
        return m_scalar.getCount() * m_rows * m_cols;
      }
      /// @brief Returns the polynomial order when it is known.
      size_t getOrder() const
      {
        return m_scalar.getOrder();
      }
      /// @brief Returns the number of matrix rows.
      size_t getRows() const
      {
        return m_rows;
      }
      /// @brief Returns the number of matrix columns.
      size_t getColumns() const
      {
        return m_cols;
      }
      /// @brief Returns the scalar element whose basis is replicated for matrix entries.
      const ScalarElement& getScalarElement() const
      {
        return m_scalar;
      }

      /// @brief Returns the scalar reference node for the selected component DOF.
      decltype(auto) getNode(size_t local) const
      {
        assert(local < getCount());
        return m_scalar.getNode(local / (m_rows * m_cols));
      }

      /// @brief Returns a basis value at the bound integration point.
      BasisFunction getBasis(size_t local) const
      {
        assert(local < getCount());
        return {m_scalar.getBasis(local / (m_rows * m_cols)), m_rows, m_cols,
          local % (m_rows * m_cols)};
      }

      /// @brief Returns the scalar nodal functional applied to the selected matrix entry.
      LinearForm getLinearForm(size_t local) const
      {
        assert(local < getCount());
        return {m_scalar.getLinearForm(local / (m_rows * m_cols)), m_rows, m_cols,
          local % (m_rows * m_cols)};
      }

    private:
      ScalarElement m_scalar;
      size_t m_rows, m_cols;
  };

  /**
   * @brief Replicates scalar DOF maps for a rectangular matrix range.
   * @par Architecture
   * The scalar family supplies reference bases, entity ordering and pullbacks.
   * Construction expands each scalar DOF into row-major matrix components and
   * prebuilds the reference elements for the mesh geometries. Global size is
   * resolved once; a distributed scalar size must not be queried in cell loops.
   * Differential evaluation uses the physical inverse Jacobian with the
   * derivative axis last. MPI ownership is supplied by DistributedMatrixSpace.
   */
  template <class Derived, class ScalarSpace, class Element, bool CellsOnly = false>
  class MatrixSpace : public FiniteElementSpace<typename ScalarSpace::MeshType, Derived>
  {
    public:
      /// @brief Scalar type of matrix or tensor entries.
      using ScalarType = typename ScalarSpace::ScalarType;
      /// @brief Evaluated matrix, tensor, or scalar range type.
      using RangeType = Math::SpatialMatrix<ScalarType>;
      /// @brief Mesh type supplying topology and physical transformations.
      using MeshType = typename ScalarSpace::MeshType;
      /// @brief Local or distributed execution context.
      using ContextType = typename ScalarSpace::ContextType;
      /// @brief Reference finite element type.
      using ElementType = Element;

      /// @brief Expands scalar DOF maps into interleaved row-major matrix components.
      MatrixSpace(ScalarSpace scalar, size_t rows, size_t cols)
        : m_scalar(std::move(scalar)),
          m_rows(rows),
          m_cols(cols)
      {
        checkMatrixDimensions(rows, cols);
        // Distributed scalar spaces may compute their size collectively.
        // Resolve it at construction, never from a local evaluation loop.
        m_size = m_scalar.getSize() * rows * cols;
        const auto& mesh = getMesh();
        m_dofs.resize(mesh.getDimension() + 1);
        for (size_t d = CellsOnly ? mesh.getDimension() : 0; d <= mesh.getDimension();
          ++d)
        {
          const size_t count = mesh.getConnectivity().getCount(d);
          m_dofs[d].reserve(count);
          for (size_t i = 0; i < count; ++i)
          {
            const auto& scalarDOFs = m_scalar.getDOFs(d, i);
            auto& dofs = m_dofs[d].emplace_back(scalarDOFs.size() * rows * cols);
            for (size_t a = 0; a < static_cast<size_t>(scalarDOFs.size()); ++a)
              for (size_t c = 0; c < rows * cols; ++c)
                dofs[a * rows * cols + c] = scalarDOFs[a] * rows * cols + c;
            const auto& scalarFE = m_scalar.getFiniteElement(d, i);
            m_elements.try_emplace(scalarFE.getGeometry(), scalarFE, rows, cols);
          }
        }
      }

      size_t getSize() const override
      {
        return m_size;
      }
      size_t getVectorDimension() const override
      {
        return m_rows * m_cols;
      }
      /// @brief Returns the number of matrix rows.
      size_t getRows() const
      {
        return m_rows;
      }
      /// @brief Returns the number of matrix columns.
      size_t getColumns() const
      {
        return m_cols;
      }
      const MeshType& getMesh() const override
      {
        return m_scalar.getMesh();
      }
      /// @brief Returns the scalar space supplying topology and component-independent maps.
      const ScalarSpace& getScalarSpace() const
      {
        return m_scalar;
      }

      /// @brief Returns the matrix reference element of a mesh entity.
      const ElementType& getFiniteElement(size_t d, Index i) const
      {
        return m_elements.at(m_scalar.getFiniteElement(d, i).getGeometry());
      }

      const IndexArray& getDOFs(size_t d, Index i) const override
      {
        return m_dofs.at(d).at(i);
      }

      Index getGlobalIndex(const std::pair<size_t, Index>& p, Index local) const override
      {
        const Index components = m_rows * m_cols;
        return m_scalar.getGlobalIndex(p, local / components) * components +
          local % components;
      }

      /**
       * @brief Selects the volume-side point used by matrix derivatives.
       * Interior faces require an explicit trace domain; boundary faces have
       * one adjacent cell. Incidences must be prepared by the caller.
       */
      Geometry::Point getDerivativePoint(const Geometry::Point& point,
        const FlatSet<Geometry::Attribute>& traceDomain) const
      {
        const auto& mesh = getMesh();
        if (!mesh.isLocalPoint(point))
        {
          if (auto included = mesh.inclusion(point))
            return getDerivativePoint(*included, traceDomain);
          if (mesh.isSubMesh())
            if (auto restricted = mesh.asSubMesh().restriction(point))
              return getDerivativePoint(*restricted, traceDomain);
          Alert::Exception() << "Derivative point is outside the finite element mesh."
                             << Alert::Raise;
        }
        const auto& poly = point.getPolytope();
        const size_t dimension = mesh.getDimension();
        if (poly.getDimension() == dimension)
          return point;
        const auto& adjacent = mesh.getConnectivity().getIncidence(
          {poly.getDimension(), dimension}, poly.getIndex());
        for (auto index : adjacent)
        {
          const auto cell = mesh.getPolytope(dimension, index);
          const auto attribute = cell->getAttribute();
          if (adjacent.size() != 1 && (!attribute || !traceDomain.contains(*attribute)))
            continue;
          Math::SpatialPoint reference;
          cell->getTransformation().inverse(reference, point.getPhysicalCoordinates());
          return Geometry::Point(*cell, reference);
        }
        Alert::Exception() << "Matrix derivative requires a determined trace domain."
                           << Alert::Raise;
        return point;
      }

      /**
       * @brief Physical gradient of a matrix basis, with the derivative axis last.
       * @f$ G_{ijk}=\sum_l\partial_l\widehat\phi_a(J^{-1})_{lk}(E_c)_{ij} @f$.
       */
      Math::SpatialTensor<ScalarType> getGradientBasis(
        size_t local, const Geometry::Point& point) const
      {
        const auto& poly = point.getPolytope();
        const auto& fe = getFiniteElement(poly.getDimension(), poly.getIndex());
        const size_t components = m_rows * m_cols;
        const auto& basis = fe.getScalarElement().getBasis(local / components);
        const auto& inverse = point.getJacobianInverse();
        Math::SpatialTensor<ScalarType> gradient(
          m_rows, m_cols, getMesh().getSpaceDimension());
        gradient.setZero();
        const size_t row = (local % components) / m_cols;
        const size_t column = local % m_cols;
        for (size_t l = 0; l < poly.getDimension(); ++l)
        {
          const auto derivative =
            basis.template getDerivative<1>(l)(point.getReferenceCoordinates());
          for (size_t k = 0; k < getMesh().getSpaceDimension(); ++k)
            gradient(row, column, k) += derivative * inverse(l, k);
        }
        return gradient;
      }

      /// @brief Pulls a physical callable back to a reference element.
      template <class Callable>
      auto getPullback(const std::pair<size_t, Index>& p, Callable&& value) const
      {
        return m_scalar.getPullback(p, std::forward<Callable>(value));
      }

      /// @brief Pushes a reference callable forward to the physical mesh.
      template <class Callable>
      auto getPushforward(const std::pair<size_t, Index>& p, Callable&& value) const
      {
        return m_scalar.getPushforward(p, std::forward<Callable>(value));
      }

    private:
      ScalarSpace m_scalar;
      size_t m_rows, m_cols;
      size_t m_size = 0;
      std::vector<std::vector<IndexArray>> m_dofs;
      std::map<Geometry::Polytope::Type, ElementType> m_elements;
  };

  /// @brief Matrix-range finite element or expression specialization.
  template <class Shape, class Derived, class FES, ShapeFunctionSpaceType Space>
  class MatrixShape : public ShapeFunctionBase<Shape, FES, Space>
  {
    public:
      /// @brief CRTP or finite element base class.
      using Parent = ShapeFunctionBase<Shape, FES, Space>;
      /// @brief Evaluated matrix, tensor, or scalar range type.
      using RangeType = typename FormLanguage::Traits<FES>::RangeType;
      /// @brief Constructs matrix basis tabulation on the supplied finite element space.
      explicit MatrixShape(const FES& fes)
        : Parent(fes)
      {}
      /// @brief Constructs matrix basis tabulation on the supplied finite element space.
      MatrixShape(const MatrixShape& other)
        : Parent(other)
      {}
      /// @brief Constructs matrix basis tabulation on the supplied finite element space.
      MatrixShape(MatrixShape&& other)
        : Parent(std::move(other))
      {}

      /// @brief Returns the local basis count for the selected polytope.
      size_t getDOFs(const Geometry::Polytope& poly) const
      {
        return this->getFiniteElementSpace()
          .getFiniteElement(poly.getDimension(), poly.getIndex())
          .getCount();
      }

      /// @brief Returns the currently bound integration point.
      const IntegrationPoint& getIntegrationPoint() const
      {
        assert(m_ip);
        return *m_ip;
      }

      /// @brief Binds the integration point and prepares local basis values.
      Shape& setIntegrationPoint(const IntegrationPoint& ip)
      {
        m_ip = &ip;
        const auto& poly = ip.getPoint().getPolytope();
        const auto& fe = this->getFiniteElementSpace().getFiniteElement(
          poly.getDimension(), poly.getIndex());
        if (!ip.getQuadratureFormula() || m_qf != ip.getQuadratureFormula() ||
          m_qp != ip.getIndex() || m_geometry != poly.getGeometry())
        {
          m_basis.resize(fe.getCount());
          const auto& scalar = fe.getScalarElement();
          const size_t components = fe.getRows() * fe.getColumns();
          const auto& point = ip.getPoint().getReferenceCoordinates();
          for (size_t a = 0; a < scalar.getCount(); ++a)
          {
            const auto value = scalar.getBasis(a)(point);
            for (size_t c = 0; c < components; ++c)
            {
              auto& basis = m_basis[a * components + c];
              basis.resize(fe.getRows(), fe.getColumns());
              basis.setZero();
              basis(c / fe.getColumns(), c % fe.getColumns()) = value;
            }
          }
          m_qf = ip.getQuadratureFormula();
          m_qp = ip.getIndex();
          m_geometry = poly.getGeometry();
        }
        return static_cast<Shape&>(*this);
      }

      /// @brief Returns a basis value at the bound integration point.
      const RangeType& getBasis(size_t local) const
      {
        return m_basis.at(local);
      }
      /// @brief Returns the leaf shape function used for assembly.
      const auto& getLeaf() const
      {
        return static_cast<const Derived&>(*this).getLeaf();
      }
      /// @brief Returns the polynomial order when it is known.
      Optional<size_t> getOrder(const Geometry::Polytope& poly) const noexcept
      {
        return this->getFiniteElementSpace()
          .getFiniteElement(poly.getDimension(), poly.getIndex())
          .getOrder();
      }

    private:
      const IntegrationPoint* m_ip = nullptr;
      const QF::QuadratureFormulaBase* m_qf = nullptr;
      size_t m_qp = 0;
      Geometry::Polytope::Type m_geometry = Geometry::Polytope::Type::Point;
      std::vector<RangeType> m_basis;
  };
}

#endif
