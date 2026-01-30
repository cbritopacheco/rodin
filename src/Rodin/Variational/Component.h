/**
 * @file
 * @brief Component extraction from vector and matrix functions.
 */

#ifndef RODIN_VARIATIONAL_COMPONENT_H
#define RODIN_VARIATIONAL_COMPONENT_H

#include "Rodin/Geometry/Point.h"

#include "Rodin/FormLanguage/Traits.h"

#include "ForwardDecls.h"
#include "Rodin/Variational/IntegrationPoint.h"

namespace Rodin::FormLanguage
{
  template <class OperandDerived, class FES, Variational::ShapeFunctionSpaceType Space>
  struct Traits<
    Variational::Component<Variational::ShapeFunctionBase<OperandDerived, FES, Space>>>
  {
    using FESType = FES;
    static constexpr const Variational::ShapeFunctionSpaceType SpaceType = Space;

    using OperandType = Variational::ShapeFunctionBase<OperandDerived, FESType, SpaceType>;

    using OperandRangeType = typename FormLanguage::Traits<OperandType>::RangeType;

    using ScalarType = typename FormLanguage::Traits<OperandRangeType>::ScalarType;

    using RangeType = ScalarType;
  };
}

namespace Rodin::Variational
{
  /**
   * @brief Extracts a single component from a vector-valued function.
   *
   * Given a vector function @f$ \mathbf{f}: \Omega \to \mathbb{R}^d @f$,
   * extracts component @f$ i @f$ such that:
   * @f[
   *    \text{Component}(\mathbf{f}, i)(x) = f_i(x)
   * @f]
   *
   * @tparam OperandDerived Type of the vector function
   *
   * @see Component<FunctionBase<OperandDerived>, size_t, size_t> for matrix components
   */
  template <class OperandDerived>
  class Component<FunctionBase<OperandDerived>, size_t> final
    : public RealFunctionBase<Component<FunctionBase<OperandDerived>, size_t>>
  {
    public:
      using OperandType = FunctionBase<OperandDerived>;

      using OperandRangeType = typename FormLanguage::Traits<OperandType>::RangeType;

      using ScalarType = typename FormLanguage::Traits<OperandRangeType>::ScalarType;

      using RangeType = ScalarType;

      using Parent = RealFunctionBase<Component<FunctionBase<OperandDerived>, size_t>>;


      /**
       * @brief Constructs component extractor for vector function.
       * @param fn Vector function to extract component from
       * @param component Zero-based component index @f$ i @f$
       */
      constexpr
      Component(const OperandType& fn, size_t component)
        : m_fn(fn.copy()), m_idx(component)
      {}

      constexpr
      Component(const Component& other)
        : Parent(other),
          m_fn(other.m_fn->copy()),
          m_idx(other.m_idx)
      {}

      constexpr
      Component(Component&& other)
        : Parent(std::move(other)),
          m_fn(std::move(other.m_fn)),
          m_idx(std::move(other.m_idx))
      {}

      /**
       * @brief Gets the component index.
       * @returns Component index @f$ i @f$
       */
      constexpr
      size_t getIndex() const
      {
        return m_idx;
      }

      /**
       * @brief Gets the underlying vector function.
       * @returns Reference to the operand function
       */
      constexpr
      const OperandType& getOperand() const
      {
        assert(m_fn);
        return *m_fn;
      }

      /**
       * @brief Evaluates the component at a point.
       * @param p Point at which to evaluate
       * @returns Scalar value @f$ f_i(p) @f$
       */
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return this->getOperand().getValue(p).coeff(m_idx);
      }

      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope& geom) const noexcept
      {
        return getOperand().getOrder(geom);
      }

      Component* copy() const noexcept override
      {
        return new Component(*this);
      }

    private:
      std::unique_ptr<OperandType> m_fn;
      const size_t m_idx;
  };

  /**
   * @brief Deduction guide for vector component extraction.
   */
  template <class OperandDerived>
  Component(const FunctionBase<OperandDerived>&, size_t) -> Component<FunctionBase<OperandDerived>, size_t>;

  /**
   * @brief Extracts a single entry from a matrix-valued function.
   *
   * Given a matrix function @f$ A: \Omega \to \mathbb{R}^{m \times n} @f$,
   * extracts entry @f$ (i,j) @f$ such that:
   * @f[
   *    \text{Component}(A, i, j)(x) = A_{ij}(x)
   * @f]
   *
   * @tparam OperandDerived Type of the matrix function
   */
  template <class OperandDerived>
  class Component<FunctionBase<OperandDerived>, size_t, size_t> final
    : public RealFunctionBase<Component<FunctionBase<OperandDerived>, size_t, size_t>>
  {
    public:
      using OperandType = FunctionBase<OperandDerived>;

      using OperandRangeType = typename FormLanguage::Traits<OperandType>::RangeType;

      using ScalarType = typename FormLanguage::Traits<OperandRangeType>::ScalarType;

      using RangeType = ScalarType;

      using Parent = RealFunctionBase<Component<FunctionBase<OperandDerived>, size_t, size_t>>;

      /**
       * @brief Constructs component extractor for matrix function.
       * @param fn Matrix function to extract entry from
       * @param i Zero-based row index
       * @param j Zero-based column index
       */
      constexpr
      Component(const OperandType& fn, size_t i, size_t j)
        : m_fn(fn.copy()), m_i(i), m_j(j)
      {}

      constexpr
      Component(const Component& other)
        : Parent(other),
          m_fn(other.m_fn->copy()),
          m_i(other.m_i),
          m_j(other.m_j)
      {}

      constexpr
      Component(Component&& other)
        : Parent(std::move(other)),
          m_fn(std::move(other.m_fn)),
          m_i(std::move(other.m_i)),
          m_j(std::move(other.m_j))
      {}

      /**
       * @brief Gets the underlying matrix function.
       * @returns Reference to the operand function
       */
      constexpr
      const OperandType& getOperand() const
      {
        assert(m_fn);
        return *m_fn;
      }

      /**
       * @brief Evaluates the matrix entry at a point.
       * @param p Point at which to evaluate
       * @returns Scalar value @f$ A_{ij}(p) @f$
       */
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return this->getOperand().getValue(p).coeff(m_i, m_j);
      }

      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope& geom) const noexcept
      {
        return getOperand().getOrder(geom);
      }

      Component* copy() const noexcept override
      {
        return new Component(*this);
      }

    private:
      std::unique_ptr<OperandType> m_fn;
      const size_t m_i;
      const size_t m_j;
  };

  /**
   * @brief Deduction guide for matrix entry extraction.
   */
  template <class OperandDerived>
  Component(const FunctionBase<OperandDerived>&, size_t, size_t)
    -> Component<FunctionBase<OperandDerived>, size_t, size_t>;

  /**
   * @brief Extracts a component from a vector-valued GridFunction.
   *
   * Specialized component extraction for discrete grid functions that
   * maintains the reference to the original GridFunction.
   *
   * @tparam FES Finite element space type
   * @tparam Data Data storage type
   */
  template <class FES, class Data>
  class Component<GridFunction<FES, Data>> final
    : public RealFunctionBase<Component<GridFunction<FES, Data>>>
  {
    public:
      using OperandType = GridFunction<FES, Data>;

      using OperandRangeType = typename FormLanguage::Traits<OperandType>::RangeType;

      using ScalarType = typename FormLanguage::Traits<OperandRangeType>::ScalarType;

      using RangeType = ScalarType;

      using Parent = RealFunctionBase<Component<OperandType>>;

      /**
       * @brief Constructs component extractor for GridFunction.
       * @param u Vector-valued GridFunction
       * @param component Zero-based component index
       */
      constexpr
      Component(OperandType& u, size_t component)
        : m_u(u), m_idx(component)
      {}

      constexpr
      Component(const Component& other)
        : Parent(other),
          m_u(other.m_u),
          m_idx(other.m_idx)
      {}

      constexpr
      Component(Component&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)),
          m_idx(std::move(other.m_idx))
      {}

      /**
       * @brief Gets the underlying GridFunction.
       * @returns Reference to the grid function
       */
      constexpr
      OperandType& getGridFunction()
      {
        return m_u.get();
      }

      /**
       * @brief Gets the underlying GridFunction (const version).
       * @returns Const reference to the grid function
       */
      constexpr
      const OperandType& getGridFunction() const
      {
        return m_u.get();
      }

      /**
       * @brief Evaluates the component at a point.
       * @param p Point at which to evaluate
       * @returns Scalar value of the component
       */
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return m_u.get().getValue(p).coeff(m_idx);
      }

      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope& geom) const noexcept
      {
        return getGridFunction().getOrder(geom);
      }

      Component* copy() const noexcept override
      {
        return new Component(*this);
      }

    private:
      std::reference_wrapper<OperandType> m_u;
      const size_t m_idx;
  };

  /**
   * @brief Deduction guide for GridFunction component extraction.
   */
  template <class FES, class Data>
  Component(GridFunction<FES, Data>&, size_t) -> Component<GridFunction<FES, Data>>;

  /**
   * @brief Extracts a component from a vector-valued ShapeFunction.
   *
   * Extracts a scalar component from vector trial or test functions in
   * finite element formulations.
   *
   * @tparam OperandDerived Type of the shape function
   * @tparam FES Finite element space type
   * @tparam Space Trial or test function space
   */
  template <class OperandDerived, class FES, ShapeFunctionSpaceType Space>
  class Component<ShapeFunctionBase<OperandDerived, FES, Space>> final
    : public ShapeFunctionBase<Component<ShapeFunctionBase<OperandDerived, FES, Space>>, FES, Space>
  {
    public:
      using FESType = FES;
      static constexpr const ShapeFunctionSpaceType SpaceType = Space;

      using OperandType = ShapeFunctionBase<OperandDerived, FESType, SpaceType>;

      using OperandRangeType = typename FormLanguage::Traits<OperandType>::RangeType;

      using ScalarType = typename FormLanguage::Traits<OperandRangeType>::ScalarType;

      using RangeType = ScalarType;

      using Parent = ShapeFunctionBase<Component<OperandType>, FES, Space>;

      static_assert(std::is_same_v<OperandRangeType, Math::Vector<ScalarType>>);

      /**
       * @brief Constructs component extractor for ShapeFunction.
       * @param u Vector-valued trial or test function
       * @param component Zero-based component index
       */
      Component(const OperandType& u, size_t component)
        : Parent(u.getFiniteElementSpace()),
          m_u(u.copy()),
          m_idx(component)
      {}

      Component(const Component& other)
        : Parent(other),
          m_u(other.m_u->copy()),
          m_idx(other.m_idx)
      {}

      Component(Component&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)),
          m_idx(std::move(other.m_idx))
      {}

      /**
       * @brief Gets the component index.
       * @returns Component index
       */
      constexpr
      size_t getIndex() const
      {
        return m_idx;
      }

      /**
       * @brief Gets the underlying shape function.
       * @returns Reference to the operand
       */
      constexpr
      const OperandType& getOperand() const
      {
        assert(m_u);
        return *m_u;
      }

      /**
       * @brief Gets the leaf (underlying trial/test function).
       * @returns Reference to the leaf function
       */
      constexpr
      const auto& getLeaf() const
      {
        return getOperand().getLeaf();
      }

      /**
       * @brief Gets number of degrees of freedom on a polytope.
       * @param polytope Mesh polytope
       * @returns Number of DOFs
       */
      constexpr
      size_t getDOFs(const Geometry::Polytope& polytope) const
      {
        const size_t d = polytope.getDimension();
        const size_t i = polytope.getIndex();
        return this->getFiniteElementSpace().getFiniteElement(d, i).getCount();
      }

      /**
       * @brief Gets the current evaluation point.
       * @returns Reference to the point
       */
      const IntegrationPoint& getIntegrationPoint() const
      {
        return m_u->getIntegrationPoint();
      }

      Component& setIntegrationPoint(const IntegrationPoint& ip)
      {
        m_u->setIntegrationPoint(ip);
        return *this;
      }

      /**
       * @brief Gets the basis function value for local DOF.
       * @param local Local DOF index
       * @returns Scalar basis function value
       */
      constexpr
      auto getBasis(size_t local) const
      {
        return this->object(this->getOperand().getBasis(local)).coeff(m_idx);
      }

      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope& geom) const noexcept
      {
        return getOperand().getOrder(geom);
      }

      Component* copy() const noexcept override
      {
        return new Component(*this);
      }

    private:
      std::unique_ptr<OperandType> m_u;
      const size_t m_idx;
  };

  /**
   * @brief Deduction guide for ShapeFunction component extraction.
   */
  template <class OperandDerived, class FES, ShapeFunctionSpaceType Space>
  Component(const ShapeFunctionBase<OperandDerived, FES, Space>&, size_t)
    -> Component<ShapeFunctionBase<OperandDerived, FES, Space>>;
}

#endif

