/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file DirichletBC.h
 * @brief Essential (Dirichlet) boundary conditions, in two flavours.
 *
 * Rodin's @ref DirichletBC encodes a strongly-imposed essential constraint on
 * a slave trial function @f$ u\in V_h^u @f$ over a subset
 * @f$ \Gamma_D \subset \mathcal{B}_h @f$ of the boundary. Two flavours share
 * the same template name and the same abstract base
 * @ref Rodin::Variational::DirichletBCBase, separated only by the type of the
 * second argument:
 *
 * 1. **Value-prescribing BC** — `DirichletBC(u, g)` with `g` a
 *    @ref FunctionBase (scalar, vector, or matrix-valued):
 *    @f[
 *      u(x) \;=\; g(x), \qquad x \in \Gamma_D.
 *    @f]
 *    This is the classical inhomogeneous Dirichlet condition, including the
 *    homogeneous case @f$ g\equiv 0 @f$.
 *
 * 2. **Identification BC** — `DirichletBC(u, A(v))` with `A(v)` a
 *    @ref ShapeFunctionBase expression linear in another (trial) shape
 *    function @f$ v\in V_h^v @f$:
 *    @f[
 *      u(x) \;=\; A(v)(x), \qquad x \in \Gamma_D,
 *    @f]
 *    where @f$ A @f$ is any operator producing a `ShapeFunctionBase`
 *    (the identity, a component selection @f$ v_x @f$, a left product
 *    @f$ f\,v @f$, a matrix product @f$ R\,v @f$, sums of these, etc.).
 *    This is *algebraic identification* of the two unknowns' boundary DOFs:
 *    both sides of the constraint involve degrees of freedom that are still
 *    unknown at solve time, not pre-evaluated values.
 *
 * # Mathematical model
 *
 * Let @f$ V_h^u @f$ have basis @f$ \{\varphi_i^u\}_{i=1}^{N_u} @f$ and dual
 * DOF functionals @f$ \{\ell_i^u\}_{i=1}^{N_u} @f$ satisfying
 * @f$ \ell_i^u(\varphi_k^u) = \delta_{ik} @f$. Likewise for @f$ V_h^v @f$.
 *
 * The DOF functionals are FE-local: for Lagrange they are nodal point
 * evaluations @f$ \ell_i^u(w) = w(x_i^u) @f$, for moment-based or
 * non-nodal elements they are integrals or pairings against dual moments.
 * Both flavours apply the same DOF functional on the slave side; they
 * differ only in what they apply it to.
 *
 * **Value-prescribing.** For each slave DOF @f$ s @f$ on @f$ \Gamma_D @f$
 * @f[
 *   \ell_s^u(u) \;=\; \ell_s^u(g) \quad\Longleftrightarrow\quad u_s = g_s,
 * @f]
 * where @f$ g_s := \ell_s^u(g) @f$. The assembled object is a map
 * @f$ s \mapsto g_s @f$ (the @ref Rodin::Variational::DirichletBCBase::ValueDOFs
 * alternative).
 *
 * **Identification.** For each slave DOF @f$ s @f$ on @f$ \Gamma_D @f$ and
 * with @f$ A @f$ linear,
 * @f[
 *   \ell_s^u(u) \;=\; \ell_s^u(A(v))
 *               \;=\; \sum_j \ell_s^u\!\bigl(A(\varphi_j^v)\bigr)\, v_j
 *               \;=\; \sum_j C_{sj}\, v_j,
 * @f]
 * with the constraint coefficients
 * @f[
 *   C_{sj} \;:=\; \ell_s^u\!\bigl(A(\varphi_j^v)\bigr).
 * @f]
 * The assembled object is a map @f$ s \mapsto (\{j_k\}, \{C_{s,j_k}\}) @f$
 * over the non-zero coefficients (the
 * @ref Rodin::Variational::DirichletBCBase::IdentifiedDOFs alternative).
 *
 * For Lagrange + same FES + @f$ A=\mathrm{id} @f$: the dual property gives
 * @f$ C_{sj}=\delta_{sj} @f$, so each slave maps to the same DOF index of
 * @f$ v @f$ with coefficient one. For @f$ A(v)=R(x)\,v @f$ on a vector FES
 * at a Lagrange node @f$ x_s @f$:
 * @f$ C_{(s,\alpha),(s,\beta)} = R_{\alpha\beta}(x_s) @f$ — a multi-master
 * row whose coefficients are exact entries of @f$ R @f$ at the slave node.
 *
 * # Exactness, by construction
 *
 * Slave/master DOF pairings are determined entirely by the FES's own
 * connectivity (`fes_u.getDOFs(faceDim, fi)` and
 * `fes_v.getDOFs(faceDim, fi)`). There is no geometric search and no
 * tolerance anywhere in the assembly — the FE's DOF-functional contract
 * (`fe_u.getLinearForm(s)` applied to a callable
 * `(Geometry::Point) -> value` through the slave-FES pullback) is the only
 * machinery used. Coefficient pruning uses strict `c != 0` so that
 * Lagrange-dual-induced sparsity is preserved bit-for-bit.
 *
 * The identification assembler invokes `Av.evaluateBasis(j, p)` on the
 * second argument: this is a direct point-evaluation path on
 * @ref ShapeFunctionBase that bypasses the
 * `setIntegrationPoint`/`getBasis` cache, so each per-point evaluation is
 * a fresh re-evaluation of the reference-element basis at
 * `p.getReferenceCoordinates()`.
 *
 * # Linear-system effect
 *
 * The same problem-level consumer loop handles both flavours by visiting
 * the variant returned from @ref DirichletBCBase::getDOFs:
 *
 * - For @ref DirichletBCBase::ValueDOFs at slave global index @f$ g_s @f$:
 *   @f[
 *     A_{g_s,g_s} \leftarrow 1, \quad
 *     A_{g_s,k}   \leftarrow 0 \;\forall k\neq g_s, \quad
 *     b_{g_s}     \leftarrow g_s\text{ value},
 *   @f]
 *   with the column at @f$ g_s @f$ moved to the RHS as
 *   @f$ b_{r} \mathrel{-}= A_{r,g_s}\cdot g_{s} @f$ before being zeroed
 *   (sparse path) or eliminated (dense path).
 *
 * - For @ref DirichletBCBase::IdentifiedDOFs at slave global index
 *   @f$ g_s @f$ with masters @f$ \{g_{m_k}\} @f$ and coefficients
 *   @f$ \{c_k\} @f$:
 *   @f[
 *     A_{g_s,g_s}     \leftarrow 1, \quad
 *     A_{g_s,g_{m_k}} \leftarrow -c_k, \quad
 *     b_{g_s}         \leftarrow 0,
 *   @f]
 *   with the column at @f$ g_s @f$ redirected to the master columns
 *   @f$ g_{m_k} @f$ scaled by @f$ c_k @f$ before being eliminated.
 *
 * # Boundary specification
 *
 * `.on(attr)` (or `.on(attr1, attr2, ...)`) selects the subset
 * @f$ \Gamma_D @f$ from the mesh's boundary attribute set. Both slave and
 * master DOFs are read from the *same* face polytopes — the assembler
 * never matches DOFs across distinct faces (no geometric pairing). For a
 * cross-face periodic relation use @ref Rodin::Variational::PeriodicBC,
 * which takes an explicit DOF adjacency map.
 *
 * # Usage examples
 *
 * @code{.cpp}
 * // Homogeneous Dirichlet BC: u = 0 on boundary
 * auto bc = DirichletBC(u, Zero());
 *
 * // Inhomogeneous Dirichlet BC on attribute 1
 * RealFunction g = [](const Point& p) { return sin(p.x()); };
 * auto bc = DirichletBC(u, g).on(1);
 *
 * // Identification: pin u on the boundary to v's DOFs
 * auto bc = DirichletBC(u, v).on(1);
 *
 * // Identification: u = R(x) v on the boundary
 * auto bc = DirichletBC(u, R * v).on(1);
 *
 * // Identification: scalar u equal to vector v's x-component
 * auto bc = DirichletBC(u, v.x()).on(1);
 * @endcode
 */
#ifndef RODIN_VARIATIONAL_DIRICHLETBC_H
#define RODIN_VARIATIONAL_DIRICHLETBC_H

#include <set>
#include <variant>

#include "Rodin/Utility.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/FormLanguage/List.h"

#include "Rodin/Assembly/ForwardDecls.h"

#include "ForwardDecls.h"

#include "Function.h"
#include "ShapeFunction.h"
#include "VectorFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup DirichletBCSpecializations DirichletBC Template Specializations
   * @brief Template specializations of the DirichletBC class.
   * @see DirichletBC
   */

  /**
   * @brief Abstract base class for a Dirichlet boundary condition.
   *
   * Encodes the essential constraint
   * @f[
   *   \mathrm{Operand} \;=\; \mathrm{Value}
   *   \quad\text{on}\quad
   *   \Gamma_D \subset \mathcal{B}_h
   * @f]
   * applied DOF-functional-wise: for each slave DOF @f$ s @f$ on
   * @f$ \Gamma_D @f$ the assembled object encodes a row of the form
   * @f[
   *   \ell_s^u(\text{Operand}) \;=\; \ell_s^u(\text{Value}),
   * @f]
   * which @ref Problem / @ref ProblemBody enforce strongly when assembling
   * the global linear system.
   *
   * The right-hand side of @f$ \ell_s^u(\text{Value}) @f$ has two
   * representations, exposed through a discriminated union via
   * @ref getDOFs:
   *
   * - **Value-prescribing** @f$ (\text{Value}=g) @f$: a scalar
   *   @f$ g_s := \ell_s^u(g) @f$ — known at assembly time.
   * - **Identification** @f$ (\text{Value}=A(v)) @f$: a sparse linear
   *   combination @f$ \sum_j C_{sj}\,v_j @f$ with coefficients
   *   @f$ C_{sj} := \ell_s^u(A(\varphi_j^v)) @f$ — coupling unknown DOFs of
   *   @f$ u @f$ to unknown DOFs of @f$ v @f$.
   *
   * Concrete subclasses fill exactly one alternative of @ref DOFs.
   *
   * @see DirichletBC
   * @see DirichletBC<TrialFunction<Sol,FES>, FunctionBase<...>>
   * @see DirichletBC<TrialFunction<Sol1,FES1>, ShapeFunctionBase<Derived2,FES2,Sp>>
   */
  template <class Scalar>
  class DirichletBCBase : public FormLanguage::Base
  {
    public:
      using ScalarType = Scalar;

      /**
       * @brief DOF map for the value-based Dirichlet condition `u = g`.
       *
       * Maps each constrained global FES DOF index to its prescribed scalar
       * value @f$ g(x_i) @f$.
       */
      using ValueDOFs = IndexMap<ScalarType>;

      /**
       * @brief DOF map for the identification-based Dirichlet condition
       *        `u = A(v)`.
       *
       * Maps each slave DOF (in @f$ u @f$'s FES) to a pair of arrays:
       *  - the master DOF indices in @f$ v @f$'s FES (FES-global, pre-offset),
       *  - the corresponding scalar coefficients @f$ C_{sj} = [A(\varphi_j^v)](x_s) @f$.
       *
       * The constraint encoded is @f$ u_s = \sum_j C_{sj}\, v_j @f$.
       */
      using IdentifiedDOFs =
        IndexMap<std::pair<IndexArray, Math::Vector<ScalarType>>>;

      /**
       * @brief Variant DOF type covering both Dirichlet semantics.
       *
       * The active alternative depends on which subclass produced it: a
       * value-prescribing BC writes the @ref ValueDOFs alternative, an
       * identification BC writes the @ref IdentifiedDOFs alternative.
       */
      using DOFs = std::variant<ValueDOFs, IdentifiedDOFs>;

      /**
       * @brief Assembles the Dirichlet boundary condition.
       *
       * Either fills the @ref ValueDOFs alternative (for value-prescribing
       * BCs) or the @ref IdentifiedDOFs alternative (for identification BCs).
       */
      virtual void assemble() = 0;

      /**
       * @brief Gets the map of constrained DOFs.
       * @return Variant holding either the ValueDOFs or the IdentifiedDOFs map.
       */
      virtual const DOFs& getDOFs() const = 0;

      /**
       * @brief Checks if this is a component-wise boundary condition.
       * @return True if BC applies to a single component of a vector function
       */
      virtual bool isComponent() const = 0;

      /**
       * @brief Gets the operand (trial function) of the BC.
       * @return Reference to the trial function being constrained
       */
      virtual const FormLanguage::Base& getOperand() const = 0;

      /**
       * @brief Gets the prescribed boundary value.
       * @return Reference to the function defining the BC value
       */
      virtual const FormLanguage::Base& getValue() const = 0;

      /**
       * @brief Returns the UUID of @f$ v @f$'s leaf trial function for an
       *        identification BC, or @c nullopt for a value-prescribing BC.
       *
       * Consumer assemblies use this to locate @f$ v @f$'s trial block and
       * apply the correct global offset when assembling identified-DOF
       * constraints.
       */
      virtual Optional<Identifiable::UUID> getValueUUID() const { return {}; }

      /**
       * @brief Creates a polymorphic copy of this BC.
       * @return Pointer to a new copy
       */
      virtual DirichletBCBase* copy() const noexcept override = 0;
  };

  /// Alias for a list of Dirichlet boundary conditions
  template <class Scalar>
  using EssentialBoundary = FormLanguage::List<DirichletBCBase<Scalar>>;

  /**
   * @ingroup DirichletBCSpecializations
   * @brief Value-prescribing Dirichlet boundary condition, @f$ u = g @f$.
   *
   * Imposes the inhomogeneous (or homogeneous, if @f$ g\equiv 0 @f$)
   * Dirichlet condition
   * @f[
   *   u(x) \;=\; g(x), \qquad x\in\Gamma_D \subset \mathcal{B}_h,
   * @f]
   * algebraically: for each slave DOF @f$ s @f$ of @f$ u @f$ on
   * @f$ \Gamma_D @f$,
   * @f[
   *   u_s \;=\; \ell_s^u(g) \;=\; g_s,
   * @f]
   * where @f$ \ell_s^u @f$ is the slave finite element's DOF functional —
   * point evaluation @f$ g(x_s) @f$ for Lagrange, an integral / moment
   * pairing for non-nodal elements. The assembler streams these scalars
   * into the @ref DirichletBCBase::ValueDOFs alternative of the variant
   * returned by @ref getDOFs.
   *
   * In the global linear system the slave row is replaced by an identity
   * row with right-hand side @f$ g_s @f$ and the slave column is moved to
   * the RHS:
   * @f[
   *   A_{g_s,g_s}\leftarrow 1,\;
   *   A_{g_s,k}\leftarrow 0\;\forall k\neq g_s,\;
   *   b_{g_s}\leftarrow g_s,\;
   *   b_r \mathrel{-}= A_{r,g_s}\cdot g_s,\;
   *   A_{r,g_s}\leftarrow 0,\; r\neq g_s.
   * @f]
   *
   * `.on(attr...)` selects the boundary attribute(s) for @f$ \Gamma_D @f$.
   *
   * @tparam Solution Solution type of the trial function being constrained
   * @tparam FES Finite element space of the trial function
   * @tparam ValueDerived CRTP-derived type of the value @f$ g @f$
   *
   * @see DirichletBCBase
   */
  template <class Solution, class FES, class ValueDerived>
  class DirichletBC<TrialFunction<Solution, FES>, FunctionBase<ValueDerived>> final
    : public DirichletBCBase<typename FormLanguage::Traits<FES>::ScalarType>
  {
    public:
      using FESType =
        FES;

      /// Operand type
      using OperandType =
        TrialFunction<Solution, FESType>;

      /// Scalar type
      using ScalarType =
        typename FormLanguage::Traits<FESType>::ScalarType;

      /// Parent class
      using Parent = DirichletBCBase<ScalarType>;

      /// Value DOF map type (the alternative populated by this specialization)
      using ValueDOFs = typename Parent::ValueDOFs;

      /// Variant DOFs type, inherited from Parent
      using DOFs = typename Parent::DOFs;

      /// Value type
      using ValueType =
        FunctionBase<ValueDerived>;

      using FESMeshType =
        typename FormLanguage::Traits<FESType>::MeshType;

      using FESRangeType =
        typename FormLanguage::Traits<FESType>::RangeType;

      using FESMeshContextType =
        typename FormLanguage::Traits<FESMeshType>::ContextType;

      using DefaultAssemblyType =
        typename Assembly::Default<FESMeshContextType>::template Type<ValueDOFs, DirichletBC>;

      using AssemblyType =
        DefaultAssemblyType;

      /**
       * @brief Constructs the object given the Operand and Value.
       * @param[in] u ShapeFunction object
       * @param[in] v Value object
       */
      DirichletBC(const OperandType& u, const ValueType& v)
        : m_u(u), m_value(v.copy())
      {}

      /**
       * @brief Copy constructor
       */
      DirichletBC(const DirichletBC& other)
        : Parent(other),
          m_u(other.m_u),
          m_value(other.m_value->copy()),
          m_essBdr(other.m_essBdr),
          m_dofs(other.m_dofs),
          m_assembly(other.m_assembly)
      {}

      /**
       * @brief Move constructor
       */
      DirichletBC(DirichletBC&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)),
          m_value(std::move(other.m_value)),
          m_essBdr(std::move(other.m_essBdr)),
          m_dofs(std::move(other.m_dofs)),
          m_assembly(std::move(other.m_assembly))
      {}

      /**
       * @brief Specifies the region of the boundary over which the condition
       * will be imposed.
       * @param[in] bdrAtr Attribute associated to the boundary region
       */
      constexpr
      DirichletBC& on(Geometry::Attribute bdrAtr)
      {
        return on(FlatSet<Geometry::Attribute>{bdrAtr});
      }

      template <class A1, class A2, class ... As>
      constexpr
      DirichletBC& on(A1 a1, A2 a2, As... as)
      {
        return on(FlatSet<Geometry::Attribute>{ a1, a2, as... });
      }

      /**
       * @brief Specifies the regions of the boundary over which the condition
       * will be imposed.
       * @param[in] bdrAttrs Attributes associated to the boundary regions
       */
      constexpr
      DirichletBC& on(const FlatSet<Geometry::Attribute>& bdrAttrs)
      {
        assert(bdrAttrs.size() > 0);
        m_essBdr = bdrAttrs;
        return *this;
      }

      /**
       * @returns Attributes over which the boundary condition is imposed.
       */
      constexpr
      const FlatSet<Geometry::Attribute>& getAttributes() const
      {
        return m_essBdr;
      }

      /**
       * @brief Computes the indices of the degrees of freedoms associated to
       * the boundary region.
       *
       * This will compute the degrees of freedom over the incidence set
       * @f[
       *  D - 1 \longrightarrow 0 ~.
       * @f]
       *
       * If the set of specified attributes is empty, this will
       * compute the degrees of freedom over the boundary, in which case the
       * incidence set
       * @f[
       *  D - 1 \longrightarrow D
       * @f]
       * is also required.
       */
      void assemble() override
      {
        // Ensure variant holds the ValueDOFs alternative.
        if (!std::holds_alternative<ValueDOFs>(m_dofs))
          m_dofs = ValueDOFs{};
        m_assembly.execute(
            std::get<ValueDOFs>(m_dofs), { m_u.get(), *m_value, m_essBdr });
      }

      bool isComponent() const override
      {
        return false;
      }

      const OperandType& getOperand() const override
      {
        return m_u.get();
      }

      const ValueType& getValue() const override
      {
        assert(m_value);
        return *m_value;
      }

      const DOFs& getDOFs() const override
      {
        return m_dofs;
      }

      const Assembly::AssemblyBase<ValueDOFs, DirichletBC>& getAssembly() const
      {
        assert(m_assembly);
        return *m_assembly;
      }

      DirichletBC* copy() const noexcept override
      {
        return new DirichletBC(*this);
      }

    private:
      std::reference_wrapper<const OperandType> m_u;
      std::unique_ptr<ValueType> m_value;
      FlatSet<Geometry::Attribute> m_essBdr;
      DOFs m_dofs{ValueDOFs{}};
      AssemblyType m_assembly;
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for DirichletBC
   * @tparam FES Type of finite element space
   * @tparam ValueDerived Derived type of FunctionBase
   */
  template <class Solution, class FES, class FunctionDerived>
  DirichletBC(const TrialFunction<Solution, FES>&, const FunctionBase<FunctionDerived>&)
    -> DirichletBC<TrialFunction<Solution, FES>, FunctionBase<FunctionDerived>>;

  /**
   * @ingroup DirichletBCSpecializations
   * @brief Identification Dirichlet boundary condition, @f$ u = A(v) @f$.
   *
   * Imposes a *linear-in-DOFs* identification of the slave trial function
   * @f$ u\in V_h^u @f$ with a shape-function expression @f$ A(v) @f$ over
   * @f$ \Gamma_D \subset \mathcal{B}_h @f$:
   * @f[
   *   u(x) \;=\; A(v)(x), \qquad x\in\Gamma_D,
   * @f]
   * where @f$ v\in V_h^v @f$ is another (trial) shape function and
   * @f$ A @f$ is any operator producing a @ref ShapeFunctionBase
   * (e.g. the identity @f$ A=\mathrm{id} @f$, a component @f$ A(v)=v_x @f$,
   * a left product @f$ A(v)=f\,v @f$, a matrix product @f$ A(v)=R\,v @f$,
   * sums of these, …).
   *
   * # Mathematical model
   *
   * Let @f$ \{\varphi_j^v\}_{j=1}^{N_v} @f$ be the basis of @f$ V_h^v @f$
   * and @f$ \{\ell_s^u\}_{s=1}^{N_u} @f$ be the dual DOF functionals of
   * @f$ V_h^u @f$ (so @f$ \ell_s^u(\varphi_k^u)=\delta_{sk} @f$). Linearity
   * of @f$ A @f$ in @f$ v @f$ gives, for any
   * @f$ v=\sum_j v_j\,\varphi_j^v @f$:
   * @f[
   *   A(v)(x) \;=\; \sum_j v_j\,A(\varphi_j^v)(x).
   * @f]
   *
   * Applying the slave DOF functional @f$ \ell_s^u @f$ on both sides of
   * @f$ u = A(v) @f$ yields the constraint row for each slave DOF
   * @f$ s @f$ on @f$ \Gamma_D @f$:
   * @f[
   *   \boxed{\;
   *     u_s \;=\; \ell_s^u(A(v))
   *           \;=\; \sum_j C_{sj}\, v_j,
   *     \quad\text{where}\quad
   *     C_{sj} \;:=\; \ell_s^u\!\bigl(A(\varphi_j^v)\bigr).
   *   \;}
   * @f]
   *
   * The assembler stores the per-slave row as the
   * @ref DirichletBCBase::IdentifiedDOFs alternative — a sparse list of
   * `(master DOF index, coefficient)` pairs, with strict @f$ C_{sj}=0 @f$
   * pruning so that Lagrange-dual sparsity (e.g. @f$ \delta_{sj} @f$ for
   * @f$ A=\mathrm{id} @f$ on the same FES) is preserved bit-for-bit.
   *
   * # Concrete examples
   *
   * - **Identity, same FES (Lagrange).**  @f$ A=\mathrm{id} @f$,
   *   @f$ V_h^u=V_h^v @f$. Then
   *   @f$ C_{sj} = \ell_s^u(\varphi_j^v) = \delta_{sj} @f$, so each slave
   *   reduces to one master with coefficient @f$ 1 @f$:
   *   @f$ u_s = v_s @f$.
   *
   * - **Scalar scaling.**  @f$ A(v)(x)=f(x)\,v(x) @f$ for a scalar field
   *   @f$ f @f$, Lagrange same FES, slave node @f$ x_s @f$:
   *   @f$ C_{sj} = f(x_s)\,\delta_{sj} @f$, hence
   *   @f$ u_s = f(x_s)\,v_s @f$.
   *
   * - **Matrix rotation, vector FES.**  @f$ A(v)(x)=R(x)\,v(x) @f$ with
   *   @f$ R\in\mathbb{R}^{d\times d} @f$, Lagrange vector @f$ V_h^v @f$,
   *   slave DOF @f$ (s,\alpha) @f$:
   *   @f[
   *     C_{(s,\alpha),(s,\beta)} \;=\; R_{\alpha\beta}(x_s),
   *   @f]
   *   so each scalar slave depends on the @f$ d @f$ component-DOFs of
   *   @f$ v @f$ at the same vertex with coefficients given by the @f$
   *   \alpha @f$-th row of @f$ R(x_s) @f$.
   *
   * - **Component extraction.**  @f$ A(v)=v_x @f$ for a vector
   *   @f$ v\in V_h^v @f$ and scalar @f$ u @f$. The
   *   @ref Component<ShapeFunctionBase<...>> propagation of
   *   @c evaluateBasis projects the chosen component, so each slave
   *   depends only on the @f$ x @f$-component DOFs of @f$ v @f$ with
   *   coefficient @f$ 1 @f$.
   *
   * # Assembly: how @f$ C_{sj} @f$ is computed exactly
   *
   * The identification assembler computes
   * @f$ C_{sj}=\ell_s^u(A(\varphi_j^v)) @f$ for every slave/master pair
   * on every face by feeding the slave finite element a callable
   * @f[
   *   p \;\mapsto\; A(\varphi_j^v)(p), \qquad p\in K \subset \Gamma_D,
   * @f]
   * via the slave-FES pullback, and then asking
   * @c fe_u.getLinearForm(s) to apply *its own* DOF functional to it. The
   * callable in turn invokes
   * @ref ShapeFunctionBase::evaluateBasis on @f$ A(v) @f$, which is the
   * direct point-evaluation path that bypasses the
   * @c setIntegrationPoint cache.
   *
   * The slave FE decides where to sample the callable — point evaluation
   * for Lagrange, an integral over the face for moment-based FEs,
   * directional samples for RT/Nedelec — and the assembler does not
   * branch on the FE type. This is what makes the implementation
   * **FES-generic**: the only contracts used are
   * `fes.getDOFs(face)`, `fes.getFiniteElement(face)`,
   * `fe.getLinearForm(s)`, and `fes.getPullback(face, callable)`.
   *
   * # Exactness, by construction
   *
   * Slave/master DOFs are paired by index on shared faces:
   * @c fes_u.getDOFs(faceDim, fi) and @c fes_v.getDOFs(faceDim, fi) on the
   * *same* face polytope. There is no geometric search, no tolerance, no
   * @f$ \varepsilon @f$ anywhere. Coefficient pruning uses strict
   * @c "c != 0" so that exact zeros from Lagrange duality (or from
   * component projections) are eliminated bit-for-bit. For a cross-face
   * periodic relation use @ref Rodin::Variational::PeriodicBC, which
   * takes an explicit DOF adjacency map.
   *
   * # Linear-system effect
   *
   * For each slave global index @f$ g_s @f$ with masters
   * @f$ \{g_{m_k}\} @f$ and coefficients @f$ \{c_k\} @f$, the consumer
   * loop replaces the slave row with the constraint row
   * @f[
   *   1\cdot u_{g_s} \;-\; \sum_k c_k\,u_{g_{m_k}} \;=\; 0,
   * @f]
   * and redirects every slave-column contribution
   * @f$ A_{r,g_s} @f$ into the master columns as
   * @f$ A_{r,g_{m_k}} \mathrel{+}= A_{r,g_s}\cdot c_k @f$, after which
   * @f$ A_{r,g_s} @f$ is zeroed.
   *
   * @tparam Solution Solution type of the slave trial @f$ u @f$
   * @tparam FES1 Finite element space of @f$ u @f$
   * @tparam Derived2 CRTP-derived type of @f$ A(v) @f$
   * @tparam FES2 Finite element space of @f$ v @f$
   * @tparam Sp Shape function space type (Trial or Test) of @f$ A(v) @f$
   *
   * @see DirichletBCBase
   * @see ShapeFunctionBase::evaluateBasis
   * @see PeriodicBC
   */
  template <class Solution, class FES1,
            class Derived2, class FES2, ShapeFunctionSpaceType Sp>
  class DirichletBC<TrialFunction<Solution, FES1>,
                    ShapeFunctionBase<Derived2, FES2, Sp>> final
    : public DirichletBCBase<typename FormLanguage::Traits<FES1>::ScalarType>
  {
    public:
      using FESType = FES1;

      /// Operand type (the slave trial function)
      using OperandType = TrialFunction<Solution, FESType>;

      /// Value type (the right-hand-side shape function expression)
      using ValueType = ShapeFunctionBase<Derived2, FES2, Sp>;

      /// Scalar type
      using ScalarType =
        typename FormLanguage::Traits<FESType>::ScalarType;

      /// Parent class
      using Parent = DirichletBCBase<ScalarType>;

      /// Identified-DOFs alternative populated by this specialization
      using IdentifiedDOFs = typename Parent::IdentifiedDOFs;

      /// Variant DOFs type
      using DOFs = typename Parent::DOFs;

      using FESMeshType =
        typename FormLanguage::Traits<FESType>::MeshType;

      using FESRangeType =
        typename FormLanguage::Traits<FESType>::RangeType;

      using FESMeshContextType =
        typename FormLanguage::Traits<FESMeshType>::ContextType;

      using DefaultAssemblyType =
        typename Assembly::Default<FESMeshContextType>::template Type<IdentifiedDOFs, DirichletBC>;

      using AssemblyType =
        DefaultAssemblyType;

      /**
       * @brief Constructs the identification Dirichlet BC.
       * @param[in] u Slave trial function
       * @param[in] v Right-hand-side shape function expression @f$ A(v) @f$
       */
      DirichletBC(const OperandType& u, const ValueType& v)
        : m_u(u), m_v(v.copy())
      {}

      /// Copy constructor
      DirichletBC(const DirichletBC& other)
        : Parent(other),
          m_u(other.m_u),
          m_v(other.m_v->copy()),
          m_essBdr(other.m_essBdr),
          m_dofs(other.m_dofs),
          m_assembly(other.m_assembly)
      {}

      /// Move constructor
      DirichletBC(DirichletBC&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)),
          m_v(std::move(other.m_v)),
          m_essBdr(std::move(other.m_essBdr)),
          m_dofs(std::move(other.m_dofs)),
          m_assembly(std::move(other.m_assembly))
      {}

      /**
       * @brief Sets @f$ \Gamma_D @f$ to a single boundary attribute.
       *
       * Selects the subset of @f$ \mathcal{B}_h @f$ over which the
       * identification @f$ u=A(v) @f$ is enforced. Both slave and master
       * DOFs are read from the *same* face polytopes — there is no
       * cross-face matching.
       */
      constexpr
      DirichletBC& on(Geometry::Attribute bdrAtr)
      {
        return on(FlatSet<Geometry::Attribute>{bdrAtr});
      }

      /**
       * @brief Sets @f$ \Gamma_D @f$ to the union of multiple attributes.
       */
      template <class A1, class A2, class ... As>
      constexpr
      DirichletBC& on(A1 a1, A2 a2, As... as)
      {
        return on(FlatSet<Geometry::Attribute>{ a1, a2, as... });
      }

      /**
       * @brief Sets @f$ \Gamma_D @f$ to a precomputed attribute set.
       */
      constexpr
      DirichletBC& on(const FlatSet<Geometry::Attribute>& bdrAttrs)
      {
        assert(bdrAttrs.size() > 0);
        m_essBdr = bdrAttrs;
        return *this;
      }

      /**
       * @returns The boundary attribute set defining @f$ \Gamma_D @f$.
       */
      constexpr
      const FlatSet<Geometry::Attribute>& getAttributes() const
      {
        return m_essBdr;
      }

      /**
       * @brief Computes the constraint rows
       *        @f$ u_s = \sum_j C_{sj}\,v_j @f$ for every slave DOF on
       *        @f$ \Gamma_D @f$.
       *
       * Iterates the boundary face polytopes whose attribute lies in
       * @ref getAttributes. For each face @f$ K @f$ and slave-local index
       * @f$ s @f$, computes
       * @f[
       *   C_{sj} \;=\; \ell_s^{u,K}\!\bigl(A(\varphi_j^{v,K})\bigr),
       *   \qquad j = 0,\dots,n_v(K)-1,
       * @f]
       * by handing the slave finite element the callable
       * @f$ p\mapsto A(\varphi_j^{v,K})(p) @f$ via the slave-FES
       * pullback and asking @c fe_u.getLinearForm(s) to apply its DOF
       * functional to it. The callable invokes
       * @ref ShapeFunctionBase::evaluateBasis on @f$ A(v) @f$, which
       * bypasses the @c setIntegrationPoint cache and re-evaluates the
       * reference-element basis at @c p.getReferenceCoordinates() each
       * call.
       *
       * Strict @f$ C_{sj}\neq 0 @f$ pruning is applied; the result is the
       * @ref DirichletBCBase::IdentifiedDOFs alternative of the variant
       * returned by @ref getDOFs.
       */
      void assemble() override
      {
        if (!std::holds_alternative<IdentifiedDOFs>(m_dofs))
          m_dofs = IdentifiedDOFs{};
        m_assembly.execute(
            std::get<IdentifiedDOFs>(m_dofs), { m_u.get(), *m_v, m_essBdr });
      }

      bool isComponent() const override
      {
        return false;
      }

      /**
       * @returns The slave trial function @f$ u @f$.
       */
      const OperandType& getOperand() const override
      {
        return m_u.get();
      }

      /**
       * @returns The right-hand-side shape-function expression @f$ A(v) @f$.
       */
      const ValueType& getValue() const override
      {
        assert(m_v);
        return *m_v;
      }

      /**
       * @returns The variant DOFs map; this specialization writes the
       *          @ref DirichletBCBase::IdentifiedDOFs alternative.
       */
      const DOFs& getDOFs() const override
      {
        return m_dofs;
      }

      /**
       * @returns The UUID of @f$ v @f$'s leaf trial function — used by
       *          consumer assemblies to locate @f$ v @f$'s trial block
       *          and apply the correct global offset to master DOFs.
       */
      Optional<Identifiable::UUID> getValueUUID() const override
      {
        assert(m_v);
        return m_v->getLeaf().getUUID();
      }

      const Assembly::AssemblyBase<IdentifiedDOFs, DirichletBC>& getAssembly() const
      {
        return m_assembly;
      }

      DirichletBC* copy() const noexcept override
      {
        return new DirichletBC(*this);
      }

    private:
      std::reference_wrapper<const OperandType> m_u;
      std::unique_ptr<ValueType> m_v;
      FlatSet<Geometry::Attribute> m_essBdr;
      DOFs m_dofs{IdentifiedDOFs{}};
      AssemblyType m_assembly;
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for the identification DirichletBC.
   */
  template <class Solution, class FES1,
            class Derived2, class FES2, ShapeFunctionSpaceType Sp>
  DirichletBC(const TrialFunction<Solution, FES1>&,
              const ShapeFunctionBase<Derived2, FES2, Sp>&)
    -> DirichletBC<TrialFunction<Solution, FES1>,
                   ShapeFunctionBase<Derived2, FES2, Sp>>;
}

#endif
