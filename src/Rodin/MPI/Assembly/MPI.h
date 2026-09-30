/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MPI_ASSEMBLY_MPI_H
#define RODIN_MPI_ASSEMBLY_MPI_H

/**
 * @file
 * @brief MPI assembly iterators and assembly specializations.
 */

#include <algorithm>
#include <utility>
#include <vector>

#include <boost/mpi/collectives.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>

#include "Rodin/Variational/Integrator.h"
#include "Rodin/Variational/IntegrationPoint.h"

#include "Rodin/MPI/Geometry/Mesh.h"
#include "Rodin/Assembly/AssemblyBase.h"
#include "BoundaryDOFs.h"

namespace Rodin::Assembly
{
  /**
   * @brief Iteration helper over a region of an MPI mesh shard.
   *
   * This type wraps a distributed mesh and a region descriptor and provides
   * a @ref Rodin::Geometry::PolytopeIterator suitable for local integration
   * loops.
   *
   * MPI assembly uses shard-local topological indices throughout. A polytope
   * returned by this iterator belongs to the distributed mesh object, but its
   * `(d, i)` pair addresses the rank-local shard. The same local pair is passed
   * to finite element spaces (`getFiniteElement()`, `getDOFs()`) and to mesh
   * geometry/quadrature routines. Owned-entity filters decide which rank
   * contributes a row/entity; ghost entities can still appear where off-process
   * trial columns are needed.
   */
  class MPIIteration
  {
    public:
      /**
       * @brief Distributed mesh type iterated by this helper.
       */
      using MeshType = Geometry::Mesh<Context::MPI>;

      /**
       * @brief Constructs an iteration helper on a mesh region.
       * @param[in] mesh Distributed mesh.
       * @param[in] region Region descriptor to iterate.
       */
      MPIIteration(const MeshType& mesh, Geometry::Region region);

      /**
       * @brief Builds an iterator over the configured region.
       * @return Iterator over local polytopes in the region.
       */
      Geometry::PolytopeIterator getIterator() const;

    private:
      std::reference_wrapper<const MeshType> m_mesh;
      Geometry::Region m_region;
  };
}

namespace Rodin::Assembly
{
  /**
   * @brief Primary template declaration of the MPI assembly executor.
   *
   * Concrete behavior is provided by partial specializations for supported
   * operand types.
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref MPI "MPI<IndexMap<Scalar>, DirichletBC<TrialFunction, FunctionBase>>" | Distributed scalar Dirichlet elimination-map assembly on MPI mesh shards. |
   * | @ref MPI "MPI<IndexMap<pair<IndexArray, Vector>>, DirichletBC<TrialFunction, ShapeFunctionBase>>" | Distributed identification-style Dirichlet constraint assembly. |
   * | @ref MPI "MPI<Vec, LinearForm<FES, Vec>>" | PETSc/MPI linear-form assembly. |
   * | @ref MPI "MPI<Mat, BilinearForm<..., Mat>>" | PETSc/MPI bilinear-form assembly. |
   * | @ref MPI "MPI<PETSc::Math::LinearSystem, Problem<PETSc::Math::LinearSystem, U, V>>" | PETSc/MPI single-field problem assembly. |
   * | @ref MPI "MPI<PETSc::Math::LinearSystem, Problem<PETSc::Math::LinearSystem, U1, U2, U3, Us...>>" | PETSc/MPI mixed multi-field problem assembly. |
   */
  template <class LinearAlgebraType, class Operand>
  class MPI;

  template <class Scalar, class Solution, class FES, class ValueDerived>
  /**
   * @brief MPI assembler for scalar-valued Dirichlet boundary data.
   *
   * Produces the distributed map of constrained indices and their prescribed
   * values for a trial function on an MPI mesh shard.
   */
  class MPI<IndexMap<Scalar>,
    Variational::DirichletBC<Variational::TrialFunction<Solution, FES>,
      Variational::FunctionBase<ValueDerived>>>
    final : public AssemblyBase<IndexMap<Scalar>,
              Variational::DirichletBC<Variational::TrialFunction<Solution, FES>,
                Variational::FunctionBase<ValueDerived>>>
  {
    public:
      /**
       * @brief Finite-element-space type attached to the trial function.
       */
      using FESType =
        FES;

      /**
       * @brief Trial-function type of the assembled Dirichlet term.
       */
      using TrialFunctionType =
        Variational::TrialFunction<Solution, FES>;

      /**
       * @brief Boundary value function base type.
       */
      using ValueType =
        Variational::FunctionBase<ValueDerived>;

      /**
       * @brief Concrete Dirichlet boundary-condition operand type.
       */
      using DirichletBCType =
        Variational::DirichletBC<TrialFunctionType, ValueType>;

      /**
       * @brief Parent assembly base specialization.
       */
      using Parent =
        AssemblyBase<IndexMap<Scalar>, DirichletBCType>;

      /**
       * @brief Value range type induced by the finite element space.
       */
      using FESRangeType =
        typename FormLanguage::Traits<FESType>::RangeType;

      /**
       * @brief Input payload type consumed by execute().
       */
      using InputType =
        typename Parent::InputType;

      /**
       * @brief Default-constructs the MPI assembler.
       */
      MPI() = default;

      /**
       * @brief Copy-constructs the MPI assembler.
       */
      MPI(const MPI& other)
        : Parent(other)
      {}

      /**
       * @brief Move-constructs the MPI assembler.
       */
      MPI(MPI&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Assembles distributed Dirichlet contributions into an index map.
       *
       * Visits selected boundary functionals in the existing MPI mesh overlap,
       * filters by essential boundary attributes, evaluates the boundary value
       * pullback, and inserts global constrained indices into @p res. P1/H1 constraints
       * are evaluated locally, including on ghost faces. Only globally
       * supported P0g requires communication.
       *
       * @param[out] res   Target distributed index map.
       * @param[in] input  Assembly input wrapper carrying operand and value.
       */
      void execute(IndexMap<Scalar>& res, const InputType& input) const override
      {
        const auto& fes = input.getOperand().getFiniteElementSpace();
        MPIBoundaryDOFs<FES>(fes, input.getEssentialBoundary())
          .assemble(res, input.getValue());
      }

      /**
       * @brief Creates a polymorphic copy of this assembler.
       * @return Heap-allocated copy.
       */
      MPI* copy() const noexcept override
      {
        return new MPI(*this);
      }
  };

  /**
   * @brief MPI assembler for the identification Dirichlet BC `u = A(v)`.
   *
   * Visits required DOFs on halo boundary entities, filters by essential boundary
   * attributes, and evaluates the slave DOF functional on each master basis
   * expression with a live @c IntegrationPoint. Coefficients are never
   * thresholded or compared to identify distributed copies.
   * The source face is selected by distributed index. No constraint exchange
   * is required for locally supported spaces.
   */
  template <class Scalar, class Sol1, class FES1, class Derived2, class FES2,
    Variational::ShapeFunctionSpaceType Sp>
  class MPI<IndexMap<std::pair<IndexArray, Math::Vector<Scalar>>>,
    Variational::DirichletBC<Variational::TrialFunction<Sol1, FES1>,
      Variational::ShapeFunctionBase<Derived2, FES2, Sp>>>
    final : public AssemblyBase<IndexMap<std::pair<IndexArray, Math::Vector<Scalar>>>,
              Variational::DirichletBC<Variational::TrialFunction<Sol1, FES1>,
                Variational::ShapeFunctionBase<Derived2, FES2, Sp>>>
  {
    public:
      /// @brief Output map containing slave DOF indices and master coefficients.
      using OutputType = IndexMap<std::pair<IndexArray, Math::Vector<Scalar>>>;
      /// @brief Trial function type receiving the Dirichlet constraint.
      using TrialFunctionType = Variational::TrialFunction<Sol1, FES1>;
      /// @brief Shape-function expression used as the boundary value.
      using ValueType = Variational::ShapeFunctionBase<Derived2, FES2, Sp>;
      /// @brief Concrete Dirichlet boundary-condition type.
      using DirichletBCType = Variational::DirichletBC<TrialFunctionType, ValueType>;
      /// @brief Parent class type.
      using Parent = AssemblyBase<OutputType, DirichletBCType>;
      /// @brief Input payload type consumed by execute().
      using InputType = typename Parent::InputType;

      /// @brief Default constructor.
      MPI() = default;
      /// @brief Copy constructor.
      MPI(const MPI& other)
        : Parent(other)
      {}
      /// @brief Move constructor.
      MPI(MPI&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Assembles distributed identification constraints.
       * @param[out] res Target map from slave DOFs to master DOFs and weights.
       * @param[in] input Assembly input wrapper carrying operand and boundary data.
       */
      void execute(OutputType& res, const InputType& input) const override
      {
        const auto& fesU = input.getOperand().getFiniteElementSpace();
        auto& Av = const_cast<ValueType&>(input.getShapeFunction());
        const auto& fesV = Av.getLeaf().getFiniteElementSpace();
        const size_t faceDim = fesU.getMesh().getDimension() - 1;
        const MPIBoundaryDOFs<FES1> boundary(fesU, input.getEssentialBoundary());
        using Entries = std::vector<std::pair<Index, Scalar>>;
        IndexMap<Entries> rows;
        for (const auto& [slave, indices] : boundary.getDOFs())
        {
          const auto [face, local] = indices;
          const auto& feU = fesU.getFiniteElement(faceDim, face);
          const auto& feV = fesV.getFiniteElement(faceDim, face);
          const auto masterDOFs = fesV.getDOFs(faceDim, face);
          Entries entries;
          for (Index j = 0; j < static_cast<Index>(feV.getCount()); ++j)
          {
            auto basis = [&Av, j](const Geometry::Point& p) {
              const Variational::IntegrationPoint ip(p);
              Av.setIntegrationPoint(ip);
              return Av.getBasis(j);
            };
            const auto mapping = fesU.getPullback({faceDim, face}, std::move(basis));
            const Scalar coefficient = feU.getLinearForm(local)(mapping);
            if (coefficient != Scalar(0))
              entries.emplace_back(masterDOFs[j], coefficient);
          }
          rows.emplace(slave, std::move(entries));
        }
        boundary.synchronize(rows);
        res.clear();
        for (const auto& [slave, entries] : rows)
        {
          const Index n = static_cast<Index>(entries.size());
          IndexArray masters(n);
          Math::Vector<Scalar> coefficients(n);
          for (Index k = 0; k < n; ++k)
          {
            masters[k] = entries[k].first;
            coefficients[k] = entries[k].second;
          }
          res.emplace(slave, std::pair{std::move(masters), std::move(coefficients)});
        }
      }

      /** Evaluates affine data through the same halo-aware functional selection. */
      template <class Function>
      void assembleValues(IndexMap<Scalar>& values, const FES1& fes,
        const FlatSet<Geometry::Attribute>& attributes, const Function& function) const
      {
        MPIBoundaryDOFs<FES1>(fes, attributes).assemble(values, function);
      }

      /// @brief Creates a polymorphic copy of this assembler.
      MPI* copy() const noexcept override
      {
        return new MPI(*this);
      }
  };
}

#endif
