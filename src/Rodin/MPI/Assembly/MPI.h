#ifndef RODIN_MPI_ASSEMBLY_MPI_H
#define RODIN_MPI_ASSEMBLY_MPI_H

/**
 * @file
 * @brief MPI assembly iterators and assembly specializations.
 */

#include "Rodin/Variational/Integrator.h"
#include "Rodin/Variational/IntegrationPoint.h"

#include "Rodin/MPI/Geometry/Mesh.h"
#include "Rodin/Assembly/AssemblyBase.h"

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
      MPIIteration(const MeshType& mesh, Geometry::Region);

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
   */
  template <class LinearAlgebraType, class Operand>
  class MPI;

  template <class Scalar, class Solution, class FES, class ValueDerived>
  class MPI<
    IndexMap<Scalar>,
    Variational::DirichletBC<
      Variational::TrialFunction<Solution, FES>, Variational::FunctionBase<ValueDerived>>> final
    : public AssemblyBase<
        IndexMap<Scalar>,
        Variational::DirichletBC<
          Variational::TrialFunction<Solution, FES>, Variational::FunctionBase<ValueDerived>>>
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
       * Iterates over locally owned boundary entities of the MPI mesh shard,
       * filters by essential boundary attributes, evaluates the boundary value
       * pullback, and inserts global constrained indices into @p res.
       *
       * @param[out] res   Target distributed index map.
       * @param[in] input  Assembly input wrapper carrying operand and value.
       */
      void execute(IndexMap<Scalar>& res, const InputType& input) const override
      {
        const auto& u = input.getOperand();
        const auto& value = input.getValue();
        const auto& essBdr = input.getEssentialBoundary();
        const auto& fes = u.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        const size_t faceDim = mesh.getDimension() - 1;

        for (auto it = mesh.getBoundary(); it; ++it)
        {
          const auto& face = *it;
          const Index i = face.getIndex();

          if (!essBdr.empty())
          {
            const auto a = face.getAttribute();
            if (!a || !essBdr.contains(*a))
              continue;
          }

          const auto& fe = fes.getFiniteElement(faceDim, i);
          const auto mapping = fes.getPullback({ faceDim, i }, value);
          const auto& dofs = fes.getDOFs(faceDim, i);

          for (Index local = 0; local < fe.getCount(); ++local)
          {
            const Index global = dofs[local];
            auto pos = res.find(global);
            if (pos == res.end())
            {
              const auto s = fe.getLinearForm(local)(mapping);
              res.emplace_hint(pos, global, s);
            }
          }
        }
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
   * Iterates the locally owned boundary entities of the MPI mesh shard,
   * filters by essential boundary attributes, and emits, for each slave DOF
   * on a shared face, the master-DOF/coefficient pair
   * @c (masterDOFs[local], 1.0) — exact, tolerance-free pairing via the FES's
   * own DOF connectivity. Generalises to non-trivial @f$ A @f$ via
   * @c Av.setIntegrationPoint(IntegrationPoint(p)).getBasis(j).
   */
  template <class Scalar, class Sol1, class FES1,
            class Derived2, class FES2,
            Variational::ShapeFunctionSpaceType Sp>
  class MPI<
    IndexMap<std::pair<IndexArray, Math::Vector<Scalar>>>,
    Variational::DirichletBC<
      Variational::TrialFunction<Sol1, FES1>,
      Variational::ShapeFunctionBase<Derived2, FES2, Sp>>> final
    : public AssemblyBase<
        IndexMap<std::pair<IndexArray, Math::Vector<Scalar>>>,
        Variational::DirichletBC<
          Variational::TrialFunction<Sol1, FES1>,
          Variational::ShapeFunctionBase<Derived2, FES2, Sp>>>
  {
    public:
      using OutputType = IndexMap<std::pair<IndexArray, Math::Vector<Scalar>>>;
      using TrialFunctionType = Variational::TrialFunction<Sol1, FES1>;
      using ValueType = Variational::ShapeFunctionBase<Derived2, FES2, Sp>;
      using DirichletBCType =
        Variational::DirichletBC<TrialFunctionType, ValueType>;
      using Parent = AssemblyBase<OutputType, DirichletBCType>;
      using InputType = typename Parent::InputType;

      MPI() = default;
      MPI(const MPI& other) : Parent(other) {}
      MPI(MPI&& other) : Parent(std::move(other)) {}

      void execute(OutputType& res, const InputType& input) const override
      {
        const auto& u = input.getOperand();
        auto& Av = const_cast<ValueType&>(input.getShapeFunction());
        const auto& essBdr = input.getEssentialBoundary();
        const auto& fes_u = u.getFiniteElementSpace();
        const auto& fes_v = Av.getLeaf().getFiniteElementSpace();
        const auto& mesh  = fes_u.getMesh();
        const size_t faceDim = mesh.getDimension() - 1;

        for (auto it = mesh.getBoundary(); it; ++it)
        {
          const auto& face = *it;
          const Index fi = face.getIndex();

          if (!essBdr.empty())
          {
            const auto a = face.getAttribute();
            if (!a || !essBdr.contains(*a)) continue;
          }

          const auto& fe_u = fes_u.getFiniteElement(faceDim, fi);
          const auto& fe_v = fes_v.getFiniteElement(faceDim, fi);
          const auto& slaveDOFs  = fes_u.getDOFs(faceDim, fi);
          const auto& masterDOFs = fes_v.getDOFs(faceDim, fi);

          const Index nMasters = static_cast<Index>(fe_v.getCount());

          for (Index s = 0;
               s < static_cast<Index>(fe_u.getCount());
               s++)
          {
            const Index slave = slaveDOFs[s];
            auto pos = res.find(slave);
            if (pos != res.end()) continue;

            std::vector<Index>  mIdx;
            std::vector<Scalar> mCoef;
            mIdx.reserve(static_cast<size_t>(nMasters));
            mCoef.reserve(static_cast<size_t>(nMasters));

            for (Index j = 0; j < nMasters; j++)
            {
              auto basisCallable = [&Av, j]
                                   (const Geometry::Point& p)
              {
                Av.setIntegrationPoint(Variational::IntegrationPoint(p));
                return Av.getBasis(static_cast<size_t>(j));
              };
              const auto mapping =
                fes_u.getPullback({ faceDim, fi }, std::move(basisCallable));
              const Scalar c =
                static_cast<Scalar>(fe_u.getLinearForm(s)(mapping));
              if (c != Scalar(0))
              {
                mIdx.push_back(masterDOFs[j]);
                mCoef.push_back(c);
              }
            }

            if (mIdx.empty()) continue;

            const Index n = static_cast<Index>(mIdx.size());
            IndexArray masters(n);
            Math::Vector<Scalar> coeffs(n);
            for (Index k = 0; k < n; k++)
            {
              masters.coeffRef(k) = mIdx[static_cast<size_t>(k)];
              coeffs.coeffRef(k)  = mCoef[static_cast<size_t>(k)];
            }
            res.emplace_hint(pos, slave,
                std::pair{ std::move(masters), std::move(coeffs) });
          }
        }
      }

      MPI* copy() const noexcept override
      {
        return new MPI(*this);
      }
  };
}

#endif
