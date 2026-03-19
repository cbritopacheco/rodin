#ifndef RODIN_MPI_ASSEMBLY_MPI_H
#define RODIN_MPI_ASSEMBLY_MPI_H

/**
 * @file
 * @brief MPI assembly iterators and assembly specializations.
 */

#include "Rodin/Variational/Integrator.h"

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
   */
  class MPIIteration
  {
    public:
      using MeshType = Geometry::Mesh<Context::MPI>;

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
      using FESType =
        FES;

      using TrialFunctionType =
        Variational::TrialFunction<Solution, FES>;

      using ValueType =
        Variational::FunctionBase<ValueDerived>;

      using DirichletBCType =
        Variational::DirichletBC<TrialFunctionType, ValueType>;

      using Parent =
        AssemblyBase<IndexMap<Scalar>, DirichletBCType>;

      using FESRangeType =
        typename FormLanguage::Traits<FESType>::RangeType;

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
}

#endif
