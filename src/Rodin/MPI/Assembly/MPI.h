#ifndef RODIN_MPI_ASSEMBLY_MPI_H
#define RODIN_MPI_ASSEMBLY_MPI_H

#include "Rodin/Variational/Integrator.h"

#include "Rodin/MPI/Geometry/Mesh.h"
#include "Rodin/Assembly/AssemblyBase.h"

namespace Rodin::Assembly
{
  class MPIIteration
  {
    public:
      using MeshType = Geometry::Mesh<Context::MPI>;

      MPIIteration(const MeshType& mesh, Variational::Integrator::Region);

      Geometry::PolytopeIterator getIterator() const;

    private:
      std::reference_wrapper<const MeshType> m_mesh;
      Variational::Integrator::Region m_region;
  };
}

namespace Rodin::Assembly
{
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

      MPI() = default;

      MPI(const MPI& other)
        : Parent(other)
      {}

      MPI(MPI&& other)
        : Parent(std::move(other))
      {}

      void execute(IndexMap<Scalar>& res, const InputType& input) const override
      {
        const auto& u = input.getOperand();
        const auto& value = input.getValue();
        const auto& essBdr = input.getEssentialBoundary();
        const auto& fes = u.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        const auto& shard = mesh.getShard();
        const size_t faceDim = shard.getDimension() - 1;
        const size_t faceCount = shard.getFaceCount();
        for (Index i = 0; i < faceCount; ++i)
        {
          if (shard.isGhost(faceDim, i))
            continue;
          if (shard.isBoundary(i))
          {
            if (essBdr.size() == 0 || essBdr.count(shard.getAttribute(faceDim, i)))
            {
              const auto& fe = fes.getShard().getFiniteElement(faceDim, i);
              const auto& mapping =
                fes.getShard().getMapping({ faceDim, i }, value.template cast<FESRangeType>());
              for (Index local = 0; local < fe.getCount(); local++)
              {
                const Index global = fes.getGlobalIndex(
                    fes.getShard().getGlobalIndex({ faceDim, i }, local));
                auto find = res.find(global);
                if (find == res.end())
                {
                  const auto& lf = fe.getLinearForm(local);
                  const auto s = lf(mapping);
                  res.insert(find, { global, s });
                }
              }
            }
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
