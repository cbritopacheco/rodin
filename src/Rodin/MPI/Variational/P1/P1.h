#ifndef RODIN_MPI_VARIATIONAL_P1_P1_H
#define RODIN_MPI_VARIATIONAL_P1_P1_H

#include "Rodin/Variational/P1/P1.h"

#include "Rodin/MPI/Variational/FiniteElementSpace.h"
#include "Rodin/MPI/MPIMesh.h"

namespace Rodin::Variational
{
  template <class Range>
  class P1<Range, Geometry::Mesh<Context::MPI>> final
    : public FiniteElementSpace<
        Geometry::Mesh<Context::MPI>, P1<Real, Geometry::Mesh<Context::MPI>>>
  {
    public:
      /// Represents the Context of the P1 space
      using ContextType = Context::MPI;

      using ShardType = P1<Range, Geometry::Mesh<Context::Local>>;

      using ScalarType = typename ShardType::ScalarType;

      /// Range type of value
      using RangeType = typename ShardType::RangeType;

      /// Type of mesh on which the finite element space is built
      using MeshType = Geometry::Mesh<ContextType>;

      /// Type of finite element
      using ElementType = P1Element<RangeType>;

      /// Parent class
      using Parent = FiniteElementSpace<MeshType, P1<RangeType, MeshType>>;

    private:
      ShardType m_fes;
  };
}

#endif
