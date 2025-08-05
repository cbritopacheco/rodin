#ifndef RODIN_VARIATIONAL_P1_MFEM_H
#define RODIN_VARIATIONAL_P1_MFEM_H

#include "Rodin/IO/MFEM.h"

namespace Rodin::IO
{
  template <class Range, class Context, class Scalar>
  class GridFunctionPrinter<
    FileFormat::MFEM, Variational::P1<Range, Geometry::Mesh<Context>>, Math::Vector<Scalar>>
    : public GridFunctionPrinterBase<FileFormat::MFEM, Variational::P1<Range, Geometry::Mesh<Context>>, Math::Vector<Scalar>>
  {
    public:
      using RangeType = Range;

      using ScalarType = Scalar;

      using ContextType = Context;

      using MeshType = Geometry::Mesh<ContextType>;

      using FESType = Variational::P1<Range, MeshType>;

      static constexpr FileFormat Format = FileFormat::MFEM;

      using DataType = Math::Vector<ScalarType>;

      using Parent = GridFunctionPrinterBase<Format, FESType, DataType>;

      using Parent::Parent;

      void printHeader(std::ostream& os) override
      {
        const auto& gf = this->getObject();
        const auto& fes = gf.getFiniteElementSpace();
        os << "FiniteElementSpace\n"
           << "FiniteElementCollection: " << "H1_" << fes.getMesh().getDimension() << "D_P1\n"
           << "VDim: " << fes.getVectorDimension() << '\n'
           << "Ordering: " << MFEM::Ordering::VectorDimension
           << "\n\n";
      }
  };

}

#endif
