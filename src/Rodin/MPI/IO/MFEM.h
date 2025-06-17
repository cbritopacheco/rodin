#ifndef RODIN_MPI_IO_MFEM_H
#define RODIN_MPI_IO_MFEM_H

#include "Rodin/IO/MFEM.h"

#include "Rodin/MPI/Context.h"

namespace Rodin::IO
{
  template <class Range, class Data>
  class GridFunctionPrinter<FileFormat::MFEM, Data, Variational::P1<Range, Geometry::Mesh<Context::MPI>>>
    : public GridFunctionPrinterBase<Variational::P1<Range, Geometry::Mesh<Context::MPI>>, Data>
  {
    public:
      using FESType = Variational::P1<Range, Geometry::Mesh<Context::MPI>>;
      using DataType = Data;

      using ObjectType = Variational::GridFunction<FESType, DataType>;

      using Parent = GridFunctionPrinterBase<FESType, DataType>;

      GridFunctionPrinter(const ObjectType& gf)
        : Parent(gf)
      {}

      void print(std::ostream& os) override
      {
        const auto& gf = this->getObject();
        const auto& fes = gf.getFiniteElementSpace();
        os << "FiniteElementSpace\n"
           << "FiniteElementCollection: " << "H1_" << fes.getMesh().getDimension() << "D_P1\n"
           << "VDim: " << fes.getVectorDimension() << '\n'
           << "Ordering: " << MFEM::Ordering::VectorDimension
           << "\n\n";
        for (size_t i = 0; i < static_cast<size_t>(gf.getSize()); i++)
          os << gf[i] << '\n';
      }
  };
}
#endif
