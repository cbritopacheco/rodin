#ifndef RODIN_MPI_IO_MFEM_H
#define RODIN_MPI_IO_MFEM_H

#include "Rodin/IO/MFEM.h"

#include "Rodin/MPI/Context.h"

namespace Rodin::IO
{
  template <class Range>
  class GridFunctionPrinter<FileFormat::MFEM, Variational::P1<Range, Geometry::Mesh<Context::MPI>>>
    : public GridFunctionPrinterBase<Variational::P1<Range, Geometry::Mesh<Context::MPI>>>
  {
    public:
      using FESType = Variational::P1<Range, Geometry::Mesh<Context::MPI>>;

      using ObjectType = Variational::GridFunction<FESType>;

      using Parent = GridFunctionPrinterBase<FESType>;

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
        const auto& matrix = gf.getData();
        const Real* data = matrix.data();
        assert(matrix.size() >= 0);
        for (size_t i = 0; i < static_cast<size_t>(matrix.size()); i++)
          os << data[i] << '\n';
      }
  };
}
#endif
