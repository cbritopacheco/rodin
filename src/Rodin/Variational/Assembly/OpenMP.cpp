#ifdef RODIN_USE_OPENMP

#include <omp.h>

#include "Rodin/Variational/FiniteElementSpace.h"
#include "Rodin/Variational/BilinearFormIntegrator.h"

#include "OpenMP.h"

namespace Rodin::Variational::Assembly
{
  mfem::SparseMatrix
  OpenMP<BilinearFormBase<mfem::SparseMatrix>>
  ::execute(const Input& input) const
  {
    OperatorType result(input.testFES.getSize(), input.trialFES.getSize());
    result = 0.0;

    for (const auto& bfi : input.bfis)
    {
      switch (bfi.getRegion())
      {
        case Geometry::Region::Domain:
        {
          std::vector<OperatorType> local;
#pragma omp parallel
          {
            const int id = omp_get_thread_num();
#pragma omp master
            {
              const int nt = omp_get_num_threads();
              local.reserve(nt);
              for (int i = 0; i < nt; i++)
                local.emplace_back(input.testFES.getSize(), input.trialFES.getSize());
            }
#pragma omp barrier
#pragma omp for
            for (int i = 0; i < input.mesh.count<Geometry::Element>(); i++)
            {
              const auto& element = input.mesh.get<Geometry::Element>(i);
              if (bfi.getAttributes().size() == 0
                  || bfi.getAttributes().count(element.getAttribute()))
              {
                local[id].AddSubMatrix(
                    input.testFES.getDOFs(element),
                    input.trialFES.getDOFs(element),
                    bfi.getElementMatrix(element));
              }
            }
#pragma omp critical
            {
              result += local[id];
            }
          }
          break;
        }
        case Geometry::Region::Boundary:
        {
          assert(false);
          break;
        }
        case Geometry::Region::Interface:
        {
          assert(false);
          break;
        }
      }
    }

    return result;
  }
}

#endif
