#ifndef RODIN_VARIATIONAL_RELATIVEERROR_H
#define RODIN_VARIATIONAL_RELATIVEERROR_H

#include "GridFunction.h"

namespace Rodin::Variational
{
    class RelativeError
    {
      public:
        enum class Norm
        {
          L1,
          L2,
          LInf
        };

        template <class FES, class FunctionType>
        static Real l1(const GridFunction<FES>& model, const FunctionType& exact)
        {
          return compute(model, exact, Norm::L1);
        }

        template <class FES, class FunctionType>
        static Real l2(const GridFunction<FES>& model, const FunctionType& exact)
        {
          return compute(model, exact, Norm::L2);
        }

        template <class FES, class FunctionType>
        static Real lInf(const GridFunction<FES>& model, const FunctionType& exact)
        {
          return compute(model, exact, Norm::LInf);
        }

        template <class FES, class FunctionType>
        static Real compute(const GridFunction<FES>& model, const FunctionType& exact, const Norm& norm)
        {
            switch (norm)
            {
                case Norm::L1:
                {
                  const auto& fes = model.getFiniteElementSpace();
                  GridFunction<FES> exactNorm(fes);
                  exactNorm = Abs(exact);
                  exactNorm.setWeights();
                  GridFunction<FES> diff(fes);
                  diff = Abs(model - exact);
                  diff.setWeights();
                  return Integral(diff).compute() / Integral(exactNorm).compute();
                }
                case Norm::L2:
                {
                  const auto& fes = model.getFiniteElementSpace();
                  GridFunction<FES> exactNorm(fes);
                  exactNorm = [&](const Geometry::Point& p) { auto v = exact(p); return Math::dot(v, v); };
                  exactNorm.setWeights();
                  GridFunction<FES> diff(fes);
                  diff = [&](const Geometry::Point& p) { auto v = exact(p) - model(p); return Math::dot(v, v); };
                  diff.setWeights();
                  return sqrt(Integral(diff).compute()) / sqrt(Integral(exactNorm).compute());
                }
                case Norm::LInf:
                {
                  const auto& fes = model.getFiniteElementSpace();
                  GridFunction<FES> exactNorm(fes);
                  exactNorm = Abs(exact);
                  GridFunction<FES> diff(fes);
                  diff = Abs(model - exact);
                  return diff.max() / exactNorm.max();
                }
            }
            return Math::Constants::nan();
        }
    };
}

#endif

