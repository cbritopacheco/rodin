/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MODELS_DISTANCE_FMM_H
#define RODIN_MODELS_DISTANCE_FMM_H

#include "Rodin/Solver/CG.h"

#include "Rodin/Variational/Problem.h"
#include "Rodin/Variational/Integral.h"
#include "Rodin/Variational/DirichletBC.h"
#include "Rodin/Variational/GridFunction.h"

#include "Rodin/Alert/MemberFunctionException.h"

namespace Rodin::Models::Distance
{
  /**
   * @brief Poisson approximation to the distance function.
   */
  class FMM
  {
    public:
      enum class Label
      {
        Far,
        Considered,
        Accepted
      };

      template <class FES>
      auto operator()(
          Geometry::Attribute interface,
          Geometry::Attribute region,
          const FES& fes) const
      {
        return operator()(
            FlatSet<Geometry::Attribute>{interface},
            FlatSet<Geometry::Attribute>{region},
            fes);
      }

      template <class FES>
      auto operator()(
          const FlatSet<Geometry::Attribute>& interface,
          const FlatSet<Geometry::Attribute>& domain,
          const FES& fes) const
      {
        using RangeType = typename FormLanguage::Traits<FES>::RangeType;
        static_assert(std::is_same_v<RangeType, Real>);
        static_assert(std::numeric_limits<Real>::has_infinity);
        Variational::GridFunction u(fes);
        std::vector<Label>  labels(fes.getSize(), Label::Far);
        u = std::numeric_limits<Real>::infinity();
        const auto& mesh = fes.getMesh();
        const auto& meshDim = mesh.getDimension();
        FlatSet<Index> considered;
        for (Index i = 0; i < mesh.getFaceCount(); i++)
        {
          if (interface.count(mesh.getAttribute(meshDim - 1, i)))
          {
            const auto& fe = fes.getFiniteElement(meshDim - 1, i);
            for (size_t local = 0; local < fe.getCount(); local++)
            {
              const Index global = fes.getGlobalIndex({ meshDim - 1, i }, local);
              labels[global] = Label::Accepted;
              u.setValue(global, 0);
            }
          }
        }
        // TODO. Still needs to be implemented.
        assert(false);
        return u;
      }
  };
}

#endif



