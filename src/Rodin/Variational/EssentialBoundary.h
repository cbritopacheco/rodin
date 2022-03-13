#ifndef RODIN_VARIATIONAL_ESSENTIALBOUNDARY_H
#define RODIN_VARIATIONAL_ESSENTIALBOUNDARY_H

#include <map>
#include <memory>
#include <variant>

#include "Component.h"
#include "TrialFunction.h"
#include "ScalarCoefficient.h"
#include "VectorCoefficient.h"

namespace Rodin::Variational
{
   class EssentialBoundary
   {
      public:
         void add(const DirichletBC<TrialFunction<H1>>& dbc);
         void add(const DirichletBC<Component<TrialFunction<H1>>>& dbc);
   };
}

#endif
