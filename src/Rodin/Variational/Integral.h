#ifndef RODIN_VARIATIONAL_INTEGRAL_H
#define RODIN_VARIATIONAL_INTEGRAL_H

#include <mfem.hpp>

#include "ForwardDecls.h"

#include "LinearForm.h"
#include "BilinearForm.h"

namespace Rodin::Variational
{
   template <IntegratorRegion Region, class T>
   class Integral;
}

#endif
