/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file LinearForm.h
 * @brief Linear form specialization for P1 functions (uses generic implementation).
 *
 * This file includes the generic LinearForm header for P1 spaces.
 * P1 linear forms are assembled using nodal integration with the
 * appropriate quadrature rules.
 *
 * ## Mathematical Background
 * A linear form @f$ L(v) @f$ for test functions @f$ v \in \mathbb{P}_1 @f$
 * is typically of the form:
 * @f[
 *   L(v) = \int_\Omega f \cdot v \, dx + \int_{\partial\Omega} g \cdot v \, ds
 * @f]
 *
 * @see LinearForm, P1
 */
#ifndef RODIN_VARIATIONAL_P1_LINEARFORM_H
#define RODIN_VARIATIONAL_P1_LINEARFORM_H

#include "Rodin/Variational/LinearForm.h"

namespace Rodin::Variational
{
}

#endif
