/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_QF_H
#define RODIN_QF_H

/**
 * @file
 * @brief Top level include for the Rodin::QF module.
 *
 * The QF (Quadrature Formula) module provides numerical integration rules
 * for computing integrals over reference polytopes. These quadrature rules
 * approximate integrals by weighted sums:
 * @f[
 *   \int_K f(x) \, dx \approx \sum_{i=1}^n w_i f(x_i)
 * @f]
 * where @f$ K @f$ is a reference polytope, @f$ x_i @f$ are quadrature points,
 * and @f$ w_i @f$ are associated weights.
 *
 * @see Rodin::QF::QuadratureFormulaBase
 * @see Rodin::QF::GaussLegendre
 * @see Rodin::QF::GrundmannMoller
 * @see Rodin::QF::Centroid
 * @see Rodin::QF::PolytopeQuadratureFormula
 * @see Rodin::QF::GaussLobatto (header-only, include QF/GaussLobato.h separately)
 */

#include "QF/QuadratureFormula.h"
#include "QF/Centroid.h"
#include "QF/GaussLegendre.h"
#include "QF/GaussLobatto.h"
#include "QF/GrundmannMoller.h"
#include "QF/PolytopeQuadratureFormula.h"

#endif
