/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_H
#define RODIN_ADAPTATION_H

/**
 * @file
 * @brief Top level include for Rodin adaptation utilities.
 *
 * Provides the LSR (level-set registration) solver, its
 * form-language integrators, the intrinsic shape quality energy, and the
 * Jacobian admissibility / nondegeneracy / singularity barrier. The current
 * implementation is the 2D triangular affine prototype; see the individual
 * headers for the FES-independence boundary.
 */

#include "Adaptation/AnalyticFunctionAdapters.h"
#include "Adaptation/CellGeomCache.h"
#include "Adaptation/LSR.h"
#include "Adaptation/LSRAdmissibility.h"
#include "Adaptation/LSRIntegrators.h"
#include "Adaptation/LSRParameters.h"
#include "Adaptation/LSRRegistration.h"
#include "Adaptation/LSRReport.h"
#include "Adaptation/LevelSetRegistration.h"
#include "Adaptation/JacobianAdmissibilityBarrier.h"

#endif
