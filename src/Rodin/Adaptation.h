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
 * Provides form-language integrators for static SDR (Surface-aware
 * Distance Reconstruction) on a P1 vector finite-element space, the
 * intrinsic shape quality energy, and the Jacobian admissibility /
 * nondegeneracy / singularity barrier. The current implementation is
 * the 2D triangular affine prototype; see the individual headers for the
 * FES-independence boundary.
 */

#include "Adaptation/AnalyticFunctionAdapters.h"
#include "Adaptation/CellGeomCache.h"
#include "Adaptation/SDRIntegrators.h"
#include "Adaptation/SignedDistanceRegistration.h"
#include "Adaptation/JacobianAdmissibilityBarrier.h"

#endif
