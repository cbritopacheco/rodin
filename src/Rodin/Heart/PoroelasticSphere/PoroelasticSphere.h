/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file PoroelasticSphere.h
 * @brief User-facing alias for the 0D poroelastic sphere stepper.
 */
#ifndef RODIN_HEART_POROELASTICSPHERE_POROELASTICSPHERE_H
#define RODIN_HEART_POROELASTICSPHERE_POROELASTICSPHERE_H

#include "Rodin/Heart/CCMLC2014/HolzapfelReducedLaw.h"
#include "Rodin/Heart/PoroelasticSphere/PassiveLaw.h"
#include "Rodin/Heart/PoroelasticSphere/Solver/Stepper.h"

// The user-facing PoroelasticSphereT alias (with its default law arguments)
// is declared exactly once, in Rodin/Heart/ForwardDecls.h — an alias template
// cannot be redeclared, and default template arguments may appear in only
// one declaration.
#include "Rodin/Heart/ForwardDecls.h"

#endif
