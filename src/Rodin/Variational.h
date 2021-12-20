/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_H
#define RODIN_VARIATIONAL_H

/**
 * @file
 * @brief Top level include for the Rodin::Variational namespace.
 */

#include "Variational/ForwardDecls.h"

// ---- Finite elements ------------------------------------------------------
#include "Variational/FiniteElementSpace.h"
#include "Variational/H1.h"
#include "Variational/GridFunction.h"

#include "Variational/Jacobian.h"
#include "Variational/Gradient.h"
#include "Variational/Derivative.h"
#include "Variational/Transpose.h"
#include "Variational/Trace.h"
#include "Variational/Dot.h"
#include "Variational/IdentityMatrix.h"

// ---- FormLanguage ---------------------------------------------------------
#include "Variational/FormLanguage.h"
#include "Variational/Problem.h"
#include "Variational/LinearForm.h"
#include "Variational/BilinearForm.h"
#include "Variational/ScalarCoefficient.h"
#include "Variational/VectorCoefficient.h"
#include "Variational/MatrixCoefficient.h"

// ---- Bilinear form integrators --------------------------------------------
#include "Variational/DiffusionIntegrator.h"
#include "Variational/ElasticityIntegrator.h"

// ---- Linear form integrators ----------------------------------------------
#include "Variational/DomainLFIntegrator.h"

// ---- Boundary conditions --------------------------------------------------
#include "Variational/DirichletBC.h"
#include "Variational/NeumannBC.h"

#endif
