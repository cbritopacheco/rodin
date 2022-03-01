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
#include "Variational/Transpose.h"
#include "Variational/Trace.h"
#include "Variational/Dot.h"
#include "Variational/IdentityMatrix.h"

#include "Variational/TrialFunction.h"
#include "Variational/TestFunction.h"

#include "Variational/Integral.h"
#include "Variational/BoundaryIntegral.h"

// ---- FormLanguage ---------------------------------------------------------
#include "Variational/FormLanguage.h"
#include "Variational/Problem.h"
#include "Variational/LinearForm.h"
#include "Variational/BilinearForm.h"
#include "Variational/Restriction.h"
#include "Variational/ScalarCoefficient.h"
#include "Variational/VectorCoefficient.h"
#include "Variational/MatrixCoefficient.h"
#include "Variational/Multiplication.h"

#include "Variational/VectorDomainLFDivIntegrator.h"

// ---- Bilinear form integrators --------------------------------------------
#include "Variational/VectorMassIntegrator.h"
#include "Variational/VectorDiffusionIntegrator.h"
#include "Variational/ElasticityIntegrator.h"
#include "Variational/DiffusionIntegrator.h"

// ---- Linear form integrators ----------------------------------------------
#include "Variational/DomainLFIntegrator.h"
#include "Variational/VectorDomainLFIntegrator.h"
#include "Variational/VectorBoundaryFluxLFIntegrator.h"

// ---- Boundary conditions --------------------------------------------------
#include "Variational/DirichletBC.h"

#endif
