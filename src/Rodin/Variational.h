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

#include "Variational/P1.h"

#include "Variational/GridFunction.h"
#include "Variational/FiniteElementSpace.h"
#include "Variational/FiniteElementCollection.h"

#include "Variational/TrialFunction.h"
#include "Variational/TestFunction.h"

#include "Variational/Component.h"
#include "Variational/Restriction.h"

#include "Variational/LinearForm.h"
#include "Variational/BilinearForm.h"

#include "Variational/Dot.h"
#include "Variational/Pow.h"
#include "Variational/Sum.h"
#include "Variational/Mult.h"
#include "Variational/Minus.h"
#include "Variational/Division.h"
#include "Variational/UnaryMinus.h"

#include "Variational/GT.h"
#include "Variational/LT.h"
#include "Variational/GEQ.h"
#include "Variational/LEQ.h"
#include "Variational/AND.h"
#include "Variational/OR.h"

#include "Variational/Div.h"
#include "Variational/Grad.h"
#include "Variational/Normal.h"
#include "Variational/FaceNormal.h"
#include "Variational/BoundaryNormal.h"
#include "Variational/Jacobian.h"
#include "Variational/Derivative.h"
#include "Variational/Jump.h"
#include "Variational/Average.h"

#include "Variational/Trace.h"
#include "Variational/Transpose.h"
#include "Variational/IdentityMatrix.h"

#include "Variational/Max.h"
#include "Variational/Min.h"

#include "Variational/Cosine.h"
#include "Variational/Tangent.h"

#include "Variational/Frobenius.h"

#include "Variational/Integral.h"
#include "Variational/FaceIntegral.h"
#include "Variational/BoundaryIntegral.h"
#include "Variational/Problem.h"

#include "Variational/ScalarFunction.h"
#include "Variational/VectorFunction.h"
#include "Variational/MatrixFunction.h"
#include "Variational/BooleanFunction.h"

#include "Variational/DirichletBC.h"

#endif
