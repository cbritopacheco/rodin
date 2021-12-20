/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FORWARDDECLS_H
#define RODIN_VARIATIONAL_FORWARDDECLS_H

namespace Rodin::Variational
{
   // ---- Problem -----------------------------------------------------------
   class ProblemBase;

   template <class FEC>
   class Problem;

   // ---- FiniteElementSpace ------------------------------------------------
   template <class Derived>
   class FiniteElementSpace;

   class H1;

   // ---- GridFunction ------------------------------------------------------
   class GridFunctionBase;

   template <class ...>
   class GridFunction;

   template <>
   class GridFunction<>;

   template <class FEC>
   class GridFunction<FEC>;

   // ---- LinearForm --------------------------------------------------------
   class LinearFormBase;

   template <class FEC>
   class LinearForm;

   class LinearFormIntegratorBase;

   class DomainLFIntegrator;

   // ---- BilinearForm ------------------------------------------------------
   class BilinearFormBase;

   template <class FEC>
   class BilinearForm;

   class BilinearFormIntegratorBase;

   class DiffusionIntegrator;

   class ElasticityIntegrator;

   // ---- Boundary Conditions -----------------------------------------------
   class BoundaryConditionBase;

   class DirichletBC;

   class NeumannBC;

   // ---- Coefficients ------------------------------------------------------
   class ScalarCoefficientBase;

   template <class T, class Enable = void>
   class ScalarCoefficient;

   class VectorCoefficientBase;

   class VectorCoefficient;

   class MatrixCoefficientBase;

   class MatrixCoefficient;

   class Jacobian;

   template <int DirectionIndex, int ComponentIndex = 1>
   class Derivative;

   class Transpose;

   template <class A, class B, class Enable = void>
   class Dot;
}

#endif
