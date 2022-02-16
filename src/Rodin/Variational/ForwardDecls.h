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
   class FiniteElementSpaceBase;

   template <class Derived>
   class FiniteElementSpace;

   class H1;

   // ---- GridFunction ------------------------------------------------------
   class GridFunctionBase;

   class GridFunctionIndexBase;

   template <class T>
   class GridFunctionIndex;

   template <class FEC>
   class GridFunction;

   class IncompleteGridFunction;

   /**
    * @brief Where the integration takes place.
    */
   enum IntegratorRegion
   {
      Domain,
      Boundary
   };

   // ---- LinearForm --------------------------------------------------------
   class LinearFormBase;

   template <class FEC>
   class LinearForm;

   class LinearFormIntegratorBase;

   class LinearFormDomainIntegrator;

   class LinearFormBoundaryIntegrator;

   class DomainLFIntegrator;

   // ---- BilinearForm ------------------------------------------------------
   class BilinearFormBase;

   template <class FEC>
   class BilinearForm;

   class BilinearFormIntegratorBase;

   class BilinearFormDomainIntegrator;

   class DiffusionIntegrator;

   class ElasticityIntegrator;

   // ---- Boundary Conditions -----------------------------------------------
   class BoundaryConditionBase;

   template <class T>
   class BoundaryCondition;

   template <class T>
   class DirichletBC;

   template <class T>
   class NeumannBC;

   // ---- Coefficients ------------------------------------------------------
   class ScalarCoefficientBase;

   template <class T>
   class ScalarCoefficient;

   class VectorCoefficientBase;

   template <class ... Values>
   class VectorCoefficient;

   class MatrixCoefficientBase;

   class MatrixCoefficient;

   class Jacobian;

   template <int DirectionIndex, int ComponentIndex = 1>
   class Derivative;

   class Transpose;

   template <class A, class B, class Enable = void>
   class Dot;

   class RestrictionBase;

   template <class T>
   class Restriction;
}

#endif
