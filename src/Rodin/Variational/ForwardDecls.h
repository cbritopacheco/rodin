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

   template <class TrialFES, class TestFES>
   class Problem;

   // ---- FiniteElementSpace ------------------------------------------------
   class FiniteElementSpaceBase;

   class H1;

   // ---- GridFunction ------------------------------------------------------
   class GridFunctionBase;

   class GridFunctionIndexBase;

   template <class T>
   class GridFunctionIndex;

   template <class FES>
   class GridFunction;

   class IncompleteGridFunction;

   enum ShapeFunctionSpaceType
   {
      Trial,
      Test
   };

   template <class FES, ShapeFunctionSpaceType Space>
   class ShapeFunction;

   template <class FES>
   class TrialFunction;

   template <class FES>
   class TestFunction;

   /**
    * @brief Where the integration takes place.
    */
   enum IntegratorRegion
   {
      Domain,
      Boundary
   };

   template <class Integrand>
   class Integral;

   template <class Integrand>
   class BoundaryIntegral;

   // ---- LinearForm --------------------------------------------------------
   class LinearFormBase;

   template <class FES>
   class LinearForm;

   class LinearFormIntegratorBase;

   class LinearFormDomainIntegrator;

   class LinearFormBoundaryIntegrator;

   class DomainLFIntegrator;

   // ---- BilinearForm ------------------------------------------------------
   class BilinearFormBase;

   template <class FES>
   class BilinearForm;

   class BilinearFormIntegratorBase;

   class BilinearFormDomainIntegrator;

   class DiffusionIntegrator;

   class ElasticityIntegrator;

   // ---- Boundary Conditions -----------------------------------------------
   class DirichletBC;

   // ---- Coefficients ------------------------------------------------------
   class ScalarCoefficientBase;

   template <class T>
   class ScalarCoefficient;

   class VectorCoefficientBase;

   template <class ... Values>
   class VectorCoefficient;

   class MatrixCoefficientBase;

   class MatrixCoefficient;

   // ---- Expressions -------------------------------------------------------
   class Jacobian;

   class Transpose;

   template <class Lhs, class Rhs>
   class Dot;

   template <class T>
   class Gradient;

   class RestrictionBase;

   template <class T>
   class Restriction;

   namespace Internal
   {
      class ScalarCoefficient;
      class VectorCoefficient;
      class MatrixCoefficient;
   }
}

#endif
