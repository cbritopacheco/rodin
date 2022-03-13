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

   /**
    * Enumeration class to indicate whether a derived instance of
    * ShapeFunctionBase is either a Trial or Test space.
    */
   enum ShapeFunctionSpaceType
   {
      Trial, ///< Trial function space
      Test ///< Test function space
   };

   template <class FES, ShapeFunctionSpaceType Space>
   class ShapeFunction;

   template <class FES>
   class TrialFunction;

   template <class FES>
   class TestFunction;

   template <class T>
   class Component;

   /**
    * Enumeration class to indicate whether the integration should be done
    * either inside the Domain or on the Boundary.
    */
   enum IntegratorRegion
   {
      Domain, ///< Perform the integration over the interior domain
      Boundary ///< Perform the integration over the boundary of the domain
   };

   /**
    * @brief Integral class
    */
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
   template <class T>
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
   template <class T>
   class Transpose;

   template <class T>
   class Grad;

   template <class T>
   class Div;

   template <class T>
   class Jacobian;

   class RestrictionBase;

   template <class T>
   class Restriction;

   template <class Operand>
   class UnaryMinus;

   template <class Lhs, class Rhs>
   class Sum;

   template <class Lhs, class Rhs>
   class Mult;

   template <class Lhs, class Rhs>
   class Dot;
}

#endif
