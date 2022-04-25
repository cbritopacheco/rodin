/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FORWARDDECLS_H
#define RODIN_VARIATIONAL_FORWARDDECLS_H

#include "Rodin/Traits.h"

namespace Rodin::Variational
{
   // ---- Problem -----------------------------------------------------------
   class ProblemBase;

   template <class TrialFEC, class TestFEC, class OperatorType, class Trait>
   class Problem;

   class ProblemBody;

   // ---- FiniteElementCollection -------------------------------------------
   class FiniteElementCollectionBase;

   class L2;

   class H1;

   // ---- FiniteElementSpace ------------------------------------------------
   class FiniteElementSpaceBase;

   template <class FEC, class Trait = Traits::Serial>
   class FiniteElementSpace;

   // ---- GridFunction ------------------------------------------------------
   class GridFunctionBase;

   template <class FEC, class Trait = Traits::Serial>
   class GridFunction;

   /**
    * Enumeration class to indicate whether a derived instance of
    * ShapeFunctionBase is either a Trial or Test space.
    */
   enum ShapeFunctionSpaceType
   {
      Trial, ///< Trial function space
      Test ///< Test function space
   };

   template <ShapeFunctionSpaceType Space>
   class ShapeFunctionBase;

   template <class FEC, ShapeFunctionSpaceType Space>
   class ShapeFunction;

   template <class FEC, class Trait>
   class TrialFunction;

   template <class FEC, class Trait>
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

   template <class FEC, class Trait>
   class LinearForm;

   class LinearFormIntegratorBase;

   class LinearFormDomainIntegrator;

   class LinearFormBoundaryIntegrator;

   // ---- BilinearForm ------------------------------------------------------
   class BilinearFormBase;

   /**
    * @brief Represents a serial bilinear form supported on two finite element
    * spaces originating from two instances of FiniteElementCollection.
    * @tparam TrialFEC Trial FiniteElementCollection
    * @tparam TestFEC Test FiniteElementCollection
    * @tparam Trait Indicates if the BilinearForm is in a parallel context. It is
    * one of Traits::Serial or Traits::Parallel.
    */
   template <class TrialFEC, class TestFEC, class Trait>
   class BilinearForm;

   class BilinearFormIntegratorBase;

   class BilinearFormDomainIntegrator;

   class DiffusionIntegrator;

   class ElasticityIntegrator;

   // ---- Boundary Conditions -----------------------------------------------
   template <class T, class Value>
   class DirichletBC;

   // ---- Coefficients ------------------------------------------------------
   class ScalarFunctionBase;

   template <class T>
   class ScalarFunction;

   class VectorFunctionBase;

   template <class ... Values>
   class VectorFunction;

   class MatrixFunctionBase;

   class MatrixFunction;

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
   class Division;

   template <class Lhs, class Rhs>
   class Dot;
}

#endif
