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
   namespace Assembly
   {
      enum class Type;

      struct Common;
      struct Device;
   }

   // ---- Problem -----------------------------------------------------------
   class ProblemBase;

   template <class TrialFES, class TestFES, class OperatorType>
   class Problem;

   class ProblemBody;

   // ---- FiniteElementCollection -------------------------------------------
   class FiniteElementCollectionBase;

   class L2;


   // ---- FiniteElementSpace ------------------------------------------------
   class FiniteElementSpaceBase;

   // template <class FEC, class Trait = Traits::Serial>
   // class FiniteElementSpace;

   template <class Trait>
   class H1;

   // ---- GridFunction ------------------------------------------------------
   class GridFunctionBase;

   template <class FES>
   class GridFunction;

   // ---- ShapeFunction -----------------------------------------------------

   /**
    * Enumeration class to indicate whether a derived instance of
    * ShapeFunctionBase is either a Trial or Test space.
    */
   enum class ShapeFunctionSpaceType
   {
      Trial, ///< Trial function space
      Test ///< Test function space
   };

   static constexpr auto TrialSpace = ShapeFunctionSpaceType::Trial;

   static constexpr auto TestSpace  = ShapeFunctionSpaceType::Test;

   /**
    * Enumeration class to indicate whether the integration should be done
    * either inside the Domain or on the Boundary.
    */
   enum class IntegratorRegion
   {
      Domain, ///< Perform the integration over the interior domain
      Boundary ///< Perform the integration over the boundary of the domain
   };

   template <ShapeFunctionSpaceType Space>
   class ShapeFunctionBase;

   template <class FES, ShapeFunctionSpaceType Space>
   class ShapeFunction;

   template <class FES>
   class TrialFunction;

   template <class FES>
   class TestFunction;

   template <class T>
   class Component;


   // ---- LinearForm --------------------------------------------------------
   class LinearFormBase;

   template <class FES>
   class LinearForm;

   class LinearFormIntegratorBase;

   class LinearFormDomainIntegrator;

   class LinearFormBoundaryIntegrator;

   // ---- BilinearForm ------------------------------------------------------
   class BilinearFormBase;

   /**
    * @brief Represents a serial bilinear form supported on two finite element
    * spaces originating from two instances of FiniteElementCollection.
    * @tparam TrialFES Trial FiniteElementCollection
    * @tparam TestFES Test FiniteElementCollection
    * @tparam Trait Indicates if the BilinearForm is in a parallel context. It is
    * one of Traits::Serial or Traits::Parallel.
    */
   template <class TrialFES, class TestFES>
   class BilinearForm;

   class BilinearFormIntegratorBase;

   class BilinearFormDomainIntegrator;

   // ---- Boundary Conditions -----------------------------------------------
   template <class T>
   class DirichletBC;

   // ---- Coefficients ------------------------------------------------------
   class FunctionBase;

   class Function;

   class ScalarFunctionBase;

   class RangeShape;

   enum class RangeType
   {
      Scalar,
      Vector,
      Matrix
   };

   template <class ... Values>
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

   template <class T>
   class Trace;

   template <class Lhs, class Rhs>
   class Composition;

   class BasisOperator;

   /**
    * @brief Integral class
    */
   template <class Integrand>
   class Integral;

   template <class Lhs, class Rhs>
   class Integral<Dot<Grad<Lhs>, Grad<Rhs>>>;

   template <class Integrand>
   class BoundaryIntegral;

   // TODO: Refactor or remove these two classes!!!!
   class LinearFormIntegratorSum;
   class BilinearFormIntegratorSum;
}

#endif
