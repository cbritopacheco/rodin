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
   class RangeShape;

   enum class RangeType
   {
      Scalar,
      Vector,
      Matrix
   };

   /**
    * Enumeration class to indicate whether the integration should be done
    * either inside the Domain or on the Boundary.
    */
   enum class IntegratorRegion
   {
      Domain, ///< Perform the integration over the interior domain
      Boundary ///< Perform the integration over the boundary of the domain
   };

   namespace Linear::Assembly
   {
      enum class Type;

      struct Common;
      struct Device;
   }

   namespace Bilinear::Assembly
   {
      enum class Type;

      struct Common;
   }

   class LinearFormBase;

   template <class FES>
   class LinearForm;

   class LinearFormIntegratorBase;

   class LinearFormDomainIntegrator;

   class LinearFormBoundaryIntegrator;

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

   class FiniteElementCollectionBase;

   class FiniteElementSpaceBase;

   class L2;

   template <class Trait>
   class H1;

   class GridFunctionBase;

   /**
    * @brief TODO
    */
   template <class FES>
   class GridFunction;

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

   class BasisOperator;

   template <ShapeFunctionSpaceType Space>
   class ShapeFunctionBase;

   template <class FES, ShapeFunctionSpaceType Space>
   class ShapeFunction;

   template <class FES>
   class TrialFunction;

   template <class FES>
   class TestFunction;

   class FunctionBase;

   class Function;

   class ScalarFunctionBase;

   template <class ... Values>
   class ScalarFunction;

   class VectorFunctionBase;

   template <class ... Values>
   class VectorFunction;

   class MatrixFunctionBase;

   class MatrixFunction;

   class BooleanFunctionBase;

   template <class T>
   class BooleanFunction;

   template <class T>
   class Component;

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

   /**
    * @ingroup DotSpecializations
    * @brief Represents the dot product between two objects
    * @tparam Lhs Type of left hand side operand
    * @tparam Rhs Type of right hand side operand
    *
    * Represents the following mathematical expression:
    * @f[
    *    \text{LHS} : \text{RHS}
    * @f]
    * where @f$ : @f$ denotes the dot product.
    *
    * For an overview of all the possible specializations of the Dot
    * class, please see @ref DotSpecializations.
    *
    * @see DotSpecializations
    */
   template <class Lhs, class Rhs>
   class Dot;

   template <class T>
   class Trace;

   template <class Lhs, class Rhs>
   class Composition;

   template <class Lhs, class Rhs>
   class LT;

   template <class Lhs, class Rhs>
   class GT;

   template <class Lhs, class Rhs>
   class EQ;

   template <class Lhs, class Rhs>
   class LEQ;

   template <class Lhs, class Rhs>
   class GEQ;

   template <class Lhs, class Rhs>
   class NEQ;

   template <class Lhs, class Rhs>
   class AND;

   template <class Lhs, class Rhs>
   class OR;

   /**
    * @ingroup IntegralSpecializations
    * @class Integral
    * @brief Represents expressions of the integral operator.
    * @tparam Integrand Type of the integrand
    *
    * Represents the integral operator with a templated integrand type:
    * @f[
    *    \int \text{Integrand} \ .
    * @f]
    *
    * For an overview of all the possible specializations of the Integral
    * class, please see @ref IntegralSpecializations.
    *
    * @see IntegralSpecializations
    */
   template <class Integrand>
   class Integral;

   template <class Integrand>
   class BoundaryIntegral;

   template <class T>
   class DirichletBC;

   class ProblemBase;

   template <class TrialFES, class TestFES, class OperatorType>
   class Problem;

   class ProblemBody;

   // TODO: Refactor or remove these two classes!!!!
   class LinearFormIntegratorSum;
   class BilinearFormIntegratorSum;
}

#endif
