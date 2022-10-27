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
   /**
    * @brief Represents the shape (dimensions) of a function.
    */
   class RangeShape;

   /**
    * @brief Represents the type of the range of a function.
    */
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

   /**
    * @brief Namespace containing utilities for the assembly of LinearForm
    * objects.
    */
   namespace Linear::Assembly
   {
      /**
       * @brief Type of assembly to be performed for LinearForm objects.
       */
      enum class Type
      {
         Common, ///< Enumerator corresponding to Linear::Assembly::Common
         Device ///< Enumerator corresponding to Linear::Assembly::Device
      };

      /**
       * @brief Struct containg the necessary data to perform a common
       * assembly.
       */
      struct Common;

      /**
       * @brief Struct containg the necessary data to perform assembly on a
       * device.
       */
      struct Device;
   }

   /**
    * @brief Namespace containing utilities for the assembly of BilinearForm
    * objects.
    */
   namespace Bilinear::Assembly
   {
      /**
       * @brief Type of assembly to be performed for BilinearForm objects.
       */
      enum class Type;

      /**
       * @brief Struct containg the necessary data to perform a common
       * assembly.
       */
      struct Common;
   }

   /**
    * @brief Base class for linear form objects.
    * @tparam VectorType Type of vector which will be assembled
    */
   template <class VectorType>
   class LinearFormBase;

   /**
    * @brief Represents a linear form on some finite element space.
    * @tparam FES Type of finite element space
    * @tparam VectorType Type of vector which will be assembled
    *
    * Represents a linear form @f$ \ell : V_h \rightarrow \mathbb{R} @f$ on a given
    * finite element space @f$ V_h @f$.
    */
   template <class FES, class VectorType>
   class LinearForm;

   /**
    * @brief Base class for linear form integrators.
    *
    * An instance of LinearFormIntegratorBase performs the assembly of the
    * element vector for each finite element.
    */
   class LinearFormIntegratorBase;

   /**
    * @brief Represents linear form integrators over a sub-domain in the mesh.
    */
   class LinearFormDomainIntegrator;

   /**
    * @brief Represents linear form integrators over the boundary of the mesh.
    */
   class LinearFormBoundaryIntegrator;

   /**
    * @brief Base class for bilinear form objects.
    * @tparam OperatorType Type of operator which will be assembled
    *
    * Represents a bilinear form @f$ a : V_h \times U_h \rightarrow \mathbb{R}
    * @f$ on given finite element spaces @f$ V_h @f$ and @f$ U_h @f$.
    */
   template <class OperatorType>
   class BilinearFormBase;

   /**
    * @brief Represents a serial bilinear form supported on two finite element
    * spaces originating from two instances of FiniteElementCollection.
    * @tparam TrialFES Trial FiniteElementCollection
    * @tparam TestFES Test FiniteElementCollection
    */
   template <class TrialFES, class TestFES, class OperatorType>
   class BilinearForm;

   /**
    * @brief Base class for bilinear form integrators.
    */
   class BilinearFormIntegratorBase;

   /**
    * @brief Represents bilinear form integrators over a sub-domain in the mesh.
    */
   class BilinearFormDomainIntegrator;

   class FiniteElementCollectionBase;

   /**
    * @brief Base class for finite element spaces.
    */
   class FiniteElementSpaceBase;

   /**
    * @brief Arbitrary order @f$ L^2(\Omega)^d @f$ conforming (continuous) finite
    * element space.
    * @tparam Type of context for the finite element space
    */
   template <class Context>
   class L2;

   /**
    * @brief Arbitrary order @f$ H^1(\Omega)^d @f$ conforming (continuous) finite
    * element space.
    * @tparam Type of context for the finite element space
    *
    * Given some discretization @f$ \mathcal{T}_h @f$ (e.g. a triangulation)
    * of @f$ \Omega @f$, instances of this class will represent the finite
    * element space
    * @f[
    *    V_h := \left\{ v : \overline{\Omega} \rightarrow \mathbb{R}^d \mid
    *       v_{|\tau} \in \mathcal{P}_\tau,
    *    \ \forall \tau \in \mathcal{T}_h \right\}
    * @f]
    * where @f$ \mathcal{P}_\tau \subset H^1(\tau) @f$ and @f$ V_h \subset
    * C^0(\Omega) @f$ so that @f$ V_h \subset H^1(\Omega)^d @f$, i.e. the
    * elements are @f$ H^1 @f$ conforming. The space @f$ P_\tau @f$ depends on
    * the kind of basis chosen.
    *
    */
   template <class Context>
   class H1;

   /**
    * @brief Base class for grid function objects.
    */
   class GridFunctionBase;

   /**
    * @brief Represents a grid function belonging to some finite element space.
    * @tparam FES Type of finite element space
    *
    * Represents a function @f$ u \in \text{FES} @f$ where FES is some discrete
    * finite element space.
    *
    * @note For an overview of all the possible specializations of the
    * GridFunction class, please see @ref GridFunctionSpecializations.
    *
    * @see GridFunctionSpecializations
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
      Test ///< %Test function space
   };

   /**
    * Shorthand variable for ShapeFunctionSpaceType::Trial.
    */
   static constexpr auto TrialSpace = ShapeFunctionSpaceType::Trial;

   /**
    * Shorthand variable for ShapeFunctionSpaceType::Test.
    */
   static constexpr auto TestSpace  = ShapeFunctionSpaceType::Test;

   class BasisOperator;

   /**
    * @brief Base class for shape function objects.
    * @tparam Space Type of shape function space (Trial or Test)
    */
   template <ShapeFunctionSpaceType Space>
   class ShapeFunctionBase;

   /**
    * @brief Base class for shape function objects.
    * @tparam FES Type of finite element space
    * @tparam Space Type of shape function space (@ref
    * ShapeFunctionSpaceType::Trial "Trial" or @ref
    * ShapeFunctionSpaceType::Test "Test")
    *
    * @note For an overview of all the possible specializations of the
    * ShapeFunction class, please see @ref ShapeFunctionSpecializations.
    *
    * @see ShapeFunctionSpecializations
    */
   template <class FES, ShapeFunctionSpaceType Space>
   class ShapeFunction;

   /**
    * @brief Represents a function which belongs to a trial space
    * @tparam FES Type of finite element space
    */
   template <class FES>
   class TrialFunction;

   /**
    * @brief Represents a function which belongs to a test space
    * @tparam FES Type of finite element space
    */
   template <class FES>
   class TestFunction;

   /**
    * @brief Base class for function objects which can be evaluated over a
    * mesh.
    *
    * Instances of FunctionBase will always have the getValue() method defined,
    * which enables the evaluation of any function on some mesh element.
    */
   class FunctionBase;

   class Function;

   /**
    * @brief Abstract base class for objects representing scalar functions.
    */
   class ScalarFunctionBase;

   /**
    * @note For an overview of all the possible specializations of the
    * ScalarFunction class, please see @ref ScalarFunctionSpecializations.
    *
    * @see ScalarFunctionSpecializations
    */
   template <class ... Values>
   class ScalarFunction;

   /**
    * @brief Abstract base class for objects representing vector functions.
    *
    * @note Vectors are zero indexed. This means that the 0-index corresponds
    * to the 1st entry of the vector.
    */
   class VectorFunctionBase;

   /**
    * @note For an overview of all the possible specializations of the
    * VectorFunction class, please see @ref VectorFunctionSpecializations.
    *
    * @see VectorFunctionSpecializations
    */
   template <class ... Values>
   class VectorFunction;

   /**
    * @brief Base class for objects representing matrix functions.
    */
   class MatrixFunctionBase;

   /**
    * @note For an overview of all the possible specializations of the
    * MatrixFunction class, please see @ref MatrixFunctionSpecializations.
    *
    * @see MatrixFunctionSpecializations
    */
   class MatrixFunction;

   /**
    * @brief Base class for objects representing boolean functions.
    */
   class BooleanFunctionBase;

   /**
    * @note For an overview of all the possible specializations of the
    * BooleanFunction class, please see @ref BooleanFunctionSpecializations.
    *
    * @see BooleanFunctionSpecializations
    */
   template <class T>
   class BooleanFunction;

   /**
    * @brief Represents the power function.
    * @tparam Base Type of the base expression
    * @tparam Exponent Type of the exponent expression
    *
    * Represents the mathematical expression:
    * @f[
    *    \text{Base}^\text{Exponent}
    * @f]
    * where Base is a type representing the base @f$ b @f$, the Exponent type
    * represents the exponent @f$ p @f$, and the exponentiation value is @f$
    * b^p @f$.
    *
    * @note For an overview of all the possible specializations of the
    * Pow class, please see @ref PowSpecializations.
    *
    * @see PowSpecializations
    */
   template <class Base, class Exponent>
   class Pow;

   template <class T>
   class Component;

   /**
    * @brief Represents the transpose matrix @f$ A^T @f$ of some matrix @f$ A
    * @f$.
    * @tparam Operand Type of operand
    *
    * Represents the mathematical expression:
    * @f[
    *    \text{Operand}^T
    * @f]
    * where Operand is a type representing
    * an @f$ n \times m @f$ matrix @f$ A @f$ and the transpose matrix @f$
    * A^T @f$ is an @f$ m \times n @f$ matrix defined by
    * @f[
    *    {A^T}_{ij} = A_{ji} \ .
    * @f]
    *
    * @note For an overview of all the possible specializations of the
    * Transpose class, please see @ref TransposeSpecializations.
    *
    * @see TransposeSpecializations
    */
   template <class Operand>
   class Transpose;

   /**
    * @brief Represents the gradient @f$ \nabla u @f$ of a scalar function
    * @f$ u @f$.
    * @tparam Operand Type of operand
    *
    * Represents the mathematical expression:
    * @f[
    *    \nabla \text{Operand}
    * @f]
    * where Operand is a type representing a scalar function
    * @f$ u : \mathbb{R}^n \rightarrow \mathbb{R} @f$ and the gradient
    * @f$ \nabla u : \mathbb{R}^n \rightarrow \mathbb{R} @f$ at the point
    * @f$ x = (x_1, \ldots, x_n) @f$ is defined by:
    * @f[
    *    \nabla u (x) =
    *    \left[
    *       \dfrac{\partial u}{\partial x_1}(x), \ldots,
    *       \dfrac{\partial u}{\partial x_n}(x)
    *    \right]^T
    * @f]
    *
    * @note For an overview of all the possible specializations of the
    * Grad class, please see @ref GradSpecializations.
    *
    * @see GradSpecializations
    */
   template <class Operand>
   class Grad;

   /**
    * @brief Represents the divergence of a vector valued function.
    * @tparam Operand Type of operand
    *
    * @note For an overview of all the possible specializations of the
    * Div class, please see @ref DivSpecializations.
    *
    * @see DivSpecializations
    */
   template <class Operand>
   class Div;

   /**
    * @brief Represents the Jacobian matrix of a type
    * @tparam Operand Type of operand
    *
    * Represents the following mathematical expression:
    * @f[
    *    \mathbf{J}_\mathrm{Operand}
    * @f]
    * where Operand is a type representing a function @f$ u : \mathbb{R}^s
    * \rightarrow \mathbb{R}^d @f$ whose Jacobian matrix @f$ \mathbf{J}_u(x)
    * @f$ at any point @f$ x = (x_1, \ldots, x_s) @f$ is defined by the @f$ s
    * \times d @f$ matrix:
    * @f[
    * \mathbf{J}_u = \begin{bmatrix}
    * \dfrac{\partial u_1}{\partial x_1} & \ldots & \dfrac{\partial u_d}{\partial x_1}\\
    * \vdots & \ddots & \vdots\\
    * \dfrac{\partial u_1}{\partial x_s} & \ldots & \dfrac{\partial u_d}{\partial x_s}
    * \end{bmatrix} .
    * @f]
    *
    * @note For an overview of all the possible specializations of the Jacobian
    * class, please see @ref JacobianSpecializations.
    *
    * @see JacobianSpecializations
    */
   template <class Operand>
   class Jacobian;

   class RestrictionBase;

   template <class Operand>
   class Restriction;

   /**
    * @brief Represent the negation of an operand.
    * @tparam Operand Type of operand
    *
    * @note For an overview of all the possible specializations of the
    * UnaryMinus class, please see @ref UnaryMinusSpecializations.
    *
    * Represents the following mathematical expression:
    * @f[
    *    - \text{Operand}
    * @f]
    *
    * @see UnaryMinusSpecializations
    */
   template <class Operand>
   class UnaryMinus;

   /**
    * @brief Represents the sum operation.
    * @tparam LHS Type of left hand side operand
    * @tparam RHS Type of right hand side operand
    *
    * @note For an overview of all the possible specializations of the
    * Sum class, please see @ref SumSpecializations.
    *
    * Represents the following mathematical expression:
    * @f[
    * \text{LHS} + \text{RHS}
    * @f]
    *
    * @see SumSpecializations
    */
   template <class LHS, class RHS>
   class Sum;

   /**
    * @brief Represents the multiplication operation.
    * @tparam LHS Type of left hand side operand
    * @tparam RHS Type of right hand side operand
    *
    * @note For an overview of all the possible specializations of the
    * Mult class, please see @ref MultSpecializations.
    *
    * Represents the following mathematical expression:
    * @f[
    * \text{LHS} * \text{RHS}
    * @f]
    *
    * @see MultSpecializations
    */
   template <class LHS, class RHS>
   class Mult;

   /**
    * @brief Represents the division operation.
    * @tparam LHS Type of left hand side operand
    * @tparam RHS Type of right hand side operand
    *
    * @note For an overview of all the possible specializations of the
    * Division class, please see @ref DivisionSpecializations.
    *
    * Represents the following mathematical expression:
    * @f[
    * \text{LHS} \div \text{RHS}
    * @f]
    *
    * @see DivisionSpecializations
    */
   template <class LHS, class RHS>
   class Division;

   /**
    * @brief Represents the dot product between two objects
    * @tparam LHS Type of left hand side operand
    * @tparam RHS Type of right hand side operand
    *
    * Represents the following mathematical expression:
    * @f[
    *    \mathrm{LHS} : \mathrm{RHS}
    * @f]
    *
    * @note For an overview of all the possible specializations of the Dot
    * class, please see @ref DotSpecializations.
    *
    * @see DotSpecializations
    */
   template <class LHS, class RHS>
   class Dot;

   /**
    * @brief Represents the trace of a matrix function
    * @tparam Operand Type of operand
    *
    * Represents the trace of a matrix valued operand:
    * @f[
    *    \mathrm{tr} \left( \mathrm{Operand} \right)
    * @f]
    * where Operand represents a matrix valued function @f$ A : \Omega
    * \rightarrow \mathbb{R}^{n \times n} @f$ and its trace @f$ \mathrm{tr} @f$
    * is defined by:
    * @f[
    *    \mathrm{tr}(A) = \sum_{i = 1}^n A_{ii} \ .
    * @f]
    *
    * @note For an overview of all the possible specializations of the Dot
    * class, please see @ref TraceSpecializations.
    *
    * @see TraceSpecializations
    */
   template <class Operand>
   class Trace;

   /**
    * @tparam LHS Type of left hand side operand
    * @tparam RHS Type of right hand side operand
    *
    * Represents the mathematical expression:
    * @f[
    * \text{LHS} \circ \text{RHS}
    * @f]
    * where LHS and RHS are types which represent, respectively, functions @f$
    * f : B \rightarrow C @f$ and @f$ g : A \rightarrow B @f$. Then their
    * composition at each point @f$ x @f$ is defined by:
    * @f[
    *    (f \circ g)(x) := f(g(x)) \ .
    * @f]
    *
    * @note For an overview of all the possible specializations of the
    * Composition class, please see @ref CompositionSpecializations.
    *
    * @see CompositionSpecializations
    */
   template <class LHS, class RHS>
   class Composition;

   /**
    * @tparam LHS Type of left hand side operand
    * @tparam RHS Type of right hand side operand
    *
    * Represents the logical LT operation between two templated types:
    * @f[
    *    \mathrm{LHS} < \mathrm{RHS} \ .
    * @f]
    *
    * @note For an overview of all the possible specializations of the
    * LT class, please see @ref LTSpecializations.
    *
    * @see LTSpecializations
    */
   template <class LHS, class RHS>
   class LT;

   /**
    * @tparam LHS Type of left hand side operand
    * @tparam RHS Type of right hand side operand
    *
    * Represents the logical GT operation between two templated types:
    * @f[
    *    \mathrm{LHS} > \mathrm{RHS} \ .
    * @f]
    *
    * @note For an overview of all the possible specializations of the
    * GT class, please see @ref GTSpecializations.
    *
    * @see GTSpecializations
    */
   template <class LHS, class RHS>
   class GT;

   /**
    * @tparam LHS Type of left hand side operand
    * @tparam RHS Type of right hand side operand
    *
    * Represents the logical EQ operation between two templated types:
    * @f[
    *    \mathrm{LHS} == \mathrm{RHS} \ .
    * @f]
    *
    * @note For an overview of all the possible specializations of the
    * EQ class, please see @ref EQSpecializations.
    *
    * @see EQSpecializations
    */
   template <class LHS, class RHS>
   class EQ;

   /**
    * @tparam LHS Type of left hand side operand
    * @tparam RHS Type of right hand side operand
    *
    * Represents the logical LEQ operation between two templated types:
    * @f[
    *    \mathrm{LHS} <= \mathrm{RHS} \ .
    * @f]
    *
    * @note For an overview of all the possible specializations of the
    * LEQ class, please see @ref LEQSpecializations.
    *
    * @see LEQSpecializations
    */
   template <class LHS, class RHS>
   class LEQ;

   /**
    * @tparam LHS Type of left hand side operand
    * @tparam RHS Type of right hand side operand
    *
    * Represents the logical GEQ operation between two templated types:
    * @f[
    *    \mathrm{LHS} >= \mathrm{RHS} \ .
    * @f]
    *
    * @note For an overview of all the possible specializations of the
    * GEQ class, please see @ref GEQSpecializations.
    *
    * @see GEQSpecializations
    */
   template <class LHS, class RHS>
   class GEQ;

   /**
    * @tparam LHS Type of left hand side operand
    * @tparam RHS Type of right hand side operand
    *
    * Represents the logical NEQ operation between two templated types:
    * @f[
    *    \mathrm{LHS} != \mathrm{RHS} \ .
    * @f]
    *
    * @note For an overview of all the possible specializations of the
    * NEQ class, please see @ref NEQSpecializations.
    *
    * @see NEQSpecializations
    */
   template <class LHS, class RHS>
   class NEQ;

   /**
    * @brief Represents the logical AND expression
    * @tparam LHS Type of left hand side operand
    * @tparam RHS Type of right hand side operand
    *
    * Represents the logical AND operation between two templated types:
    * @f[
    *    \mathrm{LHS} \land \mathrm{RHS} \ .
    * @f]
    *
    * @note For an overview of all the possible specializations of the AND class,
    * please see @ref ANDSpecializations.
    *
    * @see ANDSpecializations
    */
   template <class LHS, class RHS>
   class AND;

   /**
    * @tparam LHS Type of left hand side operand
    * @tparam RHS Type of right hand side operand
    *
    * Represents the logical OR operation between two templated types:
    * @f[
    *    \mathrm{LHS} \lor \mathrm{RHS} \ .
    * @f]
    *
    * @note For an overview of all the possible specializations of the
    * OR class, please see @ref ORSpecializations.
    *
    * @see ORSpecializations
    */
   template <class LHS, class RHS>
   class OR;

   /**
    * @brief Represents mathematical expressions of the integral operator on a
    * domain.
    * @tparam Integrand Type of the integrand
    *
    * Represents the integral operator with a templated integrand type:
    * @f[
    *    \int \mathrm{Integrand}
    * @f]
    * on a volumetric domain.
    *
    * @note For an overview of all the possible specializations of the Integral
    * class, please see @ref IntegralSpecializations.
    *
    * @see IntegralSpecializations
    */
   template <class Integrand>
   class Integral;

   /**
    * @brief Represents expressions of the integral operator on the boundary of
    * a domain.
    * @tparam Integrand Type of the integrand
    *
    * Represents the integral operator with a templated integrand type:
    * @f[
    *    \int \mathrm{Integrand}
    * @f]
    * on the boundary of a volumetric domain.
    *
    * @note For an overview of all the possible specializations of the Integral
    * class, please see @ref BoundaryIntegralSpecializations.
    *
    * @see BoundaryIntegralSpecializations
    */
   template <class Integrand>
   class BoundaryIntegral;

   /**
    * @tparam Operand Type of operand
    *
    * @note For an overview of all the possible specializations of the
    * DirichletBC class, please see @ref DirichletBCSpecializations.
    *
    * @see DirichletBCSpecializations
    */
   template <class Operand>
   class DirichletBC;

   /**
    * @brief Base class for variational problem objects.
    */
   class ProblemBase;

   /**
    * @brief Represents a variational problem.
    * @tparam Parameters Parameter pack of parameters for problem construction
    *
    * @note For an overview of all the possible specializations of the Problem
    * class, please see @ref ProblemSpecializations.
    *
    * @see ProblemSpecializations
    */
   template <class ... Parameters>
   class Problem;

   class ProblemBody;

   // TODO: Refactor or remove these two classes!!!!
   class LinearFormIntegratorSum;
   class BilinearFormIntegratorSum;
}

#endif
