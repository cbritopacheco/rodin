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
  /**
   * @brief Represents the shape (dimensions) of a function.
   */
  class RangeShape;

  /**
   * @brief Represents the possible types of ranges for a function.
   */
  enum class RangeType;

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
  template <class FES, class Context, class VectorType>
  class LinearForm;

  class Integrator;

  /**
   * @brief Base class for linear form integrators.
   *
   * An instance of LinearFormIntegratorBase performs the assembly of the
   * element vector for each finite element.
   */
  class LinearFormIntegratorBase;

  /**
   * @brief Abstract base class for objects of type BilinearForm.
   * @tparam OperatorType Type of operator which will be assembled
   *
   * Represents a bilinear form:
   * @f[
   * a : V_h \times U_h \rightarrow \mathbb{R}
   * @f]
   * on given finite element spaces @f$ V_h @f$ and @f$ U_h @f$.
   */
  template <class OperatorType>
  class BilinearFormBase;

  /**
   * @brief Represents a serial bilinear form supported on two finite element
   * spaces originating from two instances of FiniteElementCollection.
   * @tparam TrialFES Type of trial finite element space
   * @tparam TestFES Type of test finite element space
   *
   * An object of type BilinearForm represents a bilinear map:
   * @f[
   * \begin{aligned}
   *   a : U \times V &\rightarrow \mathbb{R}\\
   *      (u, v) &\mapsto a(u, v)
   * \end{aligned}
   * @f]
   * where @f$ U @f$ and @f$ V @f$ are finite element spaces.
   *
   * @see BilinearFormSpecializations
   */
  template <class TrialFES, class TestFES, class Context, class OperatorType>
  class BilinearForm;

  /**
   * @brief Base class for bilinear form integrators.
   */
  class BilinearFormIntegratorBase;

  template <class Derived>
  class FiniteElementBase;

  template <class FES>
  class FiniteElement;

  class FiniteElementCollectionBase;

  /**
   * @brief Base class for finite element spaces.
   */
  class FiniteElementSpaceBase;

  /**
   * @brief Arbitrary order @f$ H^1(\mathcal{T}_h)^d \subset L^2 (\Omega) @f$
   * broken Sobolev space.
   * @tparam Context Type of context for the finite element space
   *
   * Given some triangulation @f$ \mathcal{T}_h @f$ of @f$ \Omega @f$,
   * instances of this class will represent the finite element space:
   * @f[
   *   H^1(\mathcal{T}_h)^d := \left\{ v \in L^2(\Omega)^d : \forall \tau \in
   *   \mathcal{T}_h, \ v |_\tau \in H^1 (\tau)^d \right\}
   * @f]
   */
  template <class Range, class Context>
  class L2;

  /**
   * @brief Arbitrary order @f$ H^1(\Omega)^d @f$ continuous finite element
   * space.
   * @tparam Context Type of context for the finite element space
   *
   * Given some triangulation @f$ \mathcal{T}_h @f$ of @f$ \Omega @f$,
   * instances of this class will represent the finite element space:
   * @f[
   *   V_h := \left\{ v : \overline{\Omega} \rightarrow \mathbb{R}^d \mid
   *   \forall \tau \in \mathcal{T}_h , \ v|_\tau \in \mathcal{P}_\tau
   *   \right\}
   * @f]
   * where @f$ \mathcal{P}_\tau \subset H^1(\tau)^d @f$ such that @f$ V_h \subset
   * C^0(\Omega)^d @f$.
   *
   */
  template <class Range, class Context>
  class H1;

  template <class Range>
  class P1Element;

  template <class Range, class Context, class ... Args>
  class P1;

  template <class StrictType>
  class LazyEvaluator;

  /**
   * @brief Base class for grid function objects.
   */
  template <class Derived, class FES>
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
   * @brief Enumeration class to indicate whether a derived instance of
   * ShapeFunctionBase belongs to either a trial or test space.
   */
  enum class ShapeFunctionSpaceType
  {
    Trial, ///< Indicates that the shape function belongs to a trial space.
    Test ///< Indicates that the shape function belongs to a test space.
  };

  /**
   * @brief Shorthand variable for ShapeFunctionSpaceType::Trial.
   */
  static constexpr auto TrialSpace = ShapeFunctionSpaceType::Trial;

  /**
   * @brief Shorthand variable for ShapeFunctionSpaceType::Test.
   */
  static constexpr auto TestSpace  = ShapeFunctionSpaceType::Test;

  template <class Value>
  class TensorBasis;

  /**
   * @brief Base class for shape function objects.
   * @tparam Space Type of shape function space (Trial or Test)
   */
  template <class Derived, class FES, ShapeFunctionSpaceType Space>
  class ShapeFunctionBase;

  /**
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
  template <class Derived, class FES, ShapeFunctionSpaceType Space>
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
  template <class Derived>
  class FunctionBase;

  class Function;

  /**
   * @brief Abstract base class for objects representing scalar functions.
   */
  template <class Derived>
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
  template <class Derived>
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
  template <class Derived>
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
  template <class Derived>
  class BooleanFunctionBase;

  /**
   * @note For an overview of all the possible specializations of the
   * BooleanFunction class, please see @ref BooleanFunctionSpecializations.
   *
   * @see BooleanFunctionSpecializations
   */
  template <class T>
  class BooleanFunction;

  template <class Operand>
  class Jump;

  template <class ... Args>
  class Average;

  /**
   * @brief Represents the power function.
   * @tparam Base Type of the base expression
   * @tparam Exponent Type of the exponent expression
   *
   * Represents the mathematical expression:
   * @f[
   *   \text{Base}^\text{Exponent}
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
   *   \text{Operand}^T
   * @f]
   * where Operand is a type representing
   * an @f$ n \times m @f$ matrix @f$ A @f$ and the transpose matrix @f$
   * A^T @f$ is an @f$ m \times n @f$ matrix defined by
   * @f[
   *   {A^T}_{ij} = A_{ji} \ .
   * @f]
   *
   * @note For an overview of all the possible specializations of the
   * Transpose class, please see @ref TransposeSpecializations.
   *
   * @see TransposeSpecializations
   */
  template <class Operand>
  class Transpose;

  template <class Operand>
  class Derivative;

  /**
   * @brief Represents the gradient @f$ \nabla u @f$ of a scalar function
   * @f$ u @f$.
   * @tparam Operand Type of operand
   *
   * Represents the mathematical expression:
   * @f[
   *   \nabla \text{Operand}
   * @f]
   * where Operand is a type representing a scalar function
   * @f$ u : \mathbb{R}^n \rightarrow \mathbb{R} @f$ and the gradient
   * @f$ \nabla u : \mathbb{R}^n \rightarrow \mathbb{R} @f$ at the point
   * @f$ x = (x_1, \ldots, x_n) @f$ is defined by:
   * @f[
   *   \nabla u (x) =
   *   \left[
   *     \dfrac{\partial u}{\partial x_1}(x), \ldots,
   *     \dfrac{\partial u}{\partial x_n}(x)
   *   \right]^T
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
   *   \mathbf{J}_\mathrm{Operand}
   * @f]
   * where Operand is a type representing a function @f$ u : \mathbb{R}^s
   * \rightarrow \mathbb{R}^d @f$ whose Jacobian matrix @f$ \mathbf{J}_u(x)
   * @f$ at any point @f$ x = (x_1, \ldots, x_s) @f$ is defined by the @f$ d
   * \times s @f$ matrix:
   * @f[
   * \mathbf{J}_u = \begin{bmatrix}
   * \dfrac{\partial u_1}{\partial x_1} & \ldots & \dfrac{\partial u_s}{\partial x_1}\\
   * \vdots & \ddots & \vdots\\
   * \dfrac{\partial u_1}{\partial x_d} & \ldots & \dfrac{\partial u_s}{\partial x_d}
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
   *   - \text{Operand}
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
   * which represents the usual arithmetic addition of two operands.
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
   * @brief Represents the dot product between two objects.
   * @tparam LHS Type of left hand side operand
   * @tparam RHS Type of right hand side operand
   *
   * Represents the following mathematical expression:
   * @f[
   *   \mathrm{LHS} : \mathrm{RHS}
   * @f]
   * known as the dot product; the algebraic operation that takes two
   * equal-length sequences of numbers (usually coordinate vectors), and
   * returns a single number.
   *
   * @note For an overview of all the possible specializations of the Dot
   * class, please see @ref DotSpecializations.
   *
   * @see DotSpecializations
   */
  template <class LHS, class RHS>
  class Dot;

  /**
   * @brief Represents the Frobenius norm.
   * @tparam Operand Type of operand
   *
   * Represents the following mathematical expression:
   * @f[
   *   || \ \text{Operand} \ ||
   * @f]
   * where Operand represents a scalar, vector or matrix @f$ A @f$, whose
   * Frobenius norm is defined by:
   * @f[
   *   || A || = \sqrt{\mathrm{tr}(A^T A)} \ .
   * @f]
   *
   * @note For an overview of all the possible specializations of the
   * Frobenius class, please see @ref FrobeniusSpecializations.
   *
   * @see FrobeniusSpecializations
   */
  template <class Operand>
  class Frobenius;

  /**
   * @brief Represents the cosine function.
   * @tparam Operand Type of operand
   *
   * Represents the following mathematical expression:
   * @f[
   *   \cos{\text{Operand}}
   * @f]
   * where Operand represents a scalar.
   */
  template <class Operand>
  class Cos;

  /**
   * @brief Represents the tangent function.
   * @tparam Operand Type of operand
   *
   * Represents the following mathematical expression:
   * @f[
   *   \tan{\text{Operand}} = \dfrac{\sin{\text{Operand}}}{\cos{\text{Operand}}}
   * @f]
   * where Operand represents a scalar.
   */
  template <class Operand>
  class Tan;

  /**
   * @brief Represents the trace of a matrix function
   * @tparam Operand Type of operand
   *
   * Represents the trace of a matrix valued operand:
   * @f[
   *   \mathrm{tr} \left( \mathrm{Operand} \right)
   * @f]
   * where Operand represents a matrix valued function @f$ A : \Omega
   * \rightarrow \mathbb{R}^{n \times n} @f$ and its trace @f$ \mathrm{tr} @f$
   * is defined by:
   * @f[
   *   \mathrm{tr}(A) = \sum_{i = 1}^n A_{ii} \ .
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
   *   (f \circ g)(x) := f(g(x)) \ .
   * @f]
   *
   * @note For an overview of all the possible specializations of the
   * Composition class, please see @ref CompositionSpecializations.
   *
   * @see CompositionSpecializations
   */
  template <class LHS, class RHS>
  class Composition;

  template <class Operand>
  class TraceOperator;

  template <class ... Args>
  class Max;

  template <class ... Args>
  class Min;

  /**
   * @tparam LHS Type of left hand side operand
   * @tparam RHS Type of right hand side operand
   *
   * Represents the logical LT operation between two templated types:
   * @f[
   *   \mathrm{LHS} < \mathrm{RHS} \ .
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
   *   \mathrm{LHS} > \mathrm{RHS} \ .
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
   *   \mathrm{LHS} == \mathrm{RHS} \ .
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
   *   \mathrm{LHS} <= \mathrm{RHS} \ .
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
   *   \mathrm{LHS} >= \mathrm{RHS} \ .
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
   *   \mathrm{LHS} != \mathrm{RHS} \ .
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
   *   \mathrm{LHS} \land \mathrm{RHS} \ .
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
   *   \mathrm{LHS} \lor \mathrm{RHS} \ .
   * @f]
   *
   * @note For an overview of all the possible specializations of the
   * OR class, please see @ref ORSpecializations.
   *
   * @see ORSpecializations
   */
  template <class LHS, class RHS>
  class OR;

  template <class Integrand>
  class QuadratureRule;

  /**
   * @brief Represents mathematical expressions of the integral operator on a
   * domain.
   * @tparam Integrand Type of the integrand
   *
   * Represents the integral operator with a templated integrand type:
   * @f[
   *   \int_{\mathcal{T}_h} \mathrm{Integrand}
   * @f]
   * on a triangulation @f$ \mathcal{T}_h @f$ of a domain @f$ \Omega @f$. The
   * domain integral is defined as the sum of all members of the
   * triangulation:
   * @f[
   *   \int_{\mathcal{T}_h} \mathrm{Integrand} := \sum_{T \in \mathcal{T}_h}
   *   \int_T \mathrm{Integrand} \ .
   * @f]
   *
   * @note For an overview of all the possible specializations of the Integral
   * class, please see @ref IntegralSpecializations.
   *
   * @see IntegralSpecializations
   */
  template <class Integrand>
  class Integral;

  template <class Integrand>
  class FaceIntegral;

  /**
   * @brief Represents expressions of the integral operator on the boundary of
   * a domain.
   * @tparam Integrand Type of the integrand
   *
   * Represents the integral operator with a templated integrand type:
   * @f[
   *   \int_{\mathcal{B}_h} \mathrm{Integrand}
   * @f]
   * on the boundary @f$ \mathcal{B}_h @f$ of a triangulation @f$
   * \mathcal{T}_h @f$ a domain @f$ \Omega @f$ with boundary @f$ \partial
   * \Omega @f$. The boundary @f$ \mathcal{B}_h @f$ is defined as the union of
   * the boundary faces:
   * @f[
   *   \mathcal{B}_h := \left\{ \partial T \cap \partial \Omega : T \in
   *   \mathcal{T}_h \right\} \ .
   * @f]
   * Then the boundary integral is just the sum of integrals over each member
   * of the boundary:
   * @f[
   *   \int_{\mathcal{B}_h} \mathrm{Integrand} := \sum_{F \in \mathcal{B}_h}
   *   \int_F \mathrm{Integrand} \ .
   * @f]
   *
   * @note For an overview of all the possible specializations of the
   * BoundaryIntegral class, please see @ref BoundaryIntegralSpecializations.
   *
   * @see BoundaryIntegralSpecializations
   */
  template <class Integrand>
  class BoundaryIntegral;

  /**
   * @brief Represents expressions of the integral operator on the interface
   * of a domain.
   * @tparam Integrand Type of the integrand
   *
   * Represents the integral operator with a templated integrand type:
   * @f[
   *   \int_{\mathcal{I}_h} \mathrm{Integrand}
   * @f]
   * on the interface @f$ \mathcal{I}_h @f$ of a triangulation @f$
   * \mathcal{T}_h @f$ of a domain @f$ \Omega @f$. The interface @f$
   * \mathcal{I}_h @f$ is defined as the union of interfaces between elements:
   * @f[
   *   \mathcal{I}_h := \left\{ \partial T_1 \cap \partial T_2 : T_1, T_2 \in
   *   \mathcal{T}_h \right\} \ .
   * @f]
   * Then then the interface integral is just the sum of integrals over each
   * member of the interface:
   * @f[
   *   \int_{\mathcal{I}_h} \mathrm{Integrand} := \sum_{F \in \mathcal{I}_h}
   *   \int_F \mathrm{Integrand} \ .
   * @f]
   *
   * @note Note that an interface @f$ F \in \mathcal{I}_h @f$ can only be a
   * face that is interior to the mesh, i.e. it is not part of the boundary.
   *
   * @note For an overview of all the possible specializations of the
   * InterfaceIntegral class, please see @ref
   * InterfaceIntegralSpecializations.
   *
   * @see InterfaceIntegralSpecializations
   */
  template <class Integrand>
  class InterfaceIntegral;

  /**
   * @tparam Operand Type of operand
   * @tparam Value Type of value
   *
   * Represents a Dirichlet boundary condition:
   * @f[
   *   \mathrm{Operand} = \mathrm{Value} \ \text{ on } \ \Gamma_D
   * @f]
   * on some subset of the boundary @f$ \Gamma_D \subset \mathcal{B}_h @f$.
   *
   * @note For an overview of all the possible specializations of the
   * DirichletBC class, please see @ref DirichletBCSpecializations.
   *
   * @see DirichletBCSpecializations
   */
  template <class Operand, class Value>
  class DirichletBC;

  class ProblemBody;

  /**
   * @brief Base class for variational problem objects.
   */
  template <class OperatorType, class VectorType>
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

  namespace Assembly
  {
    template <class Operand>
    class AssemblyBase;

    template <class Operand>
    class Native;

    template <class Operand>
    class OpenMP;
  }

  class ShapeComputator;
}

#endif
