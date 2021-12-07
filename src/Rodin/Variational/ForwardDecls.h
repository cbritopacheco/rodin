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

   template <class T>
   class DomainLFIntegrator;

   // ---- BilinearForm ------------------------------------------------------
   class BilinearFormBase;

   template <class FEC>
   class BilinearForm;

   template <class T>
   class DiffusionIntegrator;

   template <class L, class M>
   class ElasticityIntegrator;

   // ---- Boundary Conditions -----------------------------------------------
   template <class T>
   class DirichletBC;

   template <class T>
   class NeumannBC;

   // ---- Coefficients ------------------------------------------------------
   template <class T, class Enable = void>
   class ScalarCoefficient;

   template <class ... Values>
   class VectorCoefficient;
}

#endif
