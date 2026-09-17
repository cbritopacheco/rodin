/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_QF_H
#define RODIN_QF_H

/**
 * @file
 * @brief Top level include for the Rodin::QF module.
 *
 * The QF (Quadrature Formula) module provides numerical integration rules
 * for computing integrals over reference polytopes. These quadrature rules
 * approximate integrals by weighted sums:
 * @f[
 *   \int_K f(x) \, dx \approx \sum_{i=1}^n w_i f(x_i)
 * @f]
 * where @f$ K @f$ is a reference polytope, @f$ x_i @f$ are quadrature points,
 * and @f$ w_i @f$ are associated weights.
 *
 * @see <a href="class_rodin_1_1_q_f_1_1_quadrature_formula_base.html">Rodin::QF::QuadratureFormulaBase</a>
 * @see <a href="class_rodin_1_1_q_f_1_1_gauss_legendre.html">Rodin::QF::GaussLegendre</a>
 * @see <a href="class_rodin_1_1_q_f_1_1_grundmann_moller.html">Rodin::QF::GrundmannMoller</a>
 * @see <a href="class_rodin_1_1_q_f_1_1_centroid.html">Rodin::QF::Centroid</a>
 * @see <a href="class_rodin_1_1_q_f_1_1_polytope_quadrature_formula.html">Rodin::QF::PolytopeQuadratureFormula</a>
 * @see <a href="class_rodin_1_1_q_f_1_1_xiao_gimbutas.html">Rodin::QF::XiaoGimbutas</a>
 * @see <a href="class_rodin_1_1_q_f_1_1_witherden_vincent.html">Rodin::QF::WitherdenVincent</a>
 * @see <a href="class_rodin_1_1_q_f_1_1_tensor_product.html">Rodin::QF::TensorProduct</a>
 * @see <a href="class_rodin_1_1_q_f_1_1_gauss_lobatto.html">Rodin::QF::GaussLobatto</a> header-only, include QF/GaussLobato.h separately
 */

#include "QF/QuadratureFormula.h"
#include "QF/Centroid.h"
#include "QF/GaussLegendre.h"
#include "QF/GaussLobatto.h"
#include "QF/GrundmannMoller.h"
#include "QF/PolytopeQuadratureFormula.h"
#include "QF/TensorProduct.h"
#include "QF/WitherdenVincent.h"
#include "QF/XiaoGimbutas.h"

#endif
