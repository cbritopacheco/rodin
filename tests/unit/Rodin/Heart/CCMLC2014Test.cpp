// /*
//  *          Copyright Carlos BRITO PACHECO 2021 - 2026.
//  * Distributed under the Boost Software License, Version 1.0.
//  *       (See accompanying file LICENSE or copy at
//  *          https://www.boost.org/LICENSE_1_0.txt)
//  */
// #include <gtest/gtest.h>
// 
// #include <cmath>
// #include <type_traits>
// #include <algorithm>
// 
// #include "Rodin/Heart/CCMLC2014.h"
// 
// using namespace Rodin;
// 
// namespace
// {
//   using Model = Heart::CCMLC2014;
//   using Scalar = Real;
//   using DenseVector = Model::DenseVector;
//   using DenseMatrix = Model::DenseMatrix;
// 
//   Model::Input makeGenericInput()
//   {
//     Model::Input input;
//     input.rho = 1.2;
//     input.d0 = 0.8;
//     input.R0 = 1.1;
//     input.eta = 0.05;
//     input.Es = 2.0;
//     input.mu = 0.9;
//     input.alpha = 0.3;
//     input.k0 = 0.4;
//     input.sigma0 = 0.35;
//     input.Cp = 1.4;
//     input.Cd = 1.8;
//     input.Rp = 1.1;
//     input.Rd = 1.6;
//     input.Kat = 0.6;
//     input.Kp = 0.4;
//     input.Kar = 0.9;
//     input.absRegularization = 1e-8;
// 
//     input.u = [](Scalar t) { return 0.2 * std::sin(t); };
//     input.pAt = [](Scalar t) { return 0.1 + 0.05 * std::cos(t); };
//     input.pSv = [](Scalar) { return 0.3; };
//     input.n0 = [](Scalar ec) { return 0.9 + 0.05 * std::tanh(ec); };
// 
//     {
//       using PassiveEnergy = std::decay_t<decltype(input.passiveEnergy)>;
//       typename PassiveEnergy::Parameters hp;
//       hp.mu1 = 0.0;
//       hp.mu2 = 0.01;
//       hp.C0 = 0.05;
//       hp.C1 = 0.5;
//       hp.C2 = 0.03;
//       hp.C3 = 0.6;
//       input.passiveEnergy = PassiveEnergy(hp);
//     }
// 
//     return input;
//   }
// 
//   DenseVector makeGenericStateVector()
//   {
//     DenseVector x(Model::NVAR);
//     x[Model::Y]    = 0.08;
//     x[Model::V]    = -0.03;
//     x[Model::PV]   = 0.12;
//     x[Model::PAR]  = 0.07;
//     x[Model::PD]   = 0.03;
//     x[Model::EC]   = 0.015;
//     x[Model::KC]   = 0.02;
//     x[Model::TAUC] = 0.04;
//     return x;
//   }
// 
//   Model::State makeGenericState(Scalar t)
//   {
//     Model::State s;
//     s.y = 0.05;
//     s.v = -0.01;
//     s.pv = 0.10;
//     s.par = 0.08;
//     s.pd = 0.04;
//     s.ec = 0.012;
//     s.kc = 0.018;
//     s.tauc = 0.035;
//     s.t = t;
//     return s;
//   }
// 
//   DenseMatrix finiteDifferenceStaticJacobian(
//       const Model::Input& input,
//       const DenseVector& x,
//       const Model::State& base)
//   {
//     DenseVector R0;
//     Model::evaluateStaticResidual(input, x, base, R0);
// 
//     DenseMatrix J(R0.size(), x.size());
//     const Scalar eps = 1e-7;
// 
//     for (Eigen::Index j = 0; j < x.size(); ++j)
//     {
//       const Scalar h = eps * std::max<Scalar>(1.0, std::abs(x[j]));
//       DenseVector xp = x;
//       DenseVector xm = x;
//       xp[j] += h;
//       xm[j] -= h;
// 
//       DenseVector Rp, Rm;
//       Model::evaluateStaticResidual(input, xp, base, Rp);
//       Model::evaluateStaticResidual(input, xm, base, Rm);
// 
//       J.col(j) = (Rp - Rm) / (2.0 * h);
//     }
// 
//     return J;
//   }
// 
//   DenseMatrix finiteDifferenceDynamicJacobian(
//       const Model::Input& input,
//       const DenseVector& x,
//       const Model::State& sn,
//       const Model::State& snm1,
//       Model::DynamicScheme scheme,
//       Scalar tnp1,
//       Scalar dt)
//   {
//     DenseVector R0;
//     Model::evaluateDynamicResidual(input, x, sn, snm1, scheme, tnp1, dt, R0);
// 
//     DenseMatrix J(R0.size(), x.size());
//     const Scalar eps = 1e-7;
// 
//     for (Eigen::Index j = 0; j < x.size(); ++j)
//     {
//       const Scalar h = eps * std::max<Scalar>(1.0, std::abs(x[j]));
//       DenseVector xp = x;
//       DenseVector xm = x;
//       xp[j] += h;
//       xm[j] -= h;
// 
//       DenseVector Rp, Rm;
//       Model::evaluateDynamicResidual(input, xp, sn, snm1, scheme, tnp1, dt, Rp);
//       Model::evaluateDynamicResidual(input, xm, sn, snm1, scheme, tnp1, dt, Rm);
// 
//       J.col(j) = (Rp - Rm) / (2.0 * h);
//     }
// 
//     return J;
//   }
// 
//   Scalar relativeMatrixError(const DenseMatrix& A, const DenseMatrix& B)
//   {
//     return (A - B).norm() / std::max<Scalar>(B.norm(), 1e-14);
//   }
// }
// 
// TEST(CCMLC2014Test, StaticJacobianMatchesFiniteDifference)
// {
//   const auto input = makeGenericInput();
//   const auto base = makeGenericState(0.2);
//   auto x = makeGenericStateVector();
// 
//   DenseMatrix Janalytic;
//   Model::evaluateStaticJacobian(input, x, base, Janalytic);
// 
//   const auto Jfd = finiteDifferenceStaticJacobian(input, x, base);
// 
//   const Scalar rel = relativeMatrixError(Janalytic, Jfd);
//   EXPECT_LT(rel, 5e-5);
// }
// 
// TEST(CCMLC2014Test, DynamicJacobianBackwardEulerMatchesFiniteDifference)
// {
//   const auto input = makeGenericInput();
// 
//   auto x = makeGenericStateVector();
//   const auto sn = makeGenericState(0.2);
//   auto snm1 = makeGenericState(0.19);
// 
//   // Make snm1 slightly different from sn
//   snm1.y -= 0.01;
//   snm1.v += 0.005;
//   snm1.pv -= 0.01;
//   snm1.par += 0.002;
//   snm1.pd -= 0.003;
//   snm1.ec -= 0.002;
//   snm1.kc += 0.001;
//   snm1.tauc -= 0.004;
// 
//   const Scalar dt = 1e-2;
//   const Scalar tnp1 = sn.t + dt;
// 
//   DenseMatrix Janalytic;
//   Model::evaluateDynamicJacobian(
//       input, x, sn, snm1, Model::DynamicScheme::BackwardEuler, tnp1, dt, Janalytic);
// 
//   const auto Jfd = finiteDifferenceDynamicJacobian(
//       input, x, sn, snm1, Model::DynamicScheme::BackwardEuler, tnp1, dt);
// 
//   const Scalar rel = relativeMatrixError(Janalytic, Jfd);
//   EXPECT_LT(rel, 5e-5);
// }
// 
// TEST(CCMLC2014Test, DynamicJacobianBDF2MatchesFiniteDifference)
// {
//   const auto input = makeGenericInput();
// 
//   auto x = makeGenericStateVector();
//   const auto sn = makeGenericState(0.2);
//   auto snm1 = makeGenericState(0.19);
// 
//   snm1.y -= 0.01;
//   snm1.v += 0.005;
//   snm1.pv -= 0.01;
//   snm1.par += 0.002;
//   snm1.pd -= 0.003;
//   snm1.ec -= 0.002;
//   snm1.kc += 0.001;
//   snm1.tauc -= 0.004;
// 
//   const Scalar dt = 1e-2;
//   const Scalar tnp1 = sn.t + dt;
// 
//   DenseMatrix Janalytic;
//   Model::evaluateDynamicJacobian(
//       input, x, sn, snm1, Model::DynamicScheme::BDF2, tnp1, dt, Janalytic);
// 
//   const auto Jfd = finiteDifferenceDynamicJacobian(
//       input, x, sn, snm1, Model::DynamicScheme::BDF2, tnp1, dt);
// 
//   const Scalar rel = relativeMatrixError(Janalytic, Jfd);
//   EXPECT_LT(rel, 5e-5);
// }
// 
// TEST(CCMLC2014Test, StaticEquilibriumConvergesForMildPassiveCase)
// {
//   Model::Input input;
//   input.rho = 1.0;
//   input.d0 = 1.0;
//   input.R0 = 1.0;
//   input.eta = 0.0;
//   input.Es = 0.5;
//   input.mu = 1.0;
//   input.alpha = 0.0;
//   input.k0 = 0.0;
//   input.sigma0 = 0.0;
// 
//   input.Cp = 1.0;
//   input.Cd = 1.0;
//   input.Rp = 1.0;
//   input.Rd = 1.0;
// 
//   input.Kat = 1.0;
//   input.Kp = 1.0;
//   input.Kar = 1.0;
// 
//   input.pAt = [](Scalar) { return 0.0; };
//   input.pSv = [](Scalar) { return 0.0; };
//   input.u   = [](Scalar) { return 0.0; };
//   input.n0  = [](Scalar) { return 0.0; };
// 
//   {
//     using PassiveEnergy = std::decay_t<decltype(input.passiveEnergy)>;
//     typename PassiveEnergy::Parameters hp;
//     hp.mu1 = 0.0;
//     hp.mu2 = 0.0;
//     hp.C0 = 0.01;
//     hp.C1 = 0.2;
//     hp.C2 = 0.01;
//     hp.C3 = 0.2;
//     input.passiveEnergy = PassiveEnergy(hp);
//   }
// 
//   Model model(input);
//   model.setMaxIterations(50)
//        .setAbsoluteTolerance(1e-10)
//        .setRelativeTolerance(1e-10)
//        .setStepTolerance(1e-12)
//        .setDampingFactor(1.0);
// 
//   Model::State initial;
//   initial.t = 0.0;
//   initial.y = 0.01;
//   initial.v = 0.0;
//   initial.pv = 0.0;
//   initial.par = 0.0;
//   initial.pd = 0.0;
//   initial.ec = 0.0;
//   initial.kc = 0.0;
//   initial.tauc = 0.0;
// 
//   model.initialize(initial);
// 
//   const auto rep = model.solveStaticEquilibrium();
//   EXPECT_TRUE(rep.converged);
// 
//   const auto& s = model.getState();
//   EXPECT_NEAR(s.v, 0.0, 1e-12);
//   EXPECT_NEAR(s.kc, 0.0, 1e-12);
//   EXPECT_NEAR(s.tauc, 0.0, 1e-12);
// 
//   const Scalar e1D =
//     0.5 * (std::pow(1.0 + s.y / input.R0, 2) - 1.0);
//   EXPECT_NEAR(s.ec, e1D, 1e-10);
// }
// 
// TEST(CCMLC2014Test, StaticThenBackwardEulerThenBDF2ConvergeForMildCase)
// {
//   Model::Input input;
//   input.rho = 1.0;
//   input.d0 = 1.0;
//   input.R0 = 1.0;
//   input.eta = 0.01;
//   input.Es = 0.2;
//   input.mu = 1.0;
//   input.alpha = 0.0;
//   input.k0 = 0.0;
//   input.sigma0 = 0.0;
// 
//   input.Cp = 1.0;
//   input.Cd = 1.0;
//   input.Rp = 10.0;
//   input.Rd = 10.0;
// 
//   input.Kat = 0.5;
//   input.Kp = 0.1;
//   input.Kar = 0.5;
// 
//   input.pAt = [](Scalar) { return 0.0; };
//   input.pSv = [](Scalar) { return 0.0; };
//   input.u   = [](Scalar) { return 0.0; };
//   input.n0  = [](Scalar) { return 0.0; };
// 
//   {
//     using PassiveEnergy = std::decay_t<decltype(input.passiveEnergy)>;
//     typename PassiveEnergy::Parameters hp;
//     hp.mu1 = 0.0;
//     hp.mu2 = 0.0;
//     hp.C0 = 0.01;
//     hp.C1 = 0.2;
//     hp.C2 = 0.01;
//     hp.C3 = 0.2;
//     input.passiveEnergy = PassiveEnergy(hp);
//   }
// 
//   Model model(input);
//   model.setMaxIterations(100)
//        .setAbsoluteTolerance(1e-10)
//        .setRelativeTolerance(1e-10)
//        .setStepTolerance(1e-12)
//        .setDampingFactor(1.0);
// 
//   Model::State initial;
//   initial.t = 0.0;
//   initial.y = 0.01;
//   initial.v = 0.0;
//   initial.pv = 0.0;
//   initial.par = 0.0;
//   initial.pd = 0.0;
//   initial.ec = 0.0;
//   initial.kc = 0.0;
//   initial.tauc = 0.0;
// 
//   model.initialize(initial);
// 
//   const auto repStatic = model.solveStaticEquilibrium();
//   ASSERT_TRUE(repStatic.converged);
// 
//   const auto repBE = model.step(1e-2);
//   EXPECT_TRUE(repBE.converged);
// 
//   const auto repBDF2 = model.step(1e-2);
//   EXPECT_TRUE(repBDF2.converged);
// }
