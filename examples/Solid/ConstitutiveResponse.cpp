/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
//! [constitutive-source]
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <Rodin/Solid.h>

using namespace Rodin;

int main()
{
  //! [constitutive-setup]
  constexpr Real lambda = 2;
  constexpr Real mu = 1;
  constexpr Real increment = 1e-6;
  constexpr Real tolerance = 1e-8;
  Solid::NeoHookean law(lambda, mu);
  std::ofstream samples("ConstitutiveResponse.csv");
  samples.precision(17);
  samples << "stretch,energy,stress,finite_difference\n";
  //! [constitutive-setup]

  Real maximumError = 0;
  for (size_t i = 0; i <= 70; ++i)
  {
    const Real stretch = 0.7 + 0.01 * i;
    std::array<Real, 3> energy{};
    Real stress = 0;
    //! [constitutive-bind]
    for (size_t sample = 0; sample < 3; ++sample)
    {
      const Real value = stretch + (Real(sample) - 1) * increment;
      Math::SpatialMatrix<Real> gradient(3, 3);
      gradient.setZero();
      gradient(0, 0) = value - 1;
      Solid::KinematicState state(3);
      state.setDisplacementGradient(gradient);
      Solid::ConstitutivePoint point(state);
      law.setConstitutivePoint(point);
      energy[sample] = law.getStrainEnergyDensity();
      if (!std::isfinite(energy[sample]))
        return EXIT_FAILURE;
      if (sample == 1)
      {
        Math::SpatialMatrix<Real> firstPiola;
        law.getFirstPiolaKirchhoffStress(firstPiola);
        stress = firstPiola(0, 0);
        if (!std::isfinite(stress))
          return EXIT_FAILURE;
      }
    }
    //! [constitutive-bind]
    const Real derivative = (energy[2] - energy[0]) / (2 * increment);
    const Real error = std::abs(derivative - stress) / std::max(Real(1), std::abs(stress));
    maximumError = std::max(maximumError, error);
    samples << stretch << ',' << energy[1] << ',' << stress << ',' << derivative << '\n';
  }
  samples.flush();
  std::cout << "Maximum scaled energy derivative error: " << maximumError << std::endl;
  return samples && maximumError < tolerance ? EXIT_SUCCESS : EXIT_FAILURE;
}
//! [constitutive-source]
