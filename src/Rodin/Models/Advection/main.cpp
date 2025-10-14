//   ______ ______ _      _  _____      ______
//  |  ____|  ____| |    (_)/ ____|    |  ____|
//  | |__  | |__  | |     _| (___   ___| |__
//  |  __| |  __| | |    | |\___ \ / __|  __|
//  | |    | |____| |____| |____) | (__| |____
//  |_|    |______|______|_|_____/ \___|______|
//  Finite Elements for Life Sciences and Engineering
//
//  License:     LGL2.1 License
//  FELiScE default license: LICENSE in root folder
//
//  Main authors:    C. Brito-Pacheco
//
//  Usage:
//
//  To customize the parameters, we have to use a combination of the data file
//  that is provided along with the code and modify the input parameters of the
//  solver. Notably the most important parameters are:
//  - The pressures of the pulmonary veins and left ventricle, which can be
//  set using data curves.
//  - The conductance of the mitral valve and pulmonary veins, which can be
//  set using the solver input parameters.
//  - The solver input parameters also allow for setting the endocardium,
//  orifices, and annuli, which are important for the simulation.
//  - Whether to use shell elements or not, which can be set using the `shell`
//  parameter.
#include "Model/Heart/LV0D_LA3D_Model.hpp"

#include "Data/Heart/Pressure.hpp"

#include "curveGen.hpp"

using namespace felisce;

// Pulmonary veins pressure
class PVPressure : public Heart::PressureBase
{
  public:
    PVPressure()
      : m_dc("LeftAtrium/DataCurve/RightVentriclePressure.dat"),
        m_cg(m_dc, 2000)
    {
      m_cg.extendByPeriodicity(5); // Extend the curve by periodicity
    }

    double operator()(felReal time) const override
    {
      if (time < 0.0015)
        return 1300;
      else
        return 0.175 * m_cg(time) + 1100;
    }

  private:
    mutable DataCurve m_dc;
    mutable CurveGen m_cg;
};

// Left ventricle pressure
class LVPressure : public Heart::PressureBase
{
  public:
    LVPressure()
      : m_dc("LeftAtrium/DataCurve/LeftVentriclePressure.dat"),
        m_cg(m_dc, 2000)
    {
      m_cg.extendByPeriodicity(5); // Extend the curve by periodicity
    }

    double operator()(felReal time) const override
    {
      return m_cg(time + 0.015) - 200;
    }

  private:
    mutable DataCurve m_dc;
    mutable CurveGen m_cg;
};

int main(const int argc, const char** argv)
{
  std::ofstream out("LV0D_LA3D_Model_Shell.log");

  CommandLineOption opt(argc, argv);
  FelisceParam::instance().initialize(opt);
  FelisceTransient fstransient;

  LVPressure left_ventricle_pressure;

  Heart::PulmonaryVeinsPressure pulmonary_veins_pressure;
  Heart::VenaCavaPressure vena_cava_pressure;

  // PVPressure pulmonary_veins_pressure;

  if (MpiInfo::rankProc() == 0)
  {
    std::cout << "Initial Left Ventricle Pressure: " << pulmonary_veins_pressure(0.0) << " Pa" << std::endl;
    std::cout << "Initial Pulmonary Veins Pressure: " << pulmonary_veins_pressure(0.0) << " Pa" << std::endl;
  }

  Heart::LV0D_LA3D_Model::Input input =
  {
    .opt = opt,
    .fstransient = fstransient,
    .base = {
      // Preloading parameters for atrium
      .load_step = {
        // Starting load
        .start = 0.0,
        // Ending load
        .end = pulmonary_veins_pressure(0.0),
        // Step size between loads
        .step = 5
      },
      .log_stream = { out },
      .active_law =
        Heart::ActiveSolver<Heart::ChapelleMoireauActiveLaw>::Input
        {
          .Es = 3e5,
          .initFibDef = 0,
          .degreeOfExactness = FelisceParam::instance().degreeOfExactness[0],
          .rheology_policy = {
            .Es = 3e5,
            .initFibDef = 0,
            .initActiveStiffness = 0,
            .initActiveStress = 0,
            .DampingParallel = 50,
            .DestructionRate = 3.0,
            .CrossBridgeStiffness = 1e5,
            .Contractility = 8.5e5,
            .alphaElecActivation = 1,
            .betaElecActivation = 0,
            .electricalActivation = {
              .typeElectricalActivation = "Prescribed",
              .initialElectricalActivation = 0,
              .heartBeatDuration = FELISCE_DATA_HEART_DEFAULT_HEARTBEAT_DURATION,
              .delay = 0,
              .depolarizationDuration = 0,
              .repolarizationDuration = 0.4,
              .plateauDuration = 0.1,
              .waveVelocity = 0.7,
              .umax = 1,
              .umin = 0,
              .elevation = 1.0,
              .alphaElecActivation = 1,
              .betaElecActivation = 0
            }
          }
        }
    },
    // Solver input parameters
    .solver = {
      // Endocardium label
      .endocardium = { 1 },
      // Orifices labels
      .orifices = { 50, 51, 52, 53, 54 },
      // Annuli labels
      .annuli = { 50, 51, 52, 53, 54 },
      // Pulmonary vein parameters
      .pulmonary_veins = {
        .conductance = 5e-7,
        .pressure = pulmonary_veins_pressure
      },
      // Mitral valve parameters
      .mitral_valve = {
        .conductance = 1.5e-6
      }
    },
    .shell = true,
    // Preloading parameters for the lumped sphere
    .load_step = {
      .start = 0.0,
      .end = pulmonary_veins_pressure(0.0),
      .step = 5
    },
    .lumped_sphere = {
      .R0 = 2.36e-2,
      .d0 = 1.42e-2,
      .mu_1 = 0,
      .mu_2 = 0,
      .C_0 = 1.9e3,
      .C_1 = 1.1e-1,
      .C_2 = 1.9e3,
      .C_3 = 1.1e-1,
      .bulk = 1e5,
      .density = 1000,

      .R_p = 8e6,
      .R_d = 1e8,
      .C_p = 2.5e-8,
      .C_d = 1e-8,
      .p_veins = vena_cava_pressure,

      .K_cavities = 9e-6,
      .K_closed = 5e-12,
      .K_artery = 1.3e-5,

      .init_atrium = pulmonary_veins_pressure(0.0),
      .init_cavity = pulmonary_veins_pressure(0.0),
      .init_artery = 11000,
      .init_distal = 10000,

      .R_artery = 1.3e-5,

      .damping = 70,

      .Es = 3e7,
      .initFibDef = 0,
      .initActiveStiffness = 0,
      .initActiveStress= 0,
      .DampingParallel = 10,
      .DestructionRate = 2.0,
      .CrossBridgeStiffness = 1e5,
      .Contractility = 1.24e5,
      .alphaElecActivation = 1,
      .betaElecActivation = 0,
      .typeElectricalActivation = "Analytical",
      .initialElectricalActivation = 0,
      .heartBeatDuration = FELISCE_DATA_HEART_DEFAULT_HEARTBEAT_DURATION,
      .delay = 0.12,
      .depolarizationDuration = 0.065,
      .repolarizationDuration = 0.07,
      .plateauDuration = 0.01,
      .waveVelocity = 2.0,
      .umax = 35,
      .umin = -12,
      .elevation = 1.0
    }
  };

  Heart::LV0D_LA3D_Model model(input);
  model.SolveDynamicProblem();

  return 0;
}


