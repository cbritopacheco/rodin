// CoupledLV0DCoronary3D.h
#ifndef EXAMPLES_HEART_CORONARYARTERY_COUPLEDLV0DCORONARY3D_H
#define EXAMPLES_HEART_CORONARYARTERY_COUPLEDLV0DCORONARY3D_H

#include <array>
#include <fstream>
#include <map>
#include <string>

#include "Rodin/Heart/CCMLC2014.h"
#include <Rodin/Types.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/MPI.h>
#include <Rodin/Variational.h>
#include <Rodin/PETSc.h>
#include <Rodin/Solver.h>

namespace Rodin::Examples::Heart
{
  class CoupledLV0DCoronary3D
  {
    public:
      using Real = Rodin::Real;
      using Model = Rodin::Heart::CCMLC2014T<>;
      using Attribute = Rodin::Geometry::Attribute;
      using MeshType = Rodin::Geometry::Mesh<Rodin::Context::MPI>;

      using VelocityFESType =
        Rodin::Variational::H1<2, Rodin::Math::SpatialVector<Real>, MeshType>;

      using PressureFESType =
        Rodin::Variational::H1<1, Real, MeshType>;

      using VelocityGridFunctionType =
        Rodin::PETSc::Variational::GridFunction<VelocityFESType>;

      using PressureGridFunctionType =
        Rodin::PETSc::Variational::GridFunction<PressureFESType>;

      using VelocityTrialFunctionType =
        Rodin::PETSc::Variational::TrialFunction<VelocityGridFunctionType, VelocityFESType>;

      using VelocityTestFunctionType =
        Rodin::PETSc::Variational::TestFunction<VelocityFESType>;

      using PressureTrialFunctionType =
        Rodin::PETSc::Variational::TrialFunction<PressureGridFunctionType, PressureFESType>;

      using PressureTestFunctionType =
        Rodin::PETSc::Variational::TestFunction<PressureFESType>;

      using FluxLinearFormType =
        Rodin::Variational::LinearForm<PressureFESType, ::Vec>;

      using LinearSystemType =
        Rodin::PETSc::Math::LinearSystem;

      using FlowProblemType =
        Rodin::Variational::Problem<
          LinearSystemType,
          VelocityTrialFunctionType,
          PressureTrialFunctionType,
          VelocityTestFunctionType,
          PressureTestFunctionType>;

      struct RCR
      {
        Real Rp = 0.0;
        Real C = 0.0;
        Real Rd = 0.0;
        Real pd = 0.0;
        Real pc = 0.0;
        Real pout = 0.0;
      };

      struct Config
      {
        std::string meshPath = "CoronaryArtery.mesh";
        std::string xdmfBasename = "CoronaryArtery";
        std::string csvPath = "CoronaryArtery.csv";

        Attribute wall = 2;
        Attribute inlet = 3;
        std::array<Attribute, 6> outlets{{4, 5, 6, 7, 8, 9}};

        Real meshScale = 1.0e-3;
        Real eps = 1.0e-12;
        Real rho = 1060.0;
        Real mu = 3.5e-3;

        Real dt = 1.0e-3;
        size_t nsteps = 3 * static_cast<int>(0.85 / 1.0e-3);

        RCR defaultRCR{5.0e8, 5.0e-9, 1.0e9, 400.0, 10500.0, 11000.0};
      };

      explicit CoupledLV0DCoronary3D(const Rodin::Context::MPI& context);
      CoupledLV0DCoronary3D(const Rodin::Context::MPI& context, const Config& cfg);
      ~CoupledLV0DCoronary3D();

      CoupledLV0DCoronary3D(const CoupledLV0DCoronary3D&) = delete;
      CoupledLV0DCoronary3D& operator=(const CoupledLV0DCoronary3D&) = delete;

      int run();
      CoupledLV0DCoronary3D& initialize();

      Config& getConfig() noexcept { return m_cfg; }
      const Config& getConfig() const noexcept { return m_cfg; }

      Model& getModel() noexcept { return m_model; }
      const Model& getModel() const noexcept { return m_model; }

    private:
      struct StepData
      {
        Real t = 0.0;
        Real pat = 0.0;
        Real psv = 0.0;

        Real y = 0.0;
        Real v = 0.0;
        Real radius = 0.0;
        Real lvVolume = 0.0;
        Real lvFlow = 0.0;
        Real pv = 0.0;
        Real par = 0.0;
        Real pd = 0.0;

        Real ec = 0.0;
        Real gamma = 0.0;
        Real beta = 0.0;
        Real kc = 0.0;
        Real tauc = 0.0;

        Real qIn = 0.0;
        Real qOutSum = 0.0;
        Real flowBalance = 0.0;

        std::map<Attribute, Real> qOut;
        std::map<Attribute, Real> pc;
        std::map<Attribute, Real> pOut;
      };

      static Model::Input makeInput();
      static MeshType makeMesh(const Rodin::Context::MPI& context, const Config& cfg);

      static void updateRCR(RCR& bc, Real Q, Real dt);
      static void updateRCRNonNew(const Model& model, RCR& bc, Real Q, Real dt);
      static Real periodic_activation(Real t);
      static Real atrial_pressure(Real t);

      void setupModel();
      void setupMeshAndSpaces();
      void setupDiagnostics();
      void printInitialState() const;
      bool isRoot() const;

      bool advance0D();
      void solve3D();
      void computeFluxes();
      void updateHistory();
      void writeOutputs();
      void writeCSVHeader();
      void writeCSVRow();
      StepData collectStepData() const;

      Config m_cfg;
      Model::Input m_input;
      Model m_model;

      MeshType m_mesh;
      Rodin::IO::XDMF m_xdmf;

      VelocityFESType m_uh;
      PressureFESType m_ph;

      VelocityTrialFunctionType m_u;
      PressureTrialFunctionType m_p;

      VelocityTestFunctionType m_v;
      PressureTestFunctionType m_q;

      VelocityGridFunctionType m_uOld;
      PressureGridFunctionType m_pOld;
      PressureGridFunctionType m_one;
      PressureTestFunctionType m_qFlux;

      FluxLinearFormType m_flux;

      FlowProblemType m_flow;
      Rodin::Solver::KSP m_flowKSP;
      Rodin::Solver::SNES m_flowSolver;

      std::map<Attribute, RCR> m_wk;

      mutable StepData m_stepData;
      std::ofstream m_csv;
      bool m_flowFieldSplitsSet = false;
      bool m_initialized = false;
  };
}

#endif
