#ifndef EXAMPLES_HEART_CORONARYARTERY_COUPLEDLV0DCORONARY3D_H
#define EXAMPLES_HEART_CORONARYARTERY_COUPLEDLV0DCORONARY3D_H

#include <array>
#include <fstream>
#include <map>
#include <memory>
#include <string>

#include "Rodin/Heart/CCMLC2014.h"
#include <Rodin/Types.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/Variational.h>
#include <Rodin/PETSc.h>

namespace Rodin::Examples::Heart
{
  class CoupledLV0DCoronary3D
  {
    public:
      using Real = Rodin::Real;
      using Model = Rodin::Heart::CCMLC2014T<>;
      using Attribute = Rodin::Geometry::Attribute;
      using MeshType = Rodin::Geometry::Mesh<Rodin::Context::Local>;
      using VelocityFESType = Rodin::Variational::H1<2, Rodin::Math::Vector<Real>, MeshType>;
      using PressureFESType = Rodin::Variational::H1<1, Real, MeshType>;
      using VelocityGridFunctionType = Rodin::PETSc::Variational::GridFunction<VelocityFESType>;
      using PressureGridFunctionType = Rodin::PETSc::Variational::GridFunction<PressureFESType>;
      using VelocityTrialFunctionType =
        Rodin::PETSc::Variational::TrialFunction<VelocityGridFunctionType, VelocityFESType>;
      using PressureTrialFunctionType =
        Rodin::PETSc::Variational::TrialFunction<PressureGridFunctionType, PressureFESType>;
      using VelocityTestFunctionType = Rodin::PETSc::Variational::TestFunction<VelocityFESType>;
      using PressureTestFunctionType = Rodin::PETSc::Variational::TestFunction<PressureFESType>;
      using FluxLinearFormType = Rodin::Variational::LinearForm<PressureFESType, ::Vec>;

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

        RCR defaultRCR{1.0e9, 1.0e-10, 5.0e9, 8000.0, 8000.0, 8000.0};
      };

      CoupledLV0DCoronary3D();

      explicit CoupledLV0DCoronary3D(const Config& cfg);

      ~CoupledLV0DCoronary3D();

      CoupledLV0DCoronary3D(const CoupledLV0DCoronary3D&) = delete;

      CoupledLV0DCoronary3D& operator=(const CoupledLV0DCoronary3D&) = delete;

      int run();

      CoupledLV0DCoronary3D& initialize();

      Config& getConfig() noexcept
      {
        return m_cfg;
      }

      const Config& getConfig() const noexcept
      {
        return m_cfg;
      }

      Model& getModel() noexcept
      {
        return *m_model;
      }

      const Model& getModel() const noexcept
      {
        return *m_model;
      }

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

      static void updateRCR(RCR& bc, Real Q, Real dt);
      static void updateRCRNonNew(const Model& model, RCR& bc, Real Q, Real dt);
      static Real periodic_activation(Real t);
      static Real atrial_pressure(Real t);

      void setupModel();
      void setupMeshAndSpaces();
      void setupDiagnostics();
      void printInitialState() const;

      bool advance0D();
      void solve3D();
      void computeFluxesAndUpdateRCR();
      void updateHistory();
      void writeOutputs() const;
      void writeCSVHeader();
      void writeCSVRow();
      StepData collectStepData() const;

      Config m_cfg;
      Model::Input m_input;
      Optional<Model> m_model;

      MeshType m_mesh;
      std::unique_ptr<Rodin::IO::XDMF> m_xdmf;

      std::unique_ptr<VelocityFESType> m_uh;
      std::unique_ptr<PressureFESType> m_ph;

      std::unique_ptr<VelocityTrialFunctionType> m_u;
      std::unique_ptr<PressureTrialFunctionType> m_p;
      std::unique_ptr<VelocityTestFunctionType> m_v;
      std::unique_ptr<PressureTestFunctionType> m_q;

      std::unique_ptr<VelocityGridFunctionType> m_uOld;
      std::unique_ptr<PressureGridFunctionType> m_pOld;
      std::unique_ptr<PressureGridFunctionType> m_one;
      std::unique_ptr<PressureTestFunctionType> m_qFlux;
      std::unique_ptr<FluxLinearFormType> m_flux;

      std::unique_ptr<PressureGridFunctionType> m_muNonNew;

      std::map<Attribute, RCR> m_wk;

      mutable StepData m_stepData;

      std::ofstream m_csv;

      bool m_initialized = false;
  };
}

#endif
