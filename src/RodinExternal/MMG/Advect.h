/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EXTERNAL_MMG_ADVECT_H
#define RODIN_EXTERNAL_MMG_ADVECT_H

#include <vector>
#include <optional>

#include "Rodin/Variational/GridFunction.h"
#include "Rodin/Variational/FiniteElementSpace.h"

#include "Common.h"
#include "ForwardDecls.h"
#include "ISCDProcess.h"

namespace Rodin::External::MMG
{
  template <class FEC>
  class Advect
  {
    public:
      /**
       * @brief Constructs an Advect object
       * @param[in, out] ls Function to advect
       * @param[in] disp Displacement velocity field
       */
      Advect(
          Variational::GridFunction<FEC, Traits::Serial>& ls,
          const Variational::GridFunction<FEC, Traits::Serial>& disp)
        : m_t(0),
          m_avoidTrunc(false),
          m_ex(true),
          m_ls(ls),
          m_disp(disp),
          m_advect(getISCDAdvectExecutable())
      {
        assert(ls.getFiniteElementSpace().getVectorDimension() == 1);
        assert(
            disp.getFiniteElementSpace().getVectorDimension()
            == ls.getFiniteElementSpace().getMesh().getSpaceDimension());
      }

      /**
       * @brief Specifies whether to enable or disable truncation of time in the CFL condition.
       *
       * If true, avoids truncation of the time period for advection due to
       * the CFL condition. Otherwise, no action is taken.
       *
       * By default it is set to false.
       */
      Advect& avoidTimeTruncation(bool avoidTrunc = true)
      {
        m_avoidTrunc = avoidTrunc;
        return *this;
      }

      /**
       * @brief Specifies whether to extrapolate characteristic lines or not.
       *
       * If true, characteristic lines are extrapolated outside the domain when
       * the input velocity field causes them to do so. Otherwise, no action is
       * taken.
       *
       * By default it is set to true.
       */
      Advect& enableExtrapolation(bool ex = true)
      {
        m_ex = ex;
        return *this;
      }

      /**
       * @brief Advances the level set function by the time step `dt`.
       * @param[in] dt Advection timestep
       */
      void step(double dt)
      {
        assert(!std::isnan(dt) && !std::isinf(dt));

        auto& mesh = m_ls.getFiniteElementSpace().getMesh();

        auto meshp = m_advect.tmpnam(".mesh", "RodinMMG");
        mesh.save(meshp, IO::MeshFormat::MEDIT);

        auto solp = m_advect.tmpnam(".sol", "RodinMMG");
        m_ls.save(solp, IO::GridFunctionFormat::MEDIT);

        auto dispp = m_advect.tmpnam(".sol", "RodinMMG");
        m_disp.save(dispp, IO::GridFunctionFormat::MEDIT);

        auto outp = m_advect.tmpnam(".sol", "RodinMMG");

        int retcode = 1;
        if (mesh.isSurface())
        {
          if (!m_avoidTrunc)
          {
            Alert::Warning()
              << "Time truncation is always avoided for surface meshes."
              << Alert::Raise;
          }

          retcode = m_advect.run(
              meshp.string(),
              "-dt", std::to_string(dt),
              m_ex ? "" : "-noex",
              "-c", solp.string(),
              "-s", dispp.string(),
              "-o", outp.string(),
              "-nocfl");
        }
        else
        {
          retcode = m_advect.run(
              meshp.string(),
              "-dt", std::to_string(dt),
              m_avoidTrunc ? "-nocfl" : "",
              m_ex ? "" : "-noex",
              "-c", solp.string(),
              "-s", dispp.string(),
              "-o", outp.string());
        }

        if (retcode != 0)
          Alert::Exception("ISCD::Advection invocation failed.").raise();

        m_ls.load(outp, IO::GridFunctionFormat::MEDIT);

        m_t += dt;
      }

    private:
      double m_t;
      bool m_avoidTrunc, m_ex;
      Variational::GridFunction<FEC, Traits::Serial>& m_ls;
      const Variational::GridFunction<FEC, Traits::Serial>& m_disp;
      ISCDProcess m_advect;
  };
}

#endif
