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
  /**
   * @brief Advection of a level set function by a velocity field.
   *
   * Solves the advection equation
   * @f[
   * \left\{
   *   \begin{aligned}
   *     \dfrac{\partial u}{\partial t} + v(x) \cdot \nabla u (t, x) &= 0
   *       && \text{ in } \Omega \times (0, + \infty) \\
   *     u(x, 0) &= u_0(x) && \text{ on } \Omega \times \{ t = 0 \}
   *   \end{aligned}
   * \right.
   * @f]
   * where @f$ u_0 : \mathbb{R}^d \rightarrow \mathbb{R} @f$ is known.
   */
  class Advect
  {
    public:
      using LevelSetFunction = Variational::GridFunction<Variational::H1<Scalar, Context::Serial>>;
      using VectorField = Variational::GridFunction<Variational::H1<Math::Vector, Context::Serial>>;

      /**
       * @brief Constructs an Advect object
       * @param[in, out] ls Function to advect
       * @param[in] disp Displacement velocity field
       */
      Advect(LevelSetFunction& ls, const VectorField& disp)
        : m_t(0),
          m_ex(true),
          m_advectTheSurface(false),
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

        auto& mesh = m_ls.get().getFiniteElementSpace().getMesh();

        auto meshp = m_advect.tmpnam(".mesh", "RodinMMG");
        mesh.save(meshp, IO::FileFormat::MEDIT);

        auto solp = m_advect.tmpnam(".sol", "RodinMMG");
        m_ls.get().save(solp, IO::FileFormat::MEDIT);

        auto dispp = m_advect.tmpnam(".sol", "RodinMMG");
        m_disp.get().save(dispp, IO::FileFormat::MEDIT);

        auto outp = m_advect.tmpnam(".sol", "RodinMMG");

        int retcode = 1;
        if (mesh.isSurface())
        {
          retcode = m_advect.run(
              meshp.string(),
              "-surf",
              m_ex ? "" : "-noex",
              "-c", solp.string(),
              "-s", dispp.string(),
              "-out", outp.string(),
              "-dt", std::to_string(dt),
              "-nocfl");
        }
        else
        {
          if (mesh.getSpaceDimension() == 2)
          {
            retcode = m_advect.run(
                meshp.string(),
                "-dt", std::to_string(dt),
                m_ex ? "" : "-noex",
                "-c", solp.string(),
                "-s", dispp.string(),
                "-o", outp.string()
                );
          }
          else if (mesh.getSpaceDimension() == 3)
          {
            retcode = m_advect.run(
                meshp.string(),
                "-dt", std::to_string(dt),
                m_ex ? "" : "-noex",
                "-c", solp.string(),
                "-s", dispp.string(),
                "-o", outp.string(),
                "-nocfl"
                );
          }
          else
          {
            Alert::Exception() << "Invalid dimension" << Alert::Raise;
          }
        }

        if (retcode != 0)
          Alert::Exception(
              "MMG::Advect: ISCD::Advection invocation failed.").raise();

        m_ls.get().load(outp, IO::FileFormat::MEDIT);

        m_t += dt;
      }

      Advect& surface(bool advectTheSurface = true)
      {
        m_advectTheSurface = advectTheSurface;
        return *this;
      }

    private:
      double m_t;
      bool m_ex;
      bool m_advectTheSurface;
      std::reference_wrapper<LevelSetFunction> m_ls;
      std::reference_wrapper<const VectorField> m_disp;
      ISCDProcess m_advect;
  };
}

#endif
