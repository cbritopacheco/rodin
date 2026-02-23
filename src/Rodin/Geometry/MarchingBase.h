/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_MARCHINGBASE_H
#define RODIN_GEOMETRY_MARCHINGBASE_H

#include <functional>
#include "Rodin/Configure.h"
#include "Rodin/Geometry/Types.h"
#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Variational/P1/P1.h"

namespace Rodin::Geometry
{
  template <class ... Params>
  class MarchingBase;

  template <class Mesh, class Data>
  class MarchingBase<Variational::GridFunction<Variational::P1<Real, Mesh>, Data>>
  {
    public:
      using MeshType = Mesh;

      using FESType = Variational::P1<Real, Mesh>;

      using GridFunctionType = Variational::GridFunction<FESType, Data>;

      struct NoSplitT {};
      static constexpr NoSplitT NoSplit{};

      struct Split
      {
        Attribute negative;
        Attribute positive;
      };

      using SplitMap = FlatMap<Attribute, std::variant<Split, NoSplitT>>;

      MarchingBase(const GridFunctionType& gf)
        : m_gf(gf),
          m_eps(1e-12)
      {
        m_split.resize(gf.getFiniteElementSpace().getMesh().getDimension() + 1);
      }

      MarchingBase& setSplit(size_t d, const Attribute& attr, const std::variant<Split, NoSplitT>& value)
      {
        assert(d < m_split.size());
        m_split[d][attr] = value;
        return *this;
      }

      MarchingBase& split(size_t d, Attribute attr, const Split& split)
      {
        return this->setSplit(d, attr, split);
      }

      MarchingBase& noSplit(size_t d, Attribute attr)
      {
        return this->setSplit(d, attr, NoSplit);
      }

      MarchingBase& setNegative(const Optional<Attribute>& attr)
      {
        m_negative = attr;
        return *this;
      }

      MarchingBase& setPositive(const Optional<Attribute>& attr)
      {
        m_positive = attr;
        return *this;
      }

      MarchingBase& setInterface(const Optional<Attribute>& attr)
      {
        m_interface = attr;
        return *this;
      }

      const SplitMap& getSplitMap(size_t d) const
      {
        assert(d < m_split.size());
        return m_split[d];
      }

      const Optional<Attribute>& getNegative() const
      {
        return m_negative;
      }

      const Optional<Attribute>& getPositive() const
      {
        return m_positive;
      }

      const Optional<Attribute>& getInterface() const
      {
        return m_interface;
      }

      MarchingBase& setTolerance(Real eps)
      {
        m_eps = eps;
        return *this;
      }

      Real getTolerance() const
      {
        return m_eps;
      }

      const GridFunctionType& getGridFunction() const
      {
        return m_gf.get();
      }

      virtual MeshType discretize() const = 0;

    private:
      std::reference_wrapper<const GridFunctionType> m_gf;

      Optional<Attribute> m_negative;
      Optional<Attribute> m_positive;
      Optional<Attribute> m_interface;

      Real m_eps;

      std::vector<SplitMap> m_split;
  };
}

#endif
