#ifndef RODIN_VARIATIONAL_QUADRATURERULE_H
#define RODIN_VARIATIONAL_QUADRATURERULE_H

#include <vector>
#include <utility>

#include "Rodin/Math.h"
#include "Rodin/Geometry.h"

#include "MFEM.h"

namespace Rodin::Variational
{
  class QuadratureRule
  {
    using Key = std::pair<Geometry::Type, size_t>;
    static std::map<Key, QuadratureRule> s_rules;

    struct ValueType
    {
      Scalar weight;
      Math::Vector point;
    };

    public:

      static const QuadratureRule& get(Geometry::Type geometry, size_t order)
      {
        Key key{geometry, order};
        auto search = s_rules.lower_bound(key);
        if (search != s_rules.end() && !(s_rules.key_comp()(key, search->first)))
        {
          // key already exists
          return search->second;
        }
        else
        {
          size_t dim;
          switch (geometry)
          {
            case Geometry::Type::Point:
            {
              dim = 0;
              break;
            }
            case Geometry::Type::Segment:
            {
              dim = 1;
              break;
            }
            case Geometry::Type::Square:
            case Geometry::Type::Triangle:
            {
              dim = 2;
              break;
            }
            case Geometry::Type::Cube:
            case Geometry::Type::Prism:
            case Geometry::Type::Pyramid:
            case Geometry::Type::Tetrahedron:
            {
              dim = 3;
              break;
            }
          }
          const mfem::IntegrationRule& ir = mfem::IntRules.Get(static_cast<int>(geometry), order);
          auto it = s_rules.insert(search, {key, QuadratureRule(ir, dim)});
          return it->second;
        }
      }

      QuadratureRule(const mfem::IntegrationRule& ir, size_t dim)
      {
        m_points.reserve(ir.GetNPoints());
        for (int i = 0; i < ir.GetNPoints(); i++)
        {
          const mfem::IntegrationPoint& ip = ir.IntPoint(i);
          m_points.push_back(ValueType{ip.weight, Internal::ip2vec(ip, dim)});
        }
      }

      QuadratureRule(const QuadratureRule&) = default;

      QuadratureRule(QuadratureRule&&) = default;

      inline
      size_t size() const
      {
        return m_points.size();
      }

      inline
      Scalar getWeight(size_t i) const
      {
        assert(i < m_points.size());
        return m_points[i].weight;
      }

      inline
      const Math::Vector& getPoint(size_t i) const
      {
        assert(i < m_points.size());
        return m_points[i].point;
      }

    private:
      std::vector<ValueType> m_points;
  };
}

#endif

