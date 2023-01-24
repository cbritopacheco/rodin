#ifndef RODIN_VARIATIONAL_ESSENTIALBOUNDARY_H
#define RODIN_VARIATIONAL_ESSENTIALBOUNDARY_H

#include <map>
#include <memory>
#include <variant>

#include "Rodin/Utility/Overloaded.h"

#include "H1.h"
#include "Component.h"
#include "DirichletBC.h"
#include "TrialFunction.h"
#include "ScalarFunction.h"
#include "VectorFunction.h"

namespace Rodin::Variational
{
  /**
   * @brief Represents the essential boundary of a variational problem.
   */
  class EssentialBoundary
  {
    struct Value
    {
      std::unique_ptr<FunctionBase> value;
      std::set<int> attributes;
    };

    public:
      EssentialBoundary() = default;

      EssentialBoundary(const EssentialBoundary& other)
      {
        for (auto it = other.m_tfVal.begin(); it != other.m_tfVal.end(); it++)
        {
          auto uuid = it->first;
          const auto& tfValue = it->second;
          m_tfVal[uuid] =
            Value{
              std::unique_ptr<FunctionBase>(tfValue.value->copy()),
              tfValue.attributes};
        }

        for (const auto& [uuid, valueMap] : other.m_tfCompVal)
        {
          for (const auto& [idx, compValue] : valueMap)
          {
            m_tfCompVal[uuid][idx] =
              Value{
                std::unique_ptr<FunctionBase>(compValue.value->copy()),
                compValue.attributes};
          }
        }
      }

      EssentialBoundary(EssentialBoundary&& other)
        : m_tfVal(std::move(other.m_tfVal)),
          m_tfCompVal(std::move(other.m_tfCompVal))
      {}

      EssentialBoundary& operator=(EssentialBoundary&&) = default;

      template <class Trait>
      void add(const DirichletBC<TrialFunction<H1<Trait>>>& dbc)
      {
        assert(dbc.getTrialFunction().getRangeType() == dbc.getValue().getRangeType());
        m_tfVal[dbc.getTrialFunction().getUUID()] =
          Value{
            std::unique_ptr<FunctionBase>(dbc.getValue().copy()),
            dbc.getBoundaryAttributes()};
      }

      template <class Trait>
      void add(const DirichletBC<Component<TrialFunction<H1<Trait>>>>& dbc)
      {
        assert(dbc.getValue().getRangeType() == RangeType::Scalar);
        auto uuid = dbc.getComponent().getTrialFunction().getUUID();
        m_tfCompVal[uuid][dbc.getComponent().getIndex()] =
          Value{
            std::unique_ptr<FunctionBase>(dbc.getValue().copy()),
            dbc.getBoundaryAttributes()};
      }

      const std::map<boost::uuids::uuid, Value>& getTFMap() const
      {
        return m_tfVal;
      }

      const std::map<boost::uuids::uuid, std::map<int, Value>>& getTFCompMap() const
      {
        return m_tfCompVal;
      }

    private:
      std::map<boost::uuids::uuid, Value> m_tfVal;
      std::map<boost::uuids::uuid, std::map<int, Value>> m_tfCompVal;
  };
}

#endif
