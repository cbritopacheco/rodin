/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_FORMLANGUAGE_BASE_H
#define RODIN_FORMLANGUAGE_BASE_H

#include <deque>
#include <memory>
#include <cassert>
#include <variant>
#include <typeinfo>

#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>

#include "Rodin/Types.h"
#include "Rodin/Math/ForwardDecls.h"
#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Variational/BasisOperator.h"

namespace Rodin::FormLanguage
{
  /**
   * @brief Base class for all classes which are part of Variational::FormLanguage.
   */
  class Base
  {
    static boost::uuids::random_generator s_gen;

    public:
      inline
      Base()
        : m_uuid(s_gen())
      {}

      inline
      Base(const Base& other)
        : m_uuid(other.m_uuid)
      {}

      inline
      Base(Base&& other)
        : m_uuid(std::move(other.m_uuid))
      {}

      Base& operator=(const Base&) = delete;

      Base& operator=(Base&&) = delete;

      inline
      const boost::uuids::uuid& getUUID() const
      {
        return m_uuid;
      }

      /**
       * @brief Virtual destructor.
       */
      virtual ~Base() = default;

      virtual const char* getName() const
      {
        return "Rodin::FormLanguage::Base";
      }

      inline
      constexpr
      const Math::Matrix& object(Math::Matrix&& obj) const
      {
        return m_mobjs.emplace_back(std::move(obj));
      }

      inline
      constexpr
      const Math::Vector& object(Math::Vector&& obj) const
      {
        return m_vobjs.emplace_back(std::move(obj));
      }

      inline
      constexpr
      const auto& object(Variational::TensorBasis<Math::Vector>&& obj) const
      {
        return m_tbvobjs.emplace_back(std::move(obj));
      }

      template <class EigenDerived>
      inline
      constexpr
      const auto& object(Variational::TensorBasis<Math::Matrix>&& obj) const
      {
        return m_tbmobjs.emplace_back(std::move(obj));
      }

      template <class ObjectType>
      inline
      constexpr
      const ObjectType& object(const ObjectType& obj) const
      {
        return obj;
      }

      inline
      constexpr
      Scalar object(Scalar s) const
      {
        return s;
      }

      inline
      constexpr
      const auto& object(const Variational::TensorBasis<Scalar>& obj) const
      {
        return obj;
      }

      template <class EigenDerived>
      inline
      constexpr
      const auto& object(const Variational::TensorBasis<Eigen::MatrixBase<EigenDerived>>& obj) const
      {
        return obj;
      }

      /**
       * @internal
       * @brief Copies the object and returns a non-owning pointer to the
       * copied object.
       * @returns Non-owning pointer to the copied object.
       */
      virtual Base* copy() const noexcept = 0;

    private:
      const boost::uuids::uuid m_uuid;

      mutable std::deque<Math::Vector> m_vobjs;
      mutable std::deque<Math::Matrix> m_mobjs;
      mutable std::deque<Variational::TensorBasis<Math::Vector>> m_tbvobjs;
      mutable std::deque<Variational::TensorBasis<Math::Matrix>> m_tbmobjs;
  };
}

#endif
