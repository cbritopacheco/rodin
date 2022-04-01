#ifndef RODIN_VARIATIONAL_FORMLANGUAGE_BASE_H
#define RODIN_VARIATIONAL_FORMLANGUAGE_BASE_H

#include <memory>
#include <cassert>

#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>

namespace Rodin::Variational::FormLanguage
{
   /**
    * @brief Base class for all classes which are part of Variational::FormLanguage.
    */
   class Base
   {
      public:
         Base()
            : m_uuid(boost::uuids::random_generator()())
         {}

         Base(const Base& other)
            : m_uuid(other.m_uuid)
         {}

         Base(Base&& other)
            : m_uuid(std::move(other.m_uuid))
         {}

         /**
          * @brief Virtual destructor.
          */
         virtual ~Base() = default;

         /**
          * @internal
          * @brief Copies the object and returns a non-owning pointer to the
          * copied object.
          * @returns Non-owning pointer to the copied object.
          */
         virtual Base* copy() const noexcept = 0;

         boost::uuids::uuid getUUID() const
         {
            return m_uuid;
         }
      private:
         const boost::uuids::uuid m_uuid;
   };

   template <class InternalValue>
   class Buildable : public Base
   {
      public:
         virtual std::unique_ptr<InternalValue> build() const = 0;
   };
}

#endif
