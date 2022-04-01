#ifndef RODIN_VARIATIONAL_TESTFUNCTION_H
#define RODIN_VARIATIONAL_TESTFUNCTION_H

#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>

#include "ShapeFunction.h"

namespace Rodin::Variational
{
   template <class FES>
   class TestFunction : public ShapeFunction<FES, Test>
   {
      public:
         TestFunction(FES& fes)
            : ShapeFunction<FES, Test>(fes),
              m_uuid(boost::uuids::random_generator()())
         {}

         TestFunction(const TestFunction& other)
            : ShapeFunction<FES, Test>(other),
              m_uuid(other.m_uuid)
         {}

         TestFunction(TestFunction&& other)
            : ShapeFunction<FES, Test>(std::move(other)),
              m_uuid(std::move(other.m_uuid))
         {}

         boost::uuids::uuid getUUID() const
         {
            return m_uuid;
         }

         TestFunction* copy() const noexcept override
         {
            return new TestFunction(*this);
         }
      private:
         const boost::uuids::uuid m_uuid;
   };
   template <class FES>
   TestFunction(FES& fes) -> TestFunction<FES>;
}
#endif
