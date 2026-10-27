/*
 * Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE or https://www.boost.org/LICENSE_1_0.txt)
 */
#include <memory>
#include <gtest/gtest.h>
#include "Rodin/IO/Printer.h"
#include "Rodin/IO/Loader.h"

namespace
{
  class TrackedPrinter final : public Rodin::IO::Printer<int>
  {
    public:
      TrackedPrinter(const int& object, bool& destroyed)
        : m_object(object),
          m_destroyed(destroyed)
      {}

      ~TrackedPrinter() override
      {
        m_destroyed = true;
      }

      void print(std::ostream&) override {}
      const int& getObject() const override
      {
        return m_object;
      }

    private:
      const int& m_object;
      bool& m_destroyed;
  };

  class TrackedLoader final : public Rodin::IO::Loader<int>
  {
    public:
      TrackedLoader(int& object, bool& destroyed)
        : m_object(object),
          m_destroyed(destroyed)
      {}

      ~TrackedLoader() override
      {
        m_destroyed = true;
      }
      void load(std::istream&) override {}

    protected:
      int& getObject() override
      {
        return m_object;
      }

    private:
      int& m_object;
      bool& m_destroyed;
  };

  TEST(Rodin_IO_Printer, DestroysDerivedPrinterThroughBase)
  {
    bool destroyed = false;
    const int value = 1;
    {
      std::unique_ptr<Rodin::IO::Printer<int>> printer =
        std::make_unique<TrackedPrinter>(value, destroyed);
      EXPECT_EQ(printer->getObject(), value);
    }
    EXPECT_TRUE(destroyed);
  }

  TEST(Rodin_IO_Loader, DestroysDerivedLoaderThroughBase)
  {
    bool destroyed = false;
    int value = 1;
    {
      std::unique_ptr<Rodin::IO::Loader<int>> loader =
        std::make_unique<TrackedLoader>(value, destroyed);
    }
    EXPECT_TRUE(destroyed);
  }
}
