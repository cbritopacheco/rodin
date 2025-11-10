/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_IO_LOADER_H
#define RODIN_IO_LOADER_H

#include <iostream>
#include <istream>
#include <fstream>

#include <boost/filesystem/path.hpp>

#include "ForwardDecls.h"

namespace Rodin::IO
{
  /**
   * @brief Abstract base class template for loading objects from streams or files.
   *
   * The Loader class provides a common interface for reading objects of type
   * @p T from input streams or files. This class serves as the foundation for
   * specialized loaders that handle different file formats and object types.
   *
   * @tparam T Type of object to load
   *
   * Derived classes must implement:
   * - load(std::istream&) - Load object from an input stream
   * - getObject() - Access the object being loaded
   *
   * ## Usage Example
   * ```cpp
   * class MyObjectLoader : public Loader<MyObject>
   * {
   * public:
   *   void load(std::istream& is) override
   *   {
   *     // Read from stream and populate object
   *     is >> getObject();
   *   }
   *
   * protected:
   *   ObjectType& getObject() override
   *   {
   *     return m_object;
   *   }
   *
   * private:
   *   MyObject m_object;
   * };
   * ```
   *
   * @see MeshLoader, GridFunctionLoader
   */
  template <class T>
  class Loader
  {
    public:
      /**
       * @brief Type of object being loaded.
       */
      using ObjectType = T;

      /**
       * @brief Loads object from an input stream.
       * @param[in] is Input stream to read from
       *
       * This pure virtual method must be implemented by derived classes to
       * read the object data from the provided stream.
       */
      virtual void load(std::istream& is) = 0;

      /**
       * @brief Loads object from a file.
       * @param[in] is Path to the file to load
       *
       * Opens the file at the specified path and delegates to load(std::istream&).
       * This provides a convenient file-based interface while reusing the stream
       * loading logic.
       */
      virtual void load(const boost::filesystem::path& is)
      {
        std::ifstream in(is.c_str());
        load(in);
      }

    protected:
      /**
       * @brief Gets a reference to the object being loaded.
       * @returns Reference to the object
       *
       * This method provides access to the object being populated during loading.
       * Derived classes must implement this to return their internal object reference.
       */
      virtual ObjectType& getObject() = 0;
  };
}

#endif
