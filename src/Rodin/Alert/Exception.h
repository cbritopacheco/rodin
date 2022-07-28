/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ALERT_EXCEPTION_H
#define RODIN_ALERT_EXCEPTION_H

#include <exception>

#include "Alert.h"

namespace Rodin::Alert
{
   /**
    * @brief Derived Alert class representing an exception.
    *
    * Represents an alert which, when raised, will terminate the program after
    * outputting a message with the reason why.
    */
   class Exception : public std::exception, public Alert
   {
      public:
         /**
          * @brief Constructs an Exception with an empty message.
          */
         Exception() = default;

         /**
          * @brief Constructs an Exception with the given message.
          * @param[in] what Description or reason for the Exception being raised.
          */
         Exception(const std::string& what);

         /**
          * @brief Copies the Exception message.
          */
         Exception(const Exception& other)
            :  std::exception(other),
               Alert(other)
         {}

         Exception(Exception&& other)
            :  std::exception(std::move(other)),
               Alert(std::move(other))
         {}

         const char* what() const noexcept override
         {
            return Alert::what();
         }

         /**
          * @brief Raises the exception to the user.
          *
          * Default behaviour is to output a formatted error message and call
          * std::abort.
          */
         virtual void raise() const override;
   };
}

#endif
