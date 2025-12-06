#ifndef RODIN_ALERT_NOTATION_H
#define RODIN_ALERT_NOTATION_H

#include "Text.h"

namespace Rodin::Alert
{
  /**
   * @brief Type alias for notation foreground color.
   *
   * Notation text uses magenta color for visual distinction in messages.
   */
  using NotationForeground = MagentaT;

  /**
   * @brief Specialized text class for mathematical and technical notation.
   *
   * Provides a colored text formatting class specifically for displaying
   * mathematical notation, variable names, and technical identifiers in
   * alert messages. Notation objects are displayed in magenta to distinguish
   * them from regular message text.
   *
   * This class includes static factory methods for common notation patterns
   * used throughout Rodin, such as polytope identifiers, incidence relations,
   * and numeric values.
   *
   * Example usage:
   * @code{.cpp}
   * Exception() << "Invalid polytope " << Notation::Polytope(2, 5) << Raise;
   * @endcode
   */
  class Notation : public Text<NotationForeground>
  {
    public:
      /// @brief Parent class type alias.
      using Parent = Text<NotationForeground>;

      /**
       * @brief Static arrow notation instance.
       *
       * Constant notation object representing an arrow symbol, commonly
       * used in incidence relations and mathematical expressions.
       */
      static const Notation Arrow;

      /**
       * @brief Constructs notation from a string.
       * @param text The notation text.
       */
      Notation(const std::string& text)
        : Parent(text)
      {}

      /**
       * @brief Constructs notation from a C-style string.
       * @param out The notation text as a null-terminated string.
       */
      Notation(const char* out)
        : Parent(std::string(out))
      {}

      /**
       * @brief Constructs notation from any type convertible to string.
       * @tparam T The type to convert (must support std::to_string).
       * @param out The value to convert to notation text.
       */
      template <class T>
      Notation(const T& out)
        : Parent(std::to_string(out))
      {}

      /**
       * @brief Copy constructor.
       */
      Notation(const Notation& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor.
       */
      Notation(Notation&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Factory method for text notation.
       * @param text The text content.
       * @return Notation object containing the text.
       */
      inline
      static Notation Text(const std::string& text)
      {
        return text;
      }

      /**
       * @brief Factory method for numeric notation.
       * @tparam T Arithmetic type.
       * @param out The numeric value.
       * @return Notation object containing the formatted number.
       */
      template <class T>
      inline
      static Notation Number(const T& out)
      {
        static_assert(std::is_arithmetic_v<T>);
        std::stringstream ss;
        ss << out;
        return ss.str();
      }

      /**
       * @brief Factory method for general printable values.
       * @tparam T Type supporting stream insertion.
       * @param out The value to print.
       * @return Notation object containing the formatted value.
       */
      template <class T>
      inline
      static Notation Print(const T& out)
      {
        std::stringstream ss;
        ss << out;
        return ss.str();
      }

      /**
       * @brief Factory method for incidence relation notation.
       * @param d First dimension.
       * @param dp Second dimension.
       * @return Notation object representing the incidence relation.
       *
       * Creates notation for mesh incidence relations in the form
       * "Incidence(d → dp)" where the arrow represents topological
       * adjacency.
       */
      inline
      static Notation Incidence(size_t d, size_t dp)
      {
        std::stringstream ss;
        ss << "Incidence(" << d << " " << Arrow << " " << dp << ")";
        return ss.str();
      }

      /**
       * @brief Factory method for polytope notation.
       * @param d Dimension of the polytope.
       * @param index Index of the polytope.
       * @return Notation object representing the polytope.
       *
       * Creates notation for polytope identifiers in the form
       * "Polytope(d, index)".
       */
      inline
      static Notation Polytope(size_t d, Index index)
      {
        std::stringstream ss;
        ss << "Polytope(" << d << ", " << index << ")";
        return ss.str();
      }

      /**
       * @brief Factory method for boolean true notation.
       * @return Notation object containing "true".
       */
      inline
      static Notation True()
      {
        return Notation("true");
      }

      /**
       * @brief Factory method for boolean false notation.
       * @return Notation object containing "false".
       */
      inline
      static Notation False()
      {
        return Notation("false");
      }

      /**
       * @brief Factory method for predicate notation.
       * @param v Boolean value of the predicate.
       * @param pred Name of the predicate.
       * @return Notation object representing the predicate evaluation.
       *
       * Creates notation showing a predicate and its boolean value in
       * the form "[pred] = true" or "[pred] = false".
       */
      inline
      static Notation Predicate(bool v, const std::string& pred)
      {
        std::stringstream ss;
        ss << "[" << pred <<  "] = " << (v ? "true" : "false");
        return ss.str();
      }

      /**
       * @brief Factory method for set notation.
       * @tparam Iterator Iterator type for the set elements.
       * @param first Iterator to the first element.
       * @param last Iterator past the last element.
       * @return Notation object representing the set.
       *
       * Creates notation for sets in the form "{ elem1, elem2, ... }".
       */
      template <class Iterator>
      static Notation Set(Iterator first, Iterator last)
      {
        std::stringstream ss;
        ss << "{ ";
        for (Iterator it = first; it != last; ++it)
        {
          ss << *it;
          if (std::next(it) != last)
              ss << ", ";
        }
        ss << " }";
        return ss.str();
    }
  };
}

#endif
