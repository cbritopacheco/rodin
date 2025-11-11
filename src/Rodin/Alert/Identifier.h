#ifndef RODIN_ALERT_IDENTIFIER_H
#define RODIN_ALERT_IDENTIFIER_H

#include "Text.h"

namespace Rodin::Alert
{
  /**
   * @brief Templated text class for code identifiers with bold formatting.
   * @ingroup AlertModule
   * @tparam Foreground The foreground color type.
   *
   * Specialized text class for displaying code identifiers (class names,
   * namespace names, function names) in alert messages. Identifiers are
   * automatically formatted with bold text for emphasis.
   */
  template <class Foreground>
  class IdentifierT : public Text<Foreground>
  {
    public:
      /// @brief Parent class type alias.
      using Parent = Text<Foreground>;

      /**
       * @brief Constructs an identifier with the given text.
       * @param text The identifier text.
       *
       * Automatically applies bold formatting to the text.
       */
      IdentifierT(const std::string& text)
        : Parent(text)
      {
        this->setBold();
      }
  };

  /**
   * @brief Namespace containing factory functions for code identifiers.
   * @ingroup AlertModule
   *
   * Provides convenient factory functions for creating colored, bold
   * identifiers for different code elements (classes, namespaces, functions).
   */
  namespace Identifier
  {
    /**
     * @brief Creates a class name identifier.
     * @param id The class name.
     * @return Cyan-colored, bold identifier for the class.
     *
     * Creates an identifier formatted for class names using cyan color
     * and bold text.
     */
    inline
    IdentifierT<CyanT> Class(const std::string& id)
    {
      return id;
    }

    /**
     * @brief Creates a namespace identifier.
     * @param id The namespace name.
     * @return Cyan-colored, bold identifier for the namespace.
     *
     * Creates an identifier formatted for namespace names using cyan color
     * and bold text.
     */
    inline
    IdentifierT<CyanT> Namespace(const std::string& id)
    {
      return id;
    }

    /**
     * @brief Creates a function name identifier.
     * @param id The function name.
     * @return Default-colored, bold identifier for the function.
     *
     * Creates an identifier formatted for function names using the default
     * terminal color and bold text.
     */
    inline
    IdentifierT<ResetT> Function(const std::string& id)
    {
      return id;
    }
  }
}

#endif

