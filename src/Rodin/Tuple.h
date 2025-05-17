/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_TUPLE_H
#define RODIN_TUPLE_H

#include <tuple>

#include "Rodin/Types.h"
#include "Rodin/Utility/Make.h"
#include "Rodin/Utility/ParameterPack.h"
#include "Rodin/Utility/UnwrapReference.h"

#include "ForwardDecls.h"

namespace Rodin
{
  /**
   * @brief Specialization of Tuple for an empty parameter pack.
   */
  template <>
  class Tuple<> : public std::tuple<>
  {
    public:
      /**
       * @brief Dummy type alias for the empty tuple.
       *
       * @tparam Index Unused index parameter.
       */
      template <std::size_t Index>
      using Type = void;

      /// @brief Parent class type.
      using Parent = std::tuple<>;

      /// @brief The number of elements in the tuple (always 0).
      static constexpr size_t Size = 0;

      /**
       * @brief Default constructor.
       */
      constexpr
      Tuple() = default;

      /**
       * @brief Copy constructor.
       */
      constexpr
      Tuple(const Tuple&) = default;

      /**
       * @brief Move constructor.
       */
      constexpr
      Tuple(Tuple&&) = default;

      /**
       * @brief Computes the Cartesian product of two empty tuples.
       *
       * @param other The other empty tuple.
       * @return An empty Tuple.
       */
      constexpr
      Tuple<> product(const Tuple<> other)
      {
        return Tuple<>{};
      }

      /**
       * @brief Concatenates this empty tuple with another tuple.
       *
       * @tparam Gs Types of elements in the other tuple.
       * @param other The tuple to concatenate.
       * @return A new tuple containing the elements of the other tuple.
       */
      template <typename ... Gs>
      constexpr
      Tuple<Gs...> concatenate(const Tuple<Gs...>& other) const
      {
        return other;
      }

      /**
       * @brief Filters the empty tuple based on a predicate.
       *
       * @tparam Predicate A template predicate class.
       * @return An empty Tuple.
       */
      template <template <class> class Predicate>
      constexpr
      Tuple<> filter() const
      {
        return Tuple{};
      }

      /**
       * @brief Applies a function to each element of the tuple.
       *
       * @tparam Function The type of the function to apply.
       * @param func The function to apply.
       */
      template <typename Function>
      constexpr
      void apply(Function&& func)
      {}

      /**
       * @brief Returns the size of the tuple.
       *
       * @return Always returns 0.
       */
      constexpr
      size_t size() const
      {
        return 0;
      }
  };

  Tuple() -> Tuple<>;

  /**
   * @brief A tuple class that extends std::tuple with additional functionality.
   *
   * @tparam T The type of the first element.
   * @tparam Ts Types of the remaining elements.
   */
  template <class T, class ... Ts>
  class Tuple<T, Ts...> : public std::tuple<T, Ts...>
  {
    public:
      /**
       * @brief Type alias for the element at a given index.
       *
       * @tparam Index The index of the element.
       */
      template <std::size_t Index>
      using Type = typename Utility::ParameterPack<T, Ts...>::template At<Index>;

      /// @brief Parent class type.
      using Parent = std::tuple<T, Ts...>;

      using Parent::Parent;

      /**
       * @brief Copy constructor.
       *
       * @param other The tuple to copy.
       */
      Tuple(const Tuple& other)
        : Parent(other)
      {}

      /**
        * @brief Move constructor.
        *
        * @param other The tuple to move.
        */
      Tuple(Tuple&& other)
        : Parent(std::move(other))
      {}

      /// @brief The number of elements in the tuple.
      static constexpr size_t Size = sizeof...(Ts) + 1;

      /**
       * @brief Copy assignment operator.
       *
       * @param other The tuple to copy from.
       * @return Reference to this tuple.
       */
      constexpr
      Tuple& operator=(const Tuple& other)
      {
        if (this != &other)
          copyImpl<0>(other);
        return *this;
      }

      /**
       * @brief Move assignment operator.
       *
       * @param other The tuple to move from.
       * @return Reference to this tuple.
       */
      constexpr
      Tuple& operator=(Tuple&& other)
      {
        moveImpl<0>(std::move(other));
        return *this;
      }

      /**
       * @brief Returns the size of the tuple.
       *
       * @return The number of elements in the tuple.
       */
      constexpr
      size_t size() const
      {
        return sizeof...(Ts) + 1;
      }

      /**
       * @brief Reduces the tuple elements using a binary function.
       *
       * @tparam Function The type of the binary function.
       * @param func The function to apply.
       * @return The result of the reduction.
       *
       * @note Requires the tuple to have at least two elements.
       */
      template <class Function>
      constexpr
      auto reduce(Function&& func) const
      {
        static_assert(sizeof...(Ts) >= 1);
        return reduceImpl<0>(func);
      }

      /**
       * @brief Computes the Cartesian product with an empty tuple.
       *
       * @param other An empty tuple.
       * @return An empty Tuple.
       */
      constexpr
      Tuple<> product(const Tuple<>& other) const
      {
        return Tuple<>{};
      }

      /**
       * @brief Computes the Cartesian product with another tuple.
       *
       * @tparam G Type of the first element of the other tuple.
       * @tparam Gs Types of the remaining elements of the other tuple.
       * @param other The other tuple.
       * @return A tuple representing the Cartesian product.
       */
      template <class G, class ... Gs>
      constexpr
      auto product(const Tuple<G, Gs...>& other) const
      {
        return product([](const auto& a, const auto& b) { return Pair(a, b); }, other);
      }

      /**
       * @brief Computes the Cartesian product with another tuple using a custom function.
       *
       * @tparam Function The type of the function to combine elements.
       * @tparam G Type of the first element of the other tuple.
       * @tparam Gs Types of the remaining elements of the other tuple.
       * @param func The function to combine elements.
       * @param other The other tuple.
       * @return A tuple representing the Cartesian product.
       */
      template <class Function, class G, class ... Gs>
      constexpr
      auto product(Function&& func, const Tuple<G, Gs...>& other) const
      {
        return productImpl<0, 0, Function, G, Gs...>(std::forward<Function>(func), other);
      }

      /**
       * @brief Zips this tuple with another tuple using a default pairing function.
       *
       * @tparam Gs Types of the elements in the other tuple.
       * @param other The other tuple.
       * @return A tuple containing pairs of elements.
       */
      template <class ... Gs>
      constexpr
      auto zip(const Tuple<Gs...>& other) const
      {
        return zip([](const auto& a, const auto& b) { return Pair(a, b); }, other);
      }

      /**
       * @brief Zips this tuple with one or more other tuples using a custom function.
       *
       * @tparam Function The type of the function to combine elements.
       * @tparam Gs Types of the elements in the first other tuple.
       * @tparam Os Types of the elements in the additional tuples.
       * @param func The function to combine elements.
       * @param t The first other tuple.
       * @param other Additional tuples to zip.
       * @return A tuple containing the zipped result.
       */
      template <class Function, class ... Gs, class ... Os>
      constexpr
      auto zip(Function&& func, const Tuple<Gs...>& t, const Os&... other) const
      {
        return zipImpl<0, Function, Tuple<Gs...>, Os...>(
            std::forward<Function>(func), t, other...);
      }

      /**
       * @brief Applies a function to each element of the tuple.
       *
       * @tparam Function The type of the function.
       * @param func The function to apply.
       * @return Reference to this tuple.
       */
      template <class Function>
      constexpr
      Tuple& apply(Function&& func)
      {
        applyImpl(std::forward<Function>(func), std::index_sequence_for<T, Ts...>{});
        return *this;
      }

      /**
       * @brief Applies a function to each element of the tuple (const version).
       *
       * @tparam Function The type of the function.
       * @param func The function to apply.
       * @return Const reference to this tuple.
       */
      template <class Function>
      constexpr
      const Tuple& apply(Function&& func) const
      {
        applyImpl(std::forward<Function>(func), std::index_sequence_for<T, Ts...>{});
        return *this;
      }

      /**
       * @brief Applies a function to each element of the tuple along with its index.
       *
       * @tparam Function The type of the function.
       * @param func The function to apply. It receives the index and the element.
       * @return Reference to this tuple.
       */
      template <class Function>
      constexpr
      Tuple& iapply(Function&& func)
      {
        iapplyImpl(std::forward<Function>(func), std::index_sequence_for<T, Ts...>{});
        return *this;
      }

      /**
       * @brief Returns a reference to the element at the given index.
       *
       * @tparam Index The index of the element.
       * @return Reference to the element.
       */
      template <std::size_t Index>
      constexpr
      auto& get()
      {
        return std::get<Index>(*this);
      }

      /**
       * @brief Returns a const reference to the element at the given index.
       *
       * @tparam Index The index of the element.
       * @return Const reference to the element.
       */
      template <std::size_t Index>
      constexpr
      const auto& get() const
      {
        return std::get<Index>(*this);
      }

      /**
       * @brief Filters the tuple based on a predicate.
       *
       * @tparam Predicate A template predicate class.
       * @return A new tuple containing only the elements that satisfy the predicate.
       */
      template <template <class> class Predicate>
      constexpr
      auto filter() const
      {
        return filterImpl<0, Predicate>();
      }

      /**
       * @brief Transforms each element of the tuple using a function (const version).
       *
       * @tparam Function The type of the function.
       * @param func The transformation function.
       * @return A new tuple with transformed elements.
       */
      template <typename Function>
      constexpr
      auto map(Function&& func) const
      {
        return mapImpl(std::index_sequence_for<T, Ts...>(), std::forward<Function>(func));
      }

      /**
       * @brief Transforms each element of the tuple using a function.
       *
       * @tparam Function The type of the function.
       * @param func The transformation function.
       * @return A new tuple with transformed elements.
       */
      template <typename Function>
      constexpr
      auto map(Function&& func)
      {
        return mapImpl(std::index_sequence_for<T, Ts...>(), std::forward<Function>(func));
      }

      /**
       * @brief Concatenates this tuple with an empty tuple.
       *
       * @param other An empty tuple.
       * @return A new tuple identical to this tuple.
       */
      constexpr
      Tuple concatenate(const Tuple<>& other) const
      {
        return *this;
      }

      /**
       * @brief Concatenates this tuple with another tuple.
       *
       * @tparam Gs Types of the elements in the other tuple.
       * @param other The other tuple.
       * @return A new tuple containing the elements of this tuple followed by those of the other tuple.
       */
      template <typename ... Gs>
      constexpr
      Tuple<T, Ts..., Gs...> concatenate(const Tuple<Gs...>& other) const
      {
        return concatenateImpl(other,
            std::index_sequence_for<T, Ts...>(), std::index_sequence_for<Gs...>());
      }

    private:
      template <class Function, std::size_t ... Indices>
      constexpr
      void applyImpl(Function&& func, std::index_sequence<Indices...>)
      {
        (func(get<Indices>()), ...);
      }

      template <class Function, std::size_t ... Indices>
      constexpr
      void applyImpl(Function&& func, std::index_sequence<Indices...>) const
      {
        (func(get<Indices>()), ...);
      }

      template <class Function, std::size_t ... Indices>
      constexpr
      void iapplyImpl(Function&& func, std::index_sequence<Indices...>)
      {
        (func(Indices, get<Indices>()), ...);
      }

      template <std::size_t ... Is, typename Func>
      constexpr
      auto mapImpl(std::index_sequence<Is...>, Func&& func) const
      {
        return Utility::Make<Tuple<decltype(func(get<Is>()))...>>()(
            func(std::get<Is>(*this))...);
      }

      template <std::size_t ... Is, typename Func>
      constexpr
      auto mapImpl(std::index_sequence<Is...>, Func&& func)
      {
        return Utility::Make<Tuple<decltype(func(get<Is>()))...>>()(
            func(std::get<Is>(*this))...);
      }

      template <typename ... Gs, std::size_t... Indices1, std::size_t... Indices2>
      constexpr
      Tuple<T, Ts..., Gs...> concatenateImpl(
          const Tuple<Gs...>& other,
          std::index_sequence<Indices1...>,
          std::index_sequence<Indices2...>) const
      {
        return Utility::Make<Tuple<T, Ts..., Gs...>>()(
            get<Indices1>()..., other.template get<Indices2>()...);
      }

      template <std::size_t Index, template <class> class Predicate>
      constexpr
      auto filterImpl() const
      {
        if constexpr(Index == sizeof...(Ts) + 1)
        {
          return Tuple<>{};
        }
        else if constexpr (Predicate<Type<Index>>::Value)
        {
          return Utility::Make<Tuple<Type<Index>>>()(
              get<Index>()).concatenate(filterImpl<Index + 1, Predicate>());
        }
        else
        {
          return filterImpl<Index + 1, Predicate>();
        }
      }

      template <std::size_t Index, class Function, class ... Os>
      constexpr
      auto zipImpl(Function&& func, const Os&... other) const
      {
        if constexpr (Index == sizeof...(Ts) + 1)
        {
          return Tuple<>{};
        }
        else
        {
          using R = decltype(func(get<Index>(), other.template get<Index>()...));
          return Utility::Make<Tuple<R>>()(func(get<Index>(), other.template get<Index>()...)
              ).concatenate(zipImpl<Index + 1, Function, Os...>(std::forward<Function>(func), other...));
        }
      }

      template <std::size_t Index, std::size_t OtherIndex, class Function, class G, class ... Gs>
      constexpr
      auto productImpl(Function&& func, const Tuple<G, Gs...>& other) const
      {
        static_assert(Index <= sizeof...(Ts) && OtherIndex <= sizeof...(Gs));
        using R = decltype(func(get<Index>(), other.template get<OtherIndex>()));
        if constexpr (Index == sizeof...(Ts) && OtherIndex == sizeof...(Gs))
        {
          return Utility::Make<Tuple<R>>()(
              func(get<Index>(), other.template get<OtherIndex>()));
        }
        else if constexpr (Index < sizeof...(Ts) && OtherIndex == sizeof...(Gs))
        {
          return Utility::Make<Tuple<R>>()(
              func(get<Index>(), other.template get<OtherIndex>()))
            .concatenate(productImpl<Index + 1, 0, Function, G, Gs...>(std::forward<Function>(func), other));
        }
        else
        {
          static_assert(Index <= sizeof...(Ts) && OtherIndex < sizeof...(Gs));
          return Utility::Make<Tuple<R>>()(
              func(get<Index>(), other.template get<OtherIndex>()))
            .concatenate(productImpl<Index, OtherIndex + 1, Function, G, Gs...>(std::forward<Function>(func), other));
        }
      }

      template <std::size_t Index, typename Function>
      constexpr
      auto reduceImpl(Function&& func) const
      {
        if constexpr (Index == sizeof...(Ts) - 1)
        {
          return func(get<Index>(), get<Index + 1>());
        }
        else
        {
          return func(get<Index>(), reduceImpl<Index + 1, Function>(std::forward<Function>(func)));
        }
      }

      template <std::size_t Index>
      constexpr
      void copyImpl(const Tuple& other)
      {
        if constexpr (Index < Size)
        {
          get<Index>() = other.get<Index>();
          copyImpl<Index + 1>(other);
        }
      }

      template <std::size_t Index>
      constexpr
      void moveImpl(Tuple&& other)
      {
        if constexpr (Index < Size)
        {
          get<Index>() = std::move(other.get<Index>());
          moveImpl<Index + 1>(std::move(other));
        }
      }
  };

  template <class ... Params>
  Tuple(Params...) -> Tuple<Params...>;

  /**
   * @brief Generates a tuple of indices.
   *
   * @tparam First The starting index.
   * @tparam Last One past the last index.
   * @return A tuple containing indices from First to Last-1.
   *
   * @note First must be less than Last.
   */
  template <Index First, Index Last>
  static constexpr auto IndexTuple()
  {
    static_assert(First < Last);
    if constexpr (First == Last)
    {
      return Tuple<>{};
    }
    else if constexpr (First == Last - 1)
    {
      return Tuple<Index>{Last - 1};
    }
    else
    {
      return Tuple<Index>{First}.concatenate(IndexTuple<First + 1, Last>());
    }
  }
}

#include "Pair.h"

#endif


