/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_TUPLE_H
#define RODIN_TUPLE_H

#include <tuple>
#include "Rodin/Utility/ParameterPack.h"

namespace Rodin
{
  template <class ... Ts>
  class Tuple;

  template <>
  class Tuple<> : public std::tuple<>
  {
    public:
      template <std::size_t Index>
      using Type = void;

      using Parent = std::tuple<>;

      constexpr
      Tuple() = default;

      constexpr
      Tuple(const Tuple&) = default;

      constexpr
      Tuple(Tuple&&) = default;

      template <typename ... Gs>
      inline
      constexpr
      Tuple<Gs...> concatenate(const Tuple<Gs...>& other) const
      {
        return other;
      }

      template <template <class> class Predicate>
      inline
      constexpr
      Tuple<> filter() const
      {
        return Tuple{};
      }

      template <typename Function>
      inline
      constexpr
      void apply(Function&& func)
      {}

      inline
      constexpr
      size_t size() const
      {
        return 0;
      }
  };

  Tuple() -> Tuple<>;

  template <class T, class ... Ts>
  class Tuple<T, Ts...> : public std::tuple<T, Ts...>
  {
    public:
      template <std::size_t Index>
      using Type = typename Utility::ParameterPack<T, Ts...>::template At<Index>;

      using Parent = std::tuple<T, Ts...>;

      template <class ... Params>
      Tuple(Params&& ... params)
        : Parent(std::forward<Params>(params)...)
      {}

      Tuple(const Tuple& other)
        : Parent(other)
      {}

      Tuple(Tuple&& other)
        : Parent(std::move(other))
      {}

      template <typename Function>
      inline
      constexpr
      void apply(Function&& func)
      {
        applyImpl(std::forward<Function>(func), std::index_sequence_for<T, Ts...>{});
      }

      template <std::size_t Index>
      inline
      constexpr
      auto& get()
      {
        return std::get<Index>(*this);
      }

      template <std::size_t Index>
      inline
      constexpr
      const auto& get() const
      {
        return std::get<Index>(*this);
      }

      template <template <class> class Predicate>
      inline
      constexpr
      auto filter() const
      {
        return filterImpl<0, Predicate>();
      }

      template <typename Func>
      inline
      constexpr
      auto map(Func&& func) const
      {
        return mapImpl(std::index_sequence_for<T, Ts...>(), std::forward<Func>(func));
      }

      inline
      constexpr
      Tuple concatenate(const Tuple<>& other) const
      {
        return *this;
      }

      template <typename ... Gs>
      inline
      constexpr
      Tuple<T, Ts..., Gs...> concatenate(const Tuple<Gs...>& other) const
      {
        return concatenateImpl(other,
            std::index_sequence_for<T, Ts...>(), std::index_sequence_for<Gs...>());
      }

      inline
      constexpr
      size_t size() const
      {
        return sizeof...(Ts) + 1;
      }

    private:
      template <class Function, std::size_t ... Indices>
      inline
      constexpr
      void applyImpl(Function&& func, std::index_sequence<Indices...>)
      {
        (func(get<Indices>()), ...);
      }

      template <size_t I0, std::size_t ... Is, typename Func>
      inline
      constexpr
      auto mapImpl(std::index_sequence<I0, Is...>, Func&& func) const
      {
        return Tuple<
          std::invoke_result_t<Func, T>,
          std::invoke_result_t<Func, Ts>...>(func(get<I0>()), func(get<Is>())...);
      }

      template <typename ... Gs, std::size_t... Indices1, std::size_t... Indices2>
      inline
      constexpr
      Tuple<T, Ts..., Gs...> concatenateImpl(
          const Tuple<Gs...>& other,
          std::index_sequence<Indices1...>,
          std::index_sequence<Indices2...>) const
      {
        return Tuple<T, Ts..., Gs...>{get<Indices1>()..., other.template get<Indices2>()...};
      }

      template <std::size_t Index, template <class> class Predicate>
      inline
      constexpr
      auto filterImpl() const
      {
        if constexpr(Index == sizeof...(Ts) + 1)
        {
          return Tuple<>{};
        }
        else if constexpr (Predicate<Type<Index>>::Value)
        {
          return Tuple<Type<Index>>{ get<Index>() }.concatenate(filterImpl<Index + 1, Predicate>());
        }
        else
        {
          return filterImpl<Index + 1, Predicate>();
        }
      }
  };

  template <class ... Params>
  Tuple(Params&&...) -> Tuple<Params...>;
}

#endif


