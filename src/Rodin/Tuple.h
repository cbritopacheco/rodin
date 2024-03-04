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
  template <>
  class Tuple<> : public std::tuple<>
  {
    public:
      template <std::size_t Index>
      using Type = void;

      using Parent = std::tuple<>;

      static constexpr size_t Size = 0;

      constexpr
      Tuple() = default;

      constexpr
      Tuple(const Tuple&) = default;

      constexpr
      Tuple(Tuple&&) = default;

      inline
      constexpr
      Tuple<> product(const Tuple<> other)
      {
        return Tuple<>{};
      }

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

      using Parent::Parent;

      Tuple(const Tuple& other)
        : Parent(other)
      {}

      Tuple(Tuple&& other)
        : Parent(std::move(other))
      {}

      static constexpr size_t Size = sizeof...(Ts) + 1;

      inline
      constexpr
      Tuple& operator=(const Tuple& other)
      {
        if (this != &other)
          copyImpl<0>(other);
        return *this;
      }

      inline
      constexpr
      Tuple& operator=(Tuple&& other)
      {
        moveImpl<0>(std::move(other));
        return *this;
      }

      inline
      constexpr
      size_t size() const
      {
        return sizeof...(Ts) + 1;
      }

      template <class Function>
      inline
      constexpr
      auto reduce(Function&& func) const
      {
        static_assert(sizeof...(Ts) >= 1);
        return reduceImpl<0>(func);
      }

      inline
      constexpr
      Tuple<> product(const Tuple<>& other) const
      {
        return Tuple<>{};
      }

      template <class G, class ... Gs>
      inline
      constexpr
      auto product(const Tuple<G, Gs...>& other) const
      {
        return product([](const auto& a, const auto& b) { return Pair(a, b); }, other);
      }

      template <class Function, class G, class ... Gs>
      inline
      constexpr
      auto product(Function&& func, const Tuple<G, Gs...>& other) const
      {
        return productImpl<0, 0, Function, G, Gs...>(std::forward<Function>(func), other);
      }

      template <class ... Gs>
      inline
      constexpr
      auto zip(const Tuple<Gs...>& other) const
      {
        return zip([](const auto& a, const auto& b) { return Pair(a, b); }, other);
      }

      template <class Function, class ... Gs, class ... Os>
      inline
      constexpr
      auto zip(Function&& func, const Tuple<Gs...>& t, const Os&... other) const
      {
        return zipImpl<0, Function, Tuple<Gs...>, Os...>(
            std::forward<Function>(func), t, other...);
      }

      template <class Function>
      inline
      constexpr
      void apply(Function&& func)
      {
        applyImpl(std::forward<Function>(func), std::index_sequence_for<T, Ts...>{});
      }

      template <class Function>
      inline
      constexpr
      void iapply(Function&& func)
      {
        iapplyImpl(std::forward<Function>(func), std::index_sequence_for<T, Ts...>{});
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

      template <typename Function>
      inline
      constexpr
      auto map(Function&& func) const
      {
        return mapImpl(std::index_sequence_for<T, Ts...>(), std::forward<Function>(func));
      }

      template <typename Function>
      inline
      constexpr
      auto map(Function&& func)
      {
        return mapImpl(std::index_sequence_for<T, Ts...>(), std::forward<Function>(func));
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

    private:
      template <class Function, std::size_t ... Indices>
      inline
      constexpr
      void applyImpl(Function&& func, std::index_sequence<Indices...>)
      {
        (func(get<Indices>()), ...);
      }

      template <class Function, std::size_t ... Indices>
      inline
      constexpr
      void iapplyImpl(Function&& func, std::index_sequence<Indices...>)
      {
        (func(Indices, get<Indices>()), ...);
      }

      template <std::size_t ... Is, typename Func>
      inline
      constexpr
      auto mapImpl(std::index_sequence<Is...>, Func&& func) const
      {
        return Utility::Make<Tuple<decltype(func(get<Is>()))...>>()(
            func(std::get<Is>(*this))...);
      }

      template <std::size_t ... Is, typename Func>
      inline
      constexpr
      auto mapImpl(std::index_sequence<Is...>, Func&& func)
      {
        return Utility::Make<Tuple<decltype(func(get<Is>()))...>>()(
            func(std::get<Is>(*this))...);
      }

      template <typename ... Gs, std::size_t... Indices1, std::size_t... Indices2>
      inline
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
          return Utility::Make<Tuple<Type<Index>>>()(
              get<Index>()).concatenate(filterImpl<Index + 1, Predicate>());
        }
        else
        {
          return filterImpl<Index + 1, Predicate>();
        }
      }

      template <std::size_t Index, class Function, class ... Os>
      inline
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
      inline
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
      inline
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
      inline
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
      inline
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

  template <Index First, Index Last>
  inline
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


