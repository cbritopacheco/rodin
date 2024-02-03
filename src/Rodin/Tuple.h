/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_TUPLE_H
#define RODIN_TUPLE_H

#include <tuple>
#include "ForwardDecls.h"
#include "Rodin/Utility/ParameterPack.h"
#include "Rodin/Utility/Make.h"
#include "Rodin/Utility/UnwrapReference.h"

namespace Rodin
{
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

      template <class ... Params>
      Tuple(Params&& ... params)
        : Parent(std::forward<Params>(params)...)
      {}

      constexpr
      Tuple(const Tuple& other)
        : Parent(other)
      {}

      constexpr
      Tuple(Tuple&& other)
        : Parent(std::move(other))
      {}

      constexpr
      Tuple& operator=(const Tuple& other)
      {
        return Parent::operator=(other);
      }

      constexpr
      Tuple& operator=(Tuple&& other)
      {
        return Parent::operator=(std::move(other));
      }

      inline
      constexpr
      Tuple<> product(const Tuple<> other)
      {
        return Tuple<>{};
      }

      template <template <class, class> class Pair, class G, class ... Gs>
      inline
      constexpr
      auto product(const Tuple<G, Gs...> other) const
      {
        return productImpl<0, 0, Pair, G, Gs...>(other);
      }

      template <template <class, class> class Pair, class G, class ... Gs>
      inline
      constexpr
      Tuple<Pair<T, G>, Pair<Ts, Gs>...> zip(const Tuple<G, Gs...> other) const
      {
        static_assert(sizeof...(Ts) == sizeof...(Gs));
        return zipImpl<0, Pair, G, Gs...>(other);
      }

      template <template <class> class External>
      inline
      constexpr
      auto wrap() const
      {
        return wrapImpl<0, External>();
      }

      template <class Function>
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
      template <template <class, class> class Pair, class ... Gs, std::size_t ... Indices>
      inline
      constexpr
      void zipImpl(const Tuple<Gs...>& other, std::index_sequence<Indices...>) const
      {
        return Tuple(
            Pair<
                decltype(get<Indices>()), decltype(other.template get<Indices>())>
                (get<Indices>(), other.template get<Indices>())...);
      }

      template <class Function, std::size_t ... Indices>
      inline
      constexpr
      void applyImpl(Function&& func, std::index_sequence<Indices...>)
      {
        (func(get<Indices>()), ...);
      }

      template <std::size_t ... Is, typename Func>
      inline
      constexpr
      auto mapImpl(std::index_sequence<Is...>, Func&& func) const
      {
        using Result =
          Tuple<typename Utility::UnwrapRefDecay<decltype(func(std::declval<Type<Is>&>()))>::Type...>;
        return Utility::Make<Result>()(func(std::get<Is>(*this))...);
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

      template <std::size_t Index, template <class, class> class Pair, class G, class ... Gs>
      inline
      constexpr
      auto zipImpl(const Tuple<G, Gs...>& other) const
      {
        using OtherTuple = Tuple<G, Gs...>;
        if constexpr (Index == sizeof...(Ts) + 1)
        {
          return Tuple<>{};
        }
        else
        {
          return
            Tuple<
              Pair<Type<Index>, typename OtherTuple::template Type<Index>>>(
                Pair(get<Index>(), other.template get<Index>()))
              .concatenate(zipImpl<Index + 1, Pair, G, Gs...>(other));
        }
      }

      template <std::size_t Index, template <class> class External>
      inline
      constexpr
      auto wrapImpl() const
      {
        if constexpr (Index == sizeof...(Ts) + 1)
        {
          return Tuple<>{};
        }
        else
        {
          return Tuple<Type<Index>>(External<Type<Index>>(get<Index>())
              ).concatenate(wrapImpl<Index + 1, External>());
        }
      }

      template <
        std::size_t Index,
        std::size_t OtherIndex,
        template <class, class> class Pair, class G, class ... Gs>
      inline
      constexpr
      auto productImpl(const Tuple<G, Gs...>& other) const
      {
        static_assert(Index <= sizeof...(Ts) && OtherIndex <= sizeof...(Gs));
        using Type = Type<Index>;
        using OtherTuple = Tuple<G, Gs...>;
        using OtherType = typename OtherTuple::template Type<OtherIndex>;
        using PairType = Pair<Type, OtherType>;
        if constexpr (Index == sizeof...(Ts) && OtherIndex == sizeof...(Gs))
        {
          return Tuple<PairType>(PairType(get<Index>(), other.template get<OtherIndex>()));
        }
        else if constexpr (Index < sizeof...(Ts) && OtherIndex == sizeof...(Gs))
        {
          return Tuple<PairType>(PairType(get<Index>(), other.template get<OtherIndex>()))
            .concatenate(productImpl<Index + 1, 0, Pair, G, Gs...>(other));
        }
        else
        {
          static_assert(Index <= sizeof...(Ts) && OtherIndex < sizeof...(Gs));
          return Tuple<PairType>(PairType(get<Index>(), other.template get<OtherIndex>()))
            .concatenate(productImpl<Index, OtherIndex + 1, Pair, G, Gs...>(other));
        }
      }
  };

  template <class ... Params>
  Tuple(Params&&...) -> Tuple<Params...>;
}

#endif


