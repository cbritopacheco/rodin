#include <vector>
#include <boost/range/adaptors.hpp>
#include <boost/unordered_map.hpp>
#include <boost/container/flat_set.hpp>
#include <boost/container/flat_map.hpp>

#include "ForwardDecls.h"

namespace Rodin::Geometry
{
  using Index = std::size_t;

  using Attribute = std::size_t;

  template <class T>
  using FlatSet = boost::container::flat_set<T>;

  template <class K, class T>
  using Map = boost::container::map<K, T>;

  template <class K, class T>
  using UnorderedMap = boost::unordered_map<K, T>;

  template <class K, class T>
  using FlatMap = boost::container::flat_map<K, T>;

  using IndexSet = FlatSet<Index>;

  using Incidence = std::vector<IndexSet>;
}
