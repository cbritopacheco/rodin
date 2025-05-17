#ifndef RODIN_SERIALIZATION_FLATSET_H
#define RODIN_SERIALIZATION_FLATSET_H

#include <boost/config.hpp>
#include <boost/version.hpp>
#include <boost/container/flat_set.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/nvp.hpp>

namespace boost::serialization
{
  template <class Archive, typename T, typename Compare, typename Allocator>
  void save(Archive & ar,
            const boost::container::flat_set<T, Compare, Allocator>& s,
            const unsigned int /*version*/)
  {
     std::size_t sz = s.size();
     ar & BOOST_SERIALIZATION_NVP(sz);
#if BOOST_VERSION >= 108100
     if(sz > 0)
     {
        const auto & seq = s.get_sequence_cref();
        ar.save_binary(reinterpret_cast<const char*>(seq.data()), sz * sizeof(T));
     }
#else
     for(const auto& elem : s)
        ar & BOOST_SERIALIZATION_NVP(elem);
#endif
  }

  template <class Archive, typename T, typename Compare, typename Allocator>
  void load(Archive & ar,
            boost::container::flat_set<T, Compare, Allocator>& s,
            const unsigned int /*version*/)
  {
     std::size_t sz;
     ar & BOOST_SERIALIZATION_NVP(sz);
     s.clear();
#if BOOST_VERSION >= 108100
     if(sz > 0)
     {
        auto & seq = s.get_sequence_ref();
        seq.resize(sz);
        ar.load_binary(reinterpret_cast<char*>(seq.data()),
                       sz * sizeof(T));
     }
#else
     for(std::size_t i = 0; i < sz; ++i)
     {
        T elem;
        ar & BOOST_SERIALIZATION_NVP(elem);
        s.insert(elem);
     }
#endif
  }

  template <class Archive, typename T, typename Compare, typename Allocator>
  void serialize(Archive & ar,
                 boost::container::flat_set<T, Compare, Allocator>& s,
                 const unsigned int version)
  {
    boost::serialization::split_free(ar, s, version);
  }
}

#endif

