#ifndef RODIN_SERIALIZATION_OPTIONAL_H
#define RODIN_SERIALIZATION_OPTIONAL_H

#include <optional>
#include <boost/serialization/version.hpp>
#include <boost/serialization/split_free.hpp>

namespace boost::serialization
{
  template <class Archive, class T>
  void save(Archive& ar, const std::optional<T>& opt, const unsigned int)
  {
      bool has_value = opt.has_value();
      ar & has_value;

      if (has_value)
          ar & *opt;
  }

  template <class Archive, class T>
  void load(Archive& ar, std::optional<T>& opt, const unsigned int)
  {
      bool has_value;
      ar & has_value;

      if (has_value)
      {
          T value;
          ar & value;
          opt = std::move(value);
      }
      else
      {
          opt.reset();
      }
  }

  template <class Archive, class T>
  void serialize(Archive& ar, std::optional<T>& opt, const unsigned int version)
  {
      split_free(ar, opt, version);
  }
}

#endif
