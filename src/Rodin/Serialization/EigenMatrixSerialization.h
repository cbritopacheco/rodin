#ifndef RODIN_SERIALIZATION_MATRIXSERIALIZATION_H
#define RODIN_SERIALIZATION_MATRIXSERIALIZATION_H

#include <Eigen/Core>

#include <boost/serialization/array.hpp>

namespace boost::serialization
{
  template <class Archive, typename Derived, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
  inline
  void serialize(
      Archive & ar,
      Eigen::Matrix<Derived, Rows, Cols, Options, MaxRows, MaxCols>& matrix,
      const unsigned int version)
  {
    const Eigen::Index rows = matrix.rows();
    const Eigen::Index cols = matrix.cols();
    ar & rows;
    ar & cols;
    ar & boost::serialization::make_array(matrix.data(), matrix.size());
  }
}

#endif
