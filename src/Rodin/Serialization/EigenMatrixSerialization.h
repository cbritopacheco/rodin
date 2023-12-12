#ifndef RODIN_SERIALIZATION_MATRIXSERIALIZATION_H
#define RODIN_SERIALIZATION_MATRIXSERIALIZATION_H

#include <Eigen/Core>

namespace boost::serialization
{
  template <class Archive, typename Derived, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
  inline
  void serialize(
      Archive & ar,
      Eigen::Matrix<Derived, Rows, Cols, Options, MaxRows, MaxCols> & matrix,
      const unsigned int version)
  {
    const Eigen::Index rows = matrix.rows();
    const Eigen::Index cols = matrix.cols();
    ar & rows;
    ar & cols;
    if (Archive::is_saving::value)
    {
      for (Eigen::Index i = 0; i < rows; ++i)
        for (Eigen::Index j = 0; j < cols; ++j)
          ar & matrix(i, j);
    }
    else
    {
      matrix.resize(rows, cols);
      for (Eigen::Index i = 0; i < rows; ++i)
        for (Eigen::Index j = 0; j < cols; ++j)
          ar & matrix(i, j);
    }
  }
}

#endif
