/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file EigenMatrix.h
 * @brief boost::serialization support for Eigen dense matrix types.
 */
#ifndef RODIN_SERIALIZATION_EIGENMATRIX_H
#define RODIN_SERIALIZATION_EIGENMATRIX_H

#include <Eigen/Core>

#include <boost/serialization/array.hpp>

namespace boost::serialization
{
  template <class Archive, typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
  void save(Archive & ar,
            const Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>& matrix,
            const unsigned int)
  {
    Eigen::Index rows = matrix.rows();
    Eigen::Index cols = matrix.cols();
    ar & rows;
    ar & cols;
    ar & boost::serialization::make_array(matrix.data(), matrix.size());
  }

  template <class Archive, typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
  void load(Archive & ar,
            Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>& matrix,
            const unsigned int)
  {
    Eigen::Index rows, cols;
    ar & rows;
    ar & cols;
    matrix.resize(rows, cols);
    ar & boost::serialization::make_array(matrix.data(), rows * cols);
  }

  template <class Archive, typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
  void serialize(
      Archive & ar,
      Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>& matrix,
      const unsigned int version)
  {
    boost::serialization::split_free(ar, matrix, version);
  }
}

#endif
