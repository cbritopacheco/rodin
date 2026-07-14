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
  /**
   * @brief Saves an Eigen matrix to an archive.
   * @tparam Archive Archive type.
   * @tparam Scalar Matrix scalar type.
   * @tparam Rows Compile-time row count.
   * @tparam Cols Compile-time column count.
   * @tparam Options Eigen storage options.
   * @tparam MaxRows Compile-time maximum row count.
   * @tparam MaxCols Compile-time maximum column count.
   * @param ar Archive used for serialization.
   * @param matrix Matrix to save.
   * @param version Serialization version.
   */
  template <class Archive, typename Scalar, int Rows, int Cols, int Options, int MaxRows,
    int MaxCols>
  void save(Archive& ar,
    const Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>& matrix,
    const unsigned int version)
  {
    (void)version;
    Eigen::Index rows = matrix.rows();
    Eigen::Index cols = matrix.cols();
    ar & rows;
    ar & cols;
    ar & boost::serialization::make_array(matrix.data(), matrix.size());
  }

  /**
   * @brief Loads an Eigen matrix from an archive.
   * @tparam Archive Archive type.
   * @tparam Scalar Matrix scalar type.
   * @tparam Rows Compile-time row count.
   * @tparam Cols Compile-time column count.
   * @tparam Options Eigen storage options.
   * @tparam MaxRows Compile-time maximum row count.
   * @tparam MaxCols Compile-time maximum column count.
   * @param ar Archive used for serialization.
   * @param matrix Matrix to load.
   * @param version Serialization version.
   */
  template <class Archive, typename Scalar, int Rows, int Cols, int Options, int MaxRows,
    int MaxCols>
  void load(Archive& ar,
    Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>& matrix,
    const unsigned int version)
  {
    (void)version;
    Eigen::Index rows, cols;
    ar & rows;
    ar & cols;
    matrix.resize(rows, cols);
    ar & boost::serialization::make_array(matrix.data(), rows * cols);
  }

  /**
   * @brief Serializes an Eigen matrix through split save/load functions.
   * @tparam Archive Archive type.
   * @tparam Scalar Matrix scalar type.
   * @tparam Rows Compile-time row count.
   * @tparam Cols Compile-time column count.
   * @tparam Options Eigen storage options.
   * @tparam MaxRows Compile-time maximum row count.
   * @tparam MaxCols Compile-time maximum column count.
   * @param ar Archive used for serialization.
   * @param matrix Matrix to serialize.
   * @param version Serialization version.
   */
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
