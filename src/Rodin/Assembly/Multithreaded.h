/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ASSEMBLY_Multithreaded_H
#define RODIN_ASSEMBLY_Multithreaded_H

#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"

#include "ForwardDecls.h"
#include "AssemblyBase.h"

namespace Rodin::Assembly
{
  template <>
  class Multithreaded<Variational::BilinearFormBase<std::vector<Eigen::Triplet<Scalar>>>>
    : public AssemblyBase<Variational::BilinearFormBase<std::vector<Eigen::Triplet<Scalar>>>>
  {
    /**
     * @internal
     */
    static void add(
        std::vector<Eigen::Triplet<Scalar>>& out, const Math::Matrix& in,
        const IndexArray& rows, const IndexArray& cols);

    public:
      using Parent =
        AssemblyBase<Variational::BilinearFormBase<std::vector<Eigen::Triplet<Scalar>>>>;
      using OperatorType = std::vector<Eigen::Triplet<Scalar>>;

      Multithreaded() = default;

      Multithreaded(const Multithreaded& other)
        : Parent(other)
      {}

      Multithreaded(Multithreaded&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Executes the assembly and returns the linear operator
       * associated to the bilinear form.
       */
      OperatorType execute(const BilinearAssemblyInput& input) const override;

      Multithreaded* copy() const noexcept override
      {
        return new Multithreaded(*this);
      }
  };

  /**
   * @brief Multithreaded assembly of the Math::SparseMatrix associated to a
   * BilinearFormBase object.
   */
  template <>
  class Multithreaded<Variational::BilinearFormBase<Math::SparseMatrix>>
    : public AssemblyBase<Variational::BilinearFormBase<Math::SparseMatrix>>
  {
    public:
      using Parent = AssemblyBase<Variational::BilinearFormBase<Math::SparseMatrix>>;
      using OperatorType = Math::SparseMatrix;

      Multithreaded() = default;

      Multithreaded(const Multithreaded& other)
        : Parent(other)
      {}

      Multithreaded(Multithreaded&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Executes the assembly and returns the linear operator
       * associated to the bilinear form.
       */
      OperatorType execute(const BilinearAssemblyInput& input) const override;

      Multithreaded* copy() const noexcept override
      {
        return new Multithreaded(*this);
      }
  };

  /**
   * @brief %Multithreaded assembly of the Math::Vector associated to a LinearFormBase
   * object.
   */
  template <>
  class Multithreaded<Variational::LinearFormBase<Math::Vector>>
    : public AssemblyBase<Variational::LinearFormBase<Math::Vector>>
  {

    static void add(Math::Vector& out, const Math::Vector& in, const IndexArray& s);

    public:
      using Parent = AssemblyBase<Variational::LinearFormBase<Math::Vector>>;
      using VectorType = Math::Vector;

      Multithreaded() = default;

      Multithreaded(const Multithreaded& other)
        : Parent(other)
      {}

      Multithreaded(Multithreaded&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Executes the assembly and returns the vector associated to the
       * linear form.
       */
      VectorType execute(const Input& input) const override;

      Multithreaded* copy() const noexcept override
      {
        return new Multithreaded(*this);
      }
  };
}

#endif

