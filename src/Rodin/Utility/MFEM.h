#ifndef RODIN_UTILITY_MFEM_H
#define RODIN_UTILITY_MFEM_H

#include <set>
#include <mfem.hpp>

#include "Rodin/Math/Vector.h"
#include "Rodin/Math/DenseMatrix.h"

#include "Wrap.h"

namespace Rodin::Utility
{
  inline mfem::Array<int> set2marker(const std::set<int>& s, int size)
  {
    mfem::Array<int> res(size);
    res = 0;
    for (const auto& v : s)
    {
      assert(v > 0);
      assert(v - 1 < size);
      res[v - 1] = 1;
    }
    return res;
  }

  template <>
  class Wrap<mfem::Vector&, Eigen::Map<Math::Vector>>
    : public Eigen::Map<Math::Vector>
  {
    public:
      inline
      Wrap(mfem::Vector& obj)
        : Eigen::Map<Math::Vector>(obj.GetData(), obj.Size())
      {}

      inline
      Wrap(const Wrap& other)
        : Eigen::Map<Math::Vector>(other)
      {}

      Wrap(Wrap&&) = delete;
      void operator=(const Wrap&) = delete;
      void operator=(Wrap&&) = delete;
  };

  template <>
  class Wrap<mfem::Vector&&, Eigen::Map<Math::Vector>>
    : public Eigen::Map<Math::Vector>
  {
    public:
      inline
      Wrap(mfem::Vector&& obj)
        : Eigen::Map<Math::Vector>(obj.GetData(), obj.Size())
      {
        m_obj.Swap(obj);
      }

      Wrap(const Wrap&) = delete;
      Wrap(Wrap&&) = delete;
      void operator=(const Wrap&) = delete;
      void operator=(Wrap&&) = delete;
    private:
      mfem::Vector m_obj;
  };

  template <>
  class Wrap<mfem::DenseMatrix&, Eigen::Map<Math::Matrix>>
    : public Eigen::Map<Math::Matrix>
  {
    public:
      inline
      Wrap(mfem::DenseMatrix& obj)
        : Eigen::Map<Math::Matrix>(obj.GetData(), obj.NumRows(), obj.NumCols())
      {}

      inline
      Wrap(const Wrap& other)
        : Eigen::Map<Math::Matrix>(other)
      {}

      Wrap(Wrap&&) = delete;
      void operator=(const Wrap&) = delete;
      void operator=(Wrap&&) = delete;
  };

  template <>
  class Wrap<mfem::DenseMatrix&&, Eigen::Map<Math::Matrix>>
    : public Eigen::Map<Math::Matrix>
  {
    public:
      inline
      Wrap(mfem::DenseMatrix&& obj)
        : Eigen::Map<Math::Matrix>(obj.GetData(), obj.NumRows(), obj.NumCols())
      {
        m_obj.Swap(obj);
      }

      Wrap(const Wrap&) = delete;
      Wrap(Wrap&&) = delete;
      void operator=(const Wrap&) = delete;
      void operator=(Wrap&&) = delete;

    private:
      mfem::DenseMatrix m_obj;
  };
}

namespace Rodin::Utility
{
  template <>
  class Wrap<Math::Vector&, mfem::Vector> : public mfem::Vector
  {
    public:
      inline
      Wrap(Math::Vector& obj)
        : m_obj(obj)
      {
        SetDataAndSize(m_obj.get().data(), m_obj.get().size());
      }

      inline
      Wrap(const Wrap& other)
        : m_obj(other.m_obj)
      {
        SetDataAndSize(m_obj.get().data(), m_obj.get().size());
      }

      Wrap(Wrap&&) = delete;
      void operator=(Wrap&&) = delete;
      void operator=(const Wrap&) = delete;

    private:
      std::reference_wrapper<Math::Vector> m_obj;
  };

  template <>
  class Wrap<Math::Vector&&, mfem::Vector> : public mfem::Vector
  {
    public:
      inline
      Wrap(Math::Vector&& obj)
        : m_obj(std::move(obj))
      {
        SetDataAndSize(m_obj.data(), m_obj.size());
      }

      Wrap(const Wrap& other) = delete;
      Wrap(Wrap&&) = delete;
      void operator=(const Wrap&) = delete;
      void operator=(Wrap&&) = delete;

    private:
      Math::Vector m_obj;
  };

  template <>
  class Wrap<Math::Matrix&, mfem::DenseMatrix> : public mfem::DenseMatrix
  {
    public:
      inline
      Wrap(Math::Matrix& obj)
        : m_obj(obj)
      {
        UseExternalData(m_obj.get().data(), m_obj.get().rows(), m_obj.get().cols());
      }

      inline
      Wrap(const Wrap& other)
        : m_obj(other.m_obj)
      {
        UseExternalData(m_obj.get().data(), m_obj.get().rows(), m_obj.get().cols());
      }

      Wrap(Wrap&&) = delete;
      void operator=(Wrap&&) = delete;
      void operator=(const Wrap&) = delete;

    private:
      std::reference_wrapper<Math::Matrix> m_obj;
  };

  template <>
  class Wrap<Math::Matrix&&, mfem::DenseMatrix> : public mfem::DenseMatrix
  {
    public:
      inline
      Wrap(Math::Matrix&& obj)
        : m_obj(std::move(obj))
      {
        UseExternalData(m_obj.data(), m_obj.rows(), m_obj.cols());
      }

      Wrap(const Wrap& other) = delete;
      Wrap(Wrap&&) = delete;
      void operator=(const Wrap&) = delete;
      void operator=(Wrap&&) = delete;

    private:
      Math::Matrix m_obj;
  };
}

#endif

