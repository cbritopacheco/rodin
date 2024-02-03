/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/Variational/LinearElasticity.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

    template <class T>
    struct IsTrialFunctionReferenceWrapper
    {
      static constexpr Boolean Value = false;
    };

    template <class T>
    struct IsTrialFunctionReferenceWrapper<std::reference_wrapper<T>>
    {
      static constexpr Boolean Value = IsTrialFunction<T>::Value;
    };

    template <class T>
    struct IsTestFunctionReferenceWrapper
    {
      static constexpr Boolean Value = false;
    };

    template <class T>
    struct IsTestFunctionReferenceWrapper<std::reference_wrapper<T>>
    {
      static constexpr Boolean Value = IsTestFunction<T>::Value;
    };

template <class T>
struct IsFloatingPoint
{
  static constexpr Boolean Value = std::is_floating_point_v<std::decay_t<T>>;
};

template <class T>
auto foo(const T& v)
{
  return 2 * v;
}

template <class T>
class Foo
{
  public:
    Foo(const T& v)
      : m_v(v)
    {}

    const T& v() const
    {
      return m_v;
    }

  private:
    T m_v;
};


int main(int argc, char** argv)
{
  // Load mesh
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, 16, 16);
  mesh.getConnectivity().compute(1, 2);

  // Functions
  P1 vh(mesh);

  TrialFunction u1(vh);
  TrialFunction u2(vh);
  TestFunction  v1(vh);
  TestFunction  v2(vh);
  TestFunction  v3(vh);

  // Tuple us{std::ref(u), std::ref(v), std::ref(u), std::ref(u)};
  // Tuple us{2, 'a', 4.0, true};
  // Tuple vs{1, 'c', 3.0, false};

  // using Miaow = Utility::Product<Tuple<int, char, double, bool>, Tuple<int, char, double, bool>>::Type<Pair>;
  // auto res = us.filter<IsTrialFunctionReferenceWrapper>();
  // auto z = us.zip<std::pair>(vs);
  // auto w = us.map([](auto& v) { return Foo(v); });

  // auto p = us.product<Pair>(vs);

  // w.apply([](auto& p){std::cout << p.v() << std::endl;});

  // Miaow t;
  // auto tt = t;
  // std::cout <<  res.size() << std::endl;
  Problem problem(u1, u2, v1, v2, v3);


  // trw.apply([](auto&& v){ std::cout << v << " "; });
  std::cout << std::endl;


  return 0;
}


