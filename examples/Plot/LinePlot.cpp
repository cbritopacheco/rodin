/**
 * @file LinePlot.cpp
 * Demonstrates a very simple use of the Rodin::Plot API by plotting a sine
 * curve.
 */
#include <Rodin/Plot.h>

using namespace std;
using namespace Rodin::Plot;

using Array = Eigen::ArrayX<float>;

const double pi = 4 * std::atan(1);

int main(int, char**)
{
  Plot plt;

  Array x = Array::LinSpaced(100, 0, 2 * pi);
  Array y = Eigen::sin(x);

  auto& fig = plt.figure(800, 600);
  auto& ax = fig.addAxes();
  ax.plot(x, y);
  plt.show();

  return 0;
}
