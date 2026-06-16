#include <cmath>
#include "utils.h"

double dist2(double x1, double y1,
             double x2, double y2)
{
  double dx = x1 - x2;
  double dy = y1 - y2;

  return std::sqrt(dx*dx + dy*dy);
}

double dist3(double x1, double y1, double z1,
             double x2, double y2, double z2)
{
  double dx = x1 - x2;
  double dy = y1 - y2;
  double dz = z1 - z2;

  return std::sqrt(dx*dx + dy*dy + dz*dz);
}
