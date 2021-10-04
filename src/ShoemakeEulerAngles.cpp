#include <ShoemakeEulerAngles.hpp>
#include <float.h>
#include <math.h>

#ifndef PI_
#include <Declarations.hpp>
#endif

/*
 * Convert matrix to Euler angles (in radians).
 * Euler angle conventions following Ken Shoemake, Graphics Gems IV, Academic
 * Press Professional, Inc., 1999, 222–229.
 *
 * https://github.com/erich666/GraphicsGems/blob/master/gemsiv/euler_angle/
 */
void
Eul_FromHMatrix(Eigen::Matrix3d &M,
                double &alpha,
                double &beta,
                double &gamma,
                const unsigned int order,
                bool positiv)
{
  int i, j, k, h __attribute__((unused)), n, s, f;
  EulGetOrd(order, i, j, k, h, n, s, f);

  if (s == EulRepYes)
  {
    double sy = sqrt(M(i, j) * M(i, j) + M(i, k) * M(i, k));
    if (sy > 16 * DBL_EPSILON)
    {
      alpha = atan2(M(i, j), M(i, k));
      beta  = atan2(sy, M(i, i));
      gamma = atan2(M(j, i), -M(k, i));
    }
    else
    {
      alpha = atan2(-M(j, k), M(j, j));
      beta  = atan2(sy, M(i, i));
      gamma = .0;
    }
  }
  else
  {
    double cy = sqrt(M(i, i) * M(i, i) + M(j, i) * M(j, i));
    if (cy > 16 * DBL_EPSILON)
    {
      alpha = atan2(M(k, j), M(k, k));
      beta  = atan2(-M(k, i), cy);
      gamma = atan2(M(j, i), M(i, i));
    }
    else
    {
      alpha = atan2(-M(j, k), M(j, j));
      beta  = atan2(-M(k, i), cy);
      gamma = .0;
    }
  }
  if (n == EulParOdd)
  {
    alpha = -alpha;
    beta  = -beta;
    gamma = -gamma;
  }
  if (f == EulFrmR)
  {
    double t = alpha;
    alpha    = gamma;
    gamma    = t;
  }
  if (positiv && beta < .0)
  {
    beta = -beta;
    alpha += PI_;
    gamma += PI_;
  }
}
