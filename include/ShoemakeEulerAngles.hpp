
#ifndef ALTERNATIVEEULERANGLES_HPP_
#define ALTERNATIVEEULERANGLES_HPP_

#include <eigen3/Eigen/Core>

#ifndef AXESDEFINITIONS_TITANIA_
#include <AxesDefinitions.hpp>
#endif

/*
 * Euler angle conventions following Ken Shoemake, Graphics Gems IV, Academic
 * Press Professional, Inc., 1999, 222–229.
 *
 * https://github.com/erich666/GraphicsGems/blob/master/gemsiv/euler_angle/
 */

/*
 * Order type constants, constructors, extractors
 *
 * There are 24 possible conventions, designated by:
 *  o EulAxI = axis used initially
 *  o EulPar = parity of axis permutation
 *  o EulRep = repetition of initial axis as last
 *  o EulFrm = frame from which axes are taken
 *
 * Axes I,J,K will be a permutation of X,Y,Z.
 *
 * Axis H will be either I or K, depending on EulRep.
 *
 * Frame S takes axes from initial static frame.
 *
 * If ord = (AxI=X, Par=Even, Rep=No, Frm=S), then {a,b,c,ord} means
 * Rz(c)Ry(b)Rx(a), where Rz(c)v rotates v around Z by c radians.
 */

#define EulFrmS     0
#define EulFrmR     1
#define EulFrm(ord) ((unsigned) (ord) &1)
#define EulRepNo    0
#define EulRepYes   1
#define EulRep(ord) (((unsigned) (ord) >> 1) & 1)
#define EulParEven  0
#define EulParOdd   1
#define EulPar(ord) (((unsigned) (ord) >> 2) & 1)

/* this code is merely a quick (and legal!) way to set arrays, EulSafe being
 * 0,1,2,0 */
#define EulSafe "\000\001\002\000"
#define EulNext "\001\002\000\001"

/* EulGetOrd unpacks all useful information about order simultaneously. */
#define EulGetOrd(ord, i, j, k, h, n, s, f) \
  { \
    unsigned o = (unsigned) ord; \
    f          = o & 1; \
    o >>= 1; \
    s = o & 1; \
    o >>= 1; \
    n = o & 1; \
    o >>= 1; \
    i = EulSafe[o & 3]; \
    j = EulNext[i + n]; \
    k = EulNext[i + 1 - n]; \
    h = s ? k : i; \
  }
/* EulOrd creates an order value between 0 and 23 from 4-tuple choices. */
#define EulOrd(i, p, r, f) (((((((i) << 1) + (p)) << 1) + (r)) << 1) + (f))
/* Static axes */
constexpr unsigned int EulOrdXYZs =
    EulOrd(STANDARD_AXIS_X_, EulParEven, EulRepNo, EulFrmS);
constexpr unsigned int EulOrdXYXs =
    EulOrd(STANDARD_AXIS_X_, EulParEven, EulRepYes, EulFrmS);
constexpr unsigned int EulOrdXZYs =
    EulOrd(STANDARD_AXIS_X_, EulParOdd, EulRepNo, EulFrmS);
constexpr unsigned int EulOrdXZXs =
    EulOrd(STANDARD_AXIS_X_, EulParOdd, EulRepYes, EulFrmS);
constexpr unsigned int EulOrdYZXs =
    EulOrd(STANDARD_AXIS_Y_, EulParEven, EulRepNo, EulFrmS);
constexpr unsigned int EulOrdYZYs =
    EulOrd(STANDARD_AXIS_Y_, EulParEven, EulRepYes, EulFrmS);
constexpr unsigned int EulOrdYXZs =
    EulOrd(STANDARD_AXIS_Y_, EulParOdd, EulRepNo, EulFrmS);
constexpr unsigned int EulOrdYXYs =
    EulOrd(STANDARD_AXIS_Y_, EulParOdd, EulRepYes, EulFrmS);
constexpr unsigned int EulOrdZXYs =
    EulOrd(STANDARD_AXIS_Z_, EulParEven, EulRepNo, EulFrmS);
constexpr unsigned int EulOrdZXZs =
    EulOrd(STANDARD_AXIS_Z_, EulParEven, EulRepYes, EulFrmS);
constexpr unsigned int EulOrdZYXs =
    EulOrd(STANDARD_AXIS_Z_, EulParOdd, EulRepNo, EulFrmS);
constexpr unsigned int EulOrdZYZs =
    EulOrd(STANDARD_AXIS_Z_, EulParOdd, EulRepYes, EulFrmS);
/* Rotating axes */
constexpr unsigned int EulOrdZYXr =
    EulOrd(STANDARD_AXIS_X_, EulParEven, EulRepNo, EulFrmR);
constexpr unsigned int EulOrdXYXr =
    EulOrd(STANDARD_AXIS_X_, EulParEven, EulRepYes, EulFrmR);
constexpr unsigned int EulOrdYZXr =
    EulOrd(STANDARD_AXIS_X_, EulParOdd, EulRepNo, EulFrmR);
constexpr unsigned int EulOrdXZXr =
    EulOrd(STANDARD_AXIS_X_, EulParOdd, EulRepYes, EulFrmR);
constexpr unsigned int EulOrdXZYr =
    EulOrd(STANDARD_AXIS_Y_, EulParEven, EulRepNo, EulFrmR);
constexpr unsigned int EulOrdYZYr =
    EulOrd(STANDARD_AXIS_Y_, EulParEven, EulRepYes, EulFrmR);
constexpr unsigned int EulOrdZXYr =
    EulOrd(STANDARD_AXIS_Y_, EulParOdd, EulRepNo, EulFrmR);
constexpr unsigned int EulOrdYXYr =
    EulOrd(STANDARD_AXIS_Y_, EulParOdd, EulRepYes, EulFrmR);
constexpr unsigned int EulOrdYXZr =
    EulOrd(STANDARD_AXIS_Z_, EulParEven, EulRepNo, EulFrmR);
constexpr unsigned int EulOrdZXZr =
    EulOrd(STANDARD_AXIS_Z_, EulParEven, EulRepYes, EulFrmR);
constexpr unsigned int EulOrdXYZr =
    EulOrd(STANDARD_AXIS_Z_, EulParOdd, EulRepNo, EulFrmR);
constexpr unsigned int EulOrdZYZr =
    EulOrd(STANDARD_AXIS_Z_, EulParOdd, EulRepYes, EulFrmR);

void Eul_FromHMatrix(Eigen::Matrix3d &M,
                     double &alpha,
                     double &beta,
                     double &gamma,
                     const unsigned int order,
                     bool positiv);

#endif
