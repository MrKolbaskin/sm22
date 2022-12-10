#include "mpi.h"
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>

using namespace std;

typedef vector<long double> Values;
typedef vector<Values> YZProjection;
typedef vector<YZProjection> Grid;

const int T = 20;

class Limits {
public:
  long double x;
  long double y;
  long double z;

  Limits(char **argv) {
    x = atof(argv[1]);
    y = atof(argv[2]);
    z = atof(argv[3]);
  }
};

class Point {
public:
  long double x;
  long double y;
  long double z;

  Point(long double _x, long double _y, long double _z) {
    x = _x;
    y = _y;
    z = _z;
  }
};

class AxisPoints {
  long double limit;
  long n;
  void fillAxisPoints() {
    theta = limit / n;
    long double tmpPoint = 0;
    int i = 0;

    while (tmpPoint < limit) {
      points.push_back(tmpPoint);
      ++i;
      tmpPoint = i * theta;
    }
    points.push_back(limit);
  }

public:
  vector<long double> points;
  long double theta;

  AxisPoints(long double _limit, long _n) {
    limit = _limit;
    n = _n;
    fillAxisPoints();
  }
};

long double phi(Point p, Limits l) {
  return sin(2 * M_PI * p.x / l.x) * sin(4 * M_PI * p.y / l.y) *
         sin(6 * M_PI * p.z / l.z);
}

long double uAnalitical(Point p, Limits l, long double t) {
  return phi(p, l) * cos(M_PI *
                         sqrt((4.0 / (l.x * l.x)) + (16.0 / (l.y * l.y)) +
                              (36.0 / (l.z * l.z))) *
                         t);
}

int main(int argc, char *argv[]) {
  if (argc < 6) {
    cout << "usage: ./run lx ly lz T N" << endl;
    return 1;
  }
  Limits l(argv);
  long double maxT = atof(argv[4]);
  int n = atoi(argv[5]);
  AxisPoints axisPointsX(l.x, n), axisPointsY(l.y, n), axisPointsZ(l.z, n),
      axisPointsT(maxT, T);

  long double h = axisPointsX.theta;
  long double tau = axisPointsT.theta;

  MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  MPI_Comm main;
  int coords[3], dims[3] = {0, 0, 0}, periods[3] = {0, 1, 1};

  double startTime = MPI_Wtime();

  MPI_Dims_create(size, 3, dims);
  MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &main);
  MPI_Cart_coords(main, rank, 3, coords);

  int blockSize[3];
  for (int i = 0; i < 3; ++i) {
    blockSize[i] = n / dims[i];
    if (coords[i] == dims[i] - 1) {
      blockSize[i] += n % dims[i];
    }
  }
  double calcValues[T + 1][blockSize[0]][blockSize[1]][blockSize[2]],
      analiticalValues[T + 1][blockSize[0]][blockSize[1]][blockSize[2]];

  size_t t, i, j, k;
  // analitical values
#pragma omp parallel for shared(analiticalValues) collapse(4);
  for (t = 0; t < T + 1; ++t) {
    for (i = 0; i < blockSize[0]; ++i) {
      for (j = 0; j < blockSize[1]; ++j) {
        for (k = 0; k < blockSize[2]; ++k) {
          long double x = (i + coords[0] * blockSize[0]) * h;
          long double y = (j + coords[1] * blockSize[1]) * h;
          long double z = (k + coords[2] * blockSize[2]) * h;
          analiticalValues[t][i][j][k] =
              uAnalitical(Point(x, y, z), l, t * tau);
        }
      }
    }
  }

#pragma omp parallel for shared(analiticalValues) collapse(3);
  for (i = 0; i < blockSize[0]; ++i) {
    for (j = 0; j < blockSize[1]; ++j) {
      for (k = 0; k < blockSize[2]; ++k) {
        long double x = (i + coords[0] * blockSize[0]) * h;
        long double y = (j + coords[1] * blockSize[1]) * h;
        long double z = (k + coords[2] * blockSize[2]) * h;
        calcValues[0][i][j][k] = phi(Point(x, y, z), l);

        long double firstFrac =
            (phi(Point(x - h, y, z), l) - 2 * calcValues[0][i][j][k] +
             phi(Point(x + h, y, z), l)) /
            (h * h);
        long double secondFrac =
            (phi(Point(x, y - h, z), l) - 2 * calcValues[0][i][j][k] +
             phi(Point(x, y + h, z), l)) /
            (h * h);
        long double thirdFrac =
            (phi(Point(x, y, z - h), l) - 2 * calcValues[0][i][j][k] +
             phi(Point(x, y, z + h), l)) /
            (h * h);

        double long laplassValue = firstFrac + secondFrac + thirdFrac;

        calcValues[1][i][j][k] =
            calcValues[0][i][j][k] + (tau * tau) / 2 * laplassValue;
      }
    }
  }

  int maxBlockSize = max(blockSize[0], max(blockSize[1], blockSize[2]));
  double outPrev[3][maxBlockSize * maxBlockSize],
      outNext[3][maxBlockSize * maxBlockSize];
  double inPrev[3][maxBlockSize * maxBlockSize],
      inNext[3][maxBlockSize * maxBlockSize];
  int srcRank, destRank;

  for (t = 2; t < T + 1; ++t) {
    /* Edge exchanging */
    for (int d = 0; d < 3; ++d) {
      int blockSizeI, blockSizeJ;

      switch (d) {
      case 0:
        blockSizeI = blockSize[1];
        blockSizeJ = blockSize[2];
        break;
      case 1:
        blockSizeI = blockSize[0];
        blockSizeJ = blockSize[2];
        break;
      case 2:
        blockSizeI = blockSize[0];
        blockSizeJ = blockSize[1];
        break;
      }

      /* Previous neighbor */
      if (coords[d] != 0) {
        for (int i = 0; i < blockSizeI; i++) {
          for (int j = 0; j < blockSizeJ; j++) {
            switch (d) {
            case 0:
              outPrev[0][j * blockSizeI + i] = calcValues[t - 1][0][i][j];
              break;
            case 1:
              outPrev[1][j * blockSizeI + i] = calcValues[t - 1][i][0][j];
              break;
            case 2:
              outPrev[2][j * blockSizeI + i] = calcValues[t - 1][i][j][0];
              break;
            }
          }
        }
        MPI_Cart_shift(main, d, -1, &srcRank, &destRank);
        MPI_Sendrecv(outPrev[d], maxBlockSize * maxBlockSize, MPI_DOUBLE,
                     destRank, 1, inPrev[d], maxBlockSize * maxBlockSize,
                     MPI_DOUBLE, destRank, 1, main, MPI_STATUS_IGNORE);
      }

      /* Next neighbor */
      if (coords[d] != dims[d] - 1) {
        for (int i = 0; i < blockSizeI; i++) {
          for (int j = 0; j < blockSizeJ; j++) {
            switch (d) {
            case 0:
              outNext[0][j * blockSizeI + i] =
                  calcValues[t - 1][blockSize[0] - 1][i][j];
              break;
            case 1:
              outNext[1][j * blockSizeI + i] =
                  calcValues[t - 1][i][blockSize[1] - 1][j];
              break;
            case 2:
              outNext[2][j * blockSizeI + i] =
                  calcValues[t - 1][i][j][blockSize[2] - 1];
              break;
            }
          }
        }
        MPI_Cart_shift(main, d, 1, &srcRank, &destRank);
        MPI_Sendrecv(outNext[d], maxBlockSize * maxBlockSize, MPI_DOUBLE,
                     destRank, 1, inNext[d], maxBlockSize * maxBlockSize,
                     MPI_DOUBLE, destRank, 1, main, MPI_STATUS_IGNORE);
      }
    }

#pragma omp parallel for shared(inPrev, inNext, calcValues) collapse(3);
    for (i = 0; i < blockSize[0]; ++i) {
      for (j = 0; j < blockSize[1]; ++j) {
        for (k = 0; k < blockSize[2]; ++k) {
          long double firstFrac = 0, secondFrac = 0, thirdFrac = 0;

          if (i == 0) {
            firstFrac = (inPrev[0][k * blockSize[1] + j] -
                         2 * calcValues[t - 1][i][j][k] +
                         calcValues[t - 1][i + 1][j][k]) /
                        (h * h);

          } else if (i == blockSize[0] - 1) {
            firstFrac = (calcValues[t - 1][i - 1][j][k] -
                         2 * calcValues[t - 1][i][j][k] +
                         inNext[0][k * blockSize[1] + j]) /
                        (h * h);

          } else {
            firstFrac = (calcValues[t - 1][i - 1][j][k] -
                         2 * calcValues[t - 1][i][j][k] +
                         calcValues[t - 1][i + 1][j][k]) /
                        (h * h);
          }

          if (j == 0) {
            secondFrac = (inPrev[1][k * blockSize[0] + i] -
                          2 * calcValues[t - 1][i][j][k] +
                          calcValues[t - 1][i][j + 1][k]) /
                         (h * h);
          } else if (j == blockSize[1] - 1) {
            secondFrac = (calcValues[t - 1][i][j - 1][k] -
                          2 * calcValues[t - 1][i][j][k] +
                          inNext[1][k * blockSize[0] + i]) /
                         (h * h);
          } else {
            secondFrac = (calcValues[t - 1][i][j - 1][k] -
                          2 * calcValues[t - 1][i][j][k] +
                          calcValues[t - 1][i][j + 1][k]) /
                         (h * h);
          }

          if (k == 0) {
            thirdFrac = (inPrev[2][j * blockSize[0] + i] -
                         2 * calcValues[t - 1][i][j][k] +
                         calcValues[t - 1][i][j][k + 1]) /
                        (h * h);
          } else if (k == blockSize[2] - 1) {
            thirdFrac = (calcValues[t - 1][i][j][k - 1] -
                         2 * calcValues[t - 1][i][j][k] +
                         inNext[2][j * blockSize[0] + i]) /
                        (h * h);
          } else {
            thirdFrac = (calcValues[t - 1][i][j][k - 1] -
                         2 * calcValues[t - 1][i][j][k] +
                         calcValues[t - 1][i][j][k + 1]) /
                        (h * h);
          }
          long double laplassValue = firstFrac + secondFrac + thirdFrac;
          calcValues[t][i][j][k] = laplassValue * tau * tau +
                                   2 * calcValues[t - 1][i][j][k] -
                                   calcValues[t - 2][i][j][k];
        }
      }
    }
  }

  long double maxDiff = 0;
  for (t = 0; t < T + 1; t++) {
    long double tmpMaxDiff =
        abs(calcValues[t][0][0][0] - analiticalValues[t][0][0][0]);
    int maxI = 0, maxJ = 0, maxK = 0;
    for (i = 1; i < blockSize[0]; ++i) {
      for (j = 1; j < blockSize[1]; ++j) {
        for (k = 1; k < blockSize[2]; ++k) {
          if (abs(calcValues[t][i][j][k] - analiticalValues[t][i][j][k]) >
              tmpMaxDiff) {
            maxI = i;
            maxJ = j;
            maxK = k;
            tmpMaxDiff =
                abs(calcValues[t][i][j][k] - analiticalValues[t][i][j][k]);
          }
        }
      }
    }
    maxDiff = max(tmpMaxDiff, maxDiff);
  }

  double endTime = MPI_Wtime() - startTime, totalTime;
  MPI_Reduce(&endTime, &totalTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  long double totalDiff;
  MPI_Reduce(&maxDiff, &totalDiff, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if (!rank) {
    cout << "Time: " << totalTime << endl;
    cout << "Diff: " << totalDiff << endl;
  }

  MPI_Finalize();
  return 0;
}
