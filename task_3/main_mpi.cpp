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

class PointCoordinate {
public:
  int x;
  int y;
  int z;

  PointCoordinate(int _x, int _y, int _z) {
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

long double laplass(PointCoordinate p, Grid &grid, long double h) {
  long double firstFrac = (grid[p.x - 1][p.y][p.z] - 2 * grid[p.x][p.y][p.z] +
                           grid[p.x + 1][p.y][p.z]) /
                          (h * h);
  long double secondFrac = (grid[p.x][p.y - 1][p.z] - 2 * grid[p.x][p.y][p.z] +
                            grid[p.x][p.y + 1][p.z]) /
                           (h * h);
  long double thirdFrac = (grid[p.x][p.y][p.z - 1] - 2 * grid[p.x][p.y][p.z] +
                           grid[p.x][p.y][p.z + 1]) /
                          (h * h);

  return firstFrac + secondFrac + thirdFrac;
}

vector<Grid> getAnaliticalValues(const AxisPoints &axisPointsX,
                                 const AxisPoints &axisPointsY,
                                 const AxisPoints &axisPointsZ,
                                 const AxisPoints &axisPointsT, Limits l,
                                 const long &n) {
  vector<Grid> analiticalValues;
  for (size_t t = 0; t < T + 1; ++t) {
    analiticalValues.push_back(Grid());
    for (size_t i = 0; i < n + 1; ++i) {
      analiticalValues[t].push_back(YZProjection());
      for (size_t j = 0; j < n + 1; ++j) {
        analiticalValues[t][i].push_back(Values());
        for (size_t k = 0; k < n + 1; ++k) {
          analiticalValues[t][i][j].push_back(
              uAnalitical(Point(axisPointsX.points[i], axisPointsY.points[j],
                                axisPointsZ.points[k]),
                          l, axisPointsT.points[t]));
        }
      }
    }
  }
  return analiticalValues;
}

void calcInnerValue(vector<Grid> &g, PointCoordinate p, int t, long double &h,
                    long double &theta) {
  // cout << "Laplass: " << laplass(p, g[t], h) << endl;
  g[t + 1][p.x][p.y][p.z] = (theta * theta) * laplass(p, g[t], h) +
                            2 * g[t][p.x][p.y][p.z] - g[t - 1][p.x][p.y][p.z];
}

void calcBoundValue(vector<Grid> &g, PointCoordinate p, int t, int &n) {
  if (p.x == 0) {
    g[t + 1][p.x][p.y][p.z] =
        (g[t + 1][1][p.y][p.z] + g[t + 1][n - 2][p.y][p.z]) / 2.0;
  }
  if (p.y == 0) {
    g[t + 1][p.x][p.y][p.z] =
        (g[t + 1][p.x][1][p.z] + g[t + 1][p.x][n - 2][p.z]) / 2.0;
  }
  if (p.z == 0) {
    g[t + 1][p.x][p.y][p.z] =
        (g[t + 1][p.x][p.y][1] + g[t + 1][p.x][p.y][n - 2]) / 2.0;
  }
}

long double calcDiff(Grid &calcValues, Grid &analiticalValues, int &n) {
  long double maxDiff = abs(calcValues[0][0][0] - analiticalValues[0][0][0]);
  int maxI = 0, maxJ = 0, maxK = 0;
  for (int i = 1; i < n; ++i) {
    for (int j = 1; j < n; ++j) {
      for (int k = 1; k < n; ++k) {
        if (abs(calcValues[i][j][k] - analiticalValues[i][j][k]) > maxDiff) {
          maxI = i;
          maxJ = j;
          maxK = k;
          maxDiff = abs(calcValues[i][j][k] - analiticalValues[i][j][k]);
        }
      }
    }
  }
  // cout << maxI << " " << maxJ << " " << maxK << endl;
  // cout << "Analitical: " << analiticalValues[maxI][maxJ][maxK] << endl;
  // cout << "Calc: " << calcValues[maxI][maxJ][maxK] << endl;
  // cout << "Diff: " << maxDiff << endl;
  return maxDiff;
}

Grid calcU0(const AxisPoints &axisPointsX, const AxisPoints &axisPointsY,
            const AxisPoints &axisPointsZ, Limits &l, int &n) {
  Grid u0;

  for (int i = 0; i < n + 1; ++i) {
    u0.push_back(YZProjection());
    for (int j = 0; j < n + 1; ++j) {
      u0[i].push_back(Values());
      for (int k = 0; k < n + 1; ++k) {
        u0[i][j].push_back(
            phi(Point(axisPointsX.points[i], axisPointsY.points[j],
                      axisPointsZ.points[k]),
                l));
      }
    }
  }
  return u0;
}

Grid calcU1(Grid &u0, Grid &u1Analitical, Limits &l, int &n, long double &theta,
            long double &h) {
  Grid u1;

  for (int i = 0; i < n + 1; ++i) {
    u1.push_back(YZProjection());
    for (int j = 0; j < n + 1; ++j) {
      u1[i].push_back(Values());
      for (int k = 0; k < n + 1; ++k) {
        if (i == 0 || i == n || j == 0 || j == n || k == 0 || k == n) {
          u1[i][j].push_back(u1Analitical[i][j][k]);
        } else {
          u1[i][j].push_back(u0[i][j][k] +
                             (theta * theta) / 2.0 *
                                 laplass(PointCoordinate(i, j, k), u0, h));
        }
      }
    }
  }
  return u1;
}

int[3] computeBlockSize(long &n, int &coords, int &dims) {
  int blockSize[3];
  for (int i = 0; i < 3; ++i) {
    blockSize[i] = n / dims[i];
    if (coords[i] == dims[i] - 1) {
      blockSize[i] += N % dims[i];
    }
  }
  return blockSize;
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

  vector<Grid> analiticalValues = getAnaliticalValues(
      axisPointsX, axisPointsY, axisPointsZ, axisPointsT, l, n);

  MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  MPI_Comm main;
  int coords[3], dims[3];

  double startTime = MPI_Wtime();

  MPI_Cart_create(MPI_COMM_WORLD, 3, dims, {0, 1, 1}, 0, &main);
  MPI_Cart_coords(main, rank, 3, coords);

  int blockSize[3] = computeBlockSize(n, coords, dims);

  long double calcValues[T + 1][blockSize[0]][blockSize[1]][blockSize[2]];

  for (int i = 0; i < blockSize[0]; ++i) {
    for (int j = 0; j < blockSize[1]; ++j) {
      for (int k = 0; k < blockSize[2]; ++k) {
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

  // cout << "theta: " << axisPointsT.theta << endl;
  // cout << "h: " << axisPointsX.theta << endl;
  for (long t = 2; t < T + 1; t++) {
    calcValues.push_back(Grid());
    for (int i = 0; i < n + 1; ++i) {
      calcValues[t].push_back(YZProjection());
      for (int j = 0; j < n + 1; ++j) {
        calcValues[t][i].push_back(Values());
        for (int k = 0; k < n + 1; ++k) {
          calcValues[t][i][j].push_back(0.0);
        }
      }
    }

    for (int i = 1; i < n; ++i) {
      for (int j = 1; j < n; ++j) {
        for (int k = 1; k < n; ++k) {
          calcInnerValue(calcValues, PointCoordinate(i, j, k), t - 1,
                         axisPointsX.theta, axisPointsT.theta);
        }
      }
    }

    for (int i = 0; i < n + 1; ++i) {
      for (int j = 0; j < n + 1; ++j) {
        for (int k = 0; k < n + 1; ++k) {
          calcBoundValue(calcValues, PointCoordinate(i, j, k), t - 1, n);
        }
      }
    }
  }

  long double maxDiff = 0;
  for (int t = 0; t < T; t++) {
    // cout << "u" << t << endl;
    // cout << "Indexes: ";
    maxDiff = max(maxDiff, calcDiff(calcValues[t], analiticalValues[t], n));
    // cout << endl;
  }

  double endTime = MPI_Wtime() - startTime, totalTime;
  MPI_Reduce(&endTime, &total_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  long double totalDiff;
  MPI_Reduce(&maxDiff, &totalDiff, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if (!rank) {
    cout << "Time: " << totalTime << endl;
    cout << "Diff: " << totalDiff << endl;
  }

  MPI_Finalize();
  return 0;
}
