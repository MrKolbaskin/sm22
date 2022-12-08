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
  return phi(p, l) *
         cos(M_PI *
             sqrt((4.0 / l.x * l.x) + (16.0 / l.y * l.y) + (36.0 / l.z * l.z)) *
             t);
}

long double laplass(PointCoordinate p, Grid &grid, long n, long double h) {
  long double firstFrac = (grid[p.x - 1][p.y][p.z] - 2 * grid[p.x][p.y][p.z] +
                           grid[p.x + 1][p.y][p.z]) /
                          (h * h);
  long double secondFrac = (grid[p.x][p.y - 1][p.z] - 2 * grid[p.x][p.y][p.z] +
                            grid[p.x][p.y + 1][p.z]) /
                           (h * h);
  long double thirdFrac = (grid[p.x][p.y][p.z - 1] - 2 * grid[p.x][p.y][p.z] +
                           grid[p.x][p.y][p.z + 1]) /
                          (h * h);

  return (firstFrac + secondFrac + thirdFrac);
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

void calcInnerValue(vector<Grid> &g, PointCoordinate p, int t, int &n,
                    long double &h, long double &theta) {
  g[t + 1][p.x][p.y][p.z] = (theta * theta) * laplass(p, g[t], n, h) +
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
        // maxDiff =
        //     max(abs(calcValues[i][j][k] - analiticalValues[i][j][k]),
        //     maxDiff);
      }
    }
  }
  cout << maxI << " " << maxJ << " " << maxK << endl;
  cout << "Analitical: " << analiticalValues[maxI][maxJ][maxK] << endl;
  cout << "Calc: " << calcValues[maxI][maxJ][maxK] << endl;
  cout << "Diff: " << maxDiff << endl;
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
                                 laplass(PointCoordinate(i, j, k), u0, n, h));
        }
      }
    }
  }
  return u1;
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

  vector<Grid> analiticalValues = getAnaliticalValues(
      axisPointsX, axisPointsY, axisPointsZ, axisPointsT, l, n);

  vector<Grid> calcValues;

  cout << "theta: " << axisPointsT.theta << endl;
  cout << "h: " << axisPointsX.theta << endl;
  for (long t = 0; t < T; t++) {
    if (t == 0) {
      calcValues.push_back(calcU0(axisPointsX, axisPointsY, axisPointsZ, l, n));
    }
    if (t == 1) {
      calcValues.push_back(calcU1(calcValues[0], analiticalValues[1], l, n,
                                  axisPointsT.theta, axisPointsX.theta));
    }

    if (t > 1) {
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
            calcInnerValue(calcValues, PointCoordinate(i, j, k), t - 1, n,
                           axisPointsX.theta, axisPointsT.theta);
          }
        }
      }

      // for (int i = 0; i < n + 1; ++i) {
      //   for (int j = 0; j < n + 1; ++j) {
      //     for (int k = 0; k < n + 1; ++k) {
      //       calcBoundValue(calcValues, PointCoordinate(i, j, k), t - 1, n);
      //     }
      //   }
      // }
    }
  }

  // cout << "u1 " << analiticalValues[10][9][5][5] << endl;
  // cout << "calc u1 " << calcValues[10][9][5][5] << endl;

  for (int t = 0; t < T; t++) {
    cout << "u" << t << endl;
    cout << "Indexes: ";
    calcDiff(calcValues[t], analiticalValues[t], n);
    // cout << "u" << t << ": " << calcDiff(calcValues[t], analiticalValues[t],
    // n)
    //      << endl;
    cout << endl;
  }

  return 0;
}
