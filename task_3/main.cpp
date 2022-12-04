#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>

using namespace std;

typedef vector<double> Values;
typedef vector<Values> YZProjection;
typedef vector<YZProjection> Grid;

const int T = 20;

class Limits {
public:
  double x;
  double y;
  double z;

  Limits(char **argv) {
    x = atof(argv[1]);
    y = atof(argv[2]);
    z = atof(argv[3]);
  }
};

class Point {
public:
  double x;
  double y;
  double z;

  Point(double _x, double _y, double _z) {
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
  double limit;
  long n;
  void fillAxisPoints() {
    double theta = limit / n, tmpPoint = 0;
    int i = 0;

    while (tmpPoint < limit) {
      points.push_back(tmpPoint);
      ++i;
      tmpPoint = i * theta;
    }
    points.push_back(limit);
  }

public:
  vector<double> points;

  AxisPoints(double _limit, long _n) {
    limit = _limit;
    n = _n;
    fillAxisPoints();
  }
};

double phi(Point p, Limits l) {
  return sin(2 * M_PI * p.x / l.x) * sin(4 * M_PI * p.y / l.y) *
         sin(6 * M_PI * p.z / l.z);
}

double uAnalitical(Point p, Limits l, double t) {
  return phi(p, l);
  // *
  //      cos(M_PI *
  //          sqrt((4.0 / l.x * l.x) + (16.0 / l.y * l.y) + (36.0 / l.z * l.z))
  //          * t);
}

double laplass(PointCoordinate p, Grid &grid, long n) {
  double denom = 0.01; // n * n;
  double firstFrac = (grid[(n + p.x - 1) % n][p.y][p.z] -
                      2 * grid[p.x][p.y][p.z] + grid[(p.x + 1) % n][p.y][p.z]);
  double secondFrac = (grid[p.x][(n + p.y - 1) % n][p.z] -
                       2 * grid[p.x][p.y][p.z] + grid[p.x][(p.y + 1) % n][p.z]);
  double thirdFrac = (grid[p.x][p.y][(n + p.z - 1) % n] -
                      2 * grid[p.x][p.y][p.z] + grid[p.x][p.y][(p.z + 1) % n]);
  return (firstFrac + secondFrac + thirdFrac) / denom;
}

vector<Grid> getAnaliticalValues(const AxisPoints &axisPointsX,
                                 const AxisPoints &axisPointsY,
                                 const AxisPoints &axisPointsZ, Limits l,
                                 const long &n) {
  vector<Grid> analiticalValues;
  for (size_t t = 0; t < T; ++t) {
    analiticalValues.push_back(Grid());
    for (size_t i = 0; i < n; ++i) {
      analiticalValues[t].push_back(YZProjection());
      for (size_t j = 0; j < n; ++j) {
        analiticalValues[t][i].push_back(Values());
        for (size_t k = 0; k < n; ++k) {
          analiticalValues[t][i][j].push_back(
              uAnalitical(Point(axisPointsX.points[i], axisPointsY.points[j],
                                axisPointsZ.points[k]),
                          l, t));
        }
      }
    }
  }
  return analiticalValues;
}

void calcNextValue(vector<Grid> &g, PointCoordinate p, int t, long &n) {
  g[t + 1][p.x][p.y][p.z] =
      laplass(p, g[t], n) + 2 * g[t][p.x][p.y][p.z] - g[t - 1][p.x][p.y][p.z];
}

double calcDiff(Grid &calcValues, Grid &analiticalValues, long &n) {
  double maxDiff = abs(calcValues[0][0][0] - analiticalValues[0][0][0]);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int k = 0; k < n; ++k) {
        maxDiff =
            max(abs(calcValues[i][j][k] - analiticalValues[i][j][k]), maxDiff);
      }
    }
  }
  return maxDiff;
}

Grid calcU0(const AxisPoints &axisPointsX, const AxisPoints &axisPointsY,
            const AxisPoints &axisPointsZ, Limits &l, long &n) {
  Grid u0;

  for (int i = 0; i < n; ++i) {
    u0.push_back(YZProjection());
    for (int j = 0; j < n; ++j) {
      u0[i].push_back(Values());
      for (int k = 0; k < n; ++k) {
        u0[i][j].push_back(
            phi(Point(axisPointsX.points[i], axisPointsY.points[j],
                      axisPointsZ.points[k]),
                l));
      }
    }
  }
  return u0;
}

Grid calcU1(Grid &u0, Limits &l, long &n) {
  Grid u1;

  for (int i = 0; i < n; ++i) {
    u1.push_back(YZProjection());
    for (int j = 0; j < n; ++j) {
      u1[i].push_back(Values());
      for (int k = 0; k < n; ++k) {
        u1[i][j].push_back(u0[i][j][k] +
                           0.5 * laplass(PointCoordinate(i, j, k), u0, n));
      }
    }
  }
  return u1;
}

int main(int argc, char *argv[]) {
  if (argc < 5) {
    cout << "usage: ./run lx ly lz N" << endl;
    return 1;
  }
  Limits l(argv);
  long n = atol(argv[4]);
  AxisPoints axisPointsX(l.x, n), axisPointsY(l.y, n), axisPointsZ(l.z, n);

  vector<Grid> analiticalValues =
      getAnaliticalValues(axisPointsX, axisPointsY, axisPointsZ, l, n);

  vector<Grid> calcValues;

  for (int t = 0; t < T; t++) {
    if (t == 0) {
      calcValues.push_back(calcU0(axisPointsX, axisPointsY, axisPointsZ, l, n));
    }
    if (t == 1) {
      calcValues.push_back(calcU1(calcValues[0], l, n));
    }

    if (t > 1) {
      calcValues.push_back(Grid());
      for (int i = 0; i < n; ++i) {
        calcValues[t].push_back(YZProjection());
        for (int j = 0; j < n; ++j) {
          calcValues[t][i].push_back(Values());
          for (int k = 0; k < n; ++k) {
            calcValues[t][i][j].push_back(0.0);
            calcNextValue(calcValues, PointCoordinate(i, j, k), t - 1, n);
          }
        }
      }
    }
  }

  for (int t = 0; t < T; t++) {
    cout << calcDiff(calcValues[t], analiticalValues[0], n) << endl;
  }

  return 0;
}
