#include "geometry.hpp"
#include <boost/math/differentiation/autodiff.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>
#include <fstream>
#include <iostream>
#include <memory>
#include <unordered_map>

#ifdef USE_CACHE
#define CACHE_DECL mutable std::unordered_map<double, Point> cache;
#define CACHE_CHECK(t)                                                         \
  auto cached = cache.find(t);                                                 \
  if (cached != cache.end()) {                                                 \
    return cached->second;                                                     \
  }
#define CACHE_STORE(t, result) cache[t] = result;
#else
#define CACHE_DECL
#define CACHE_CHECK(t)
#define CACHE_STORE(t, result)
#endif


Point operator*(double factor, const Point &p) { return p * factor; }

class Curve {

public:
  virtual ~Curve() = default;
  virtual Point at(double t) const = 0;
};

class StraightLine : public Curve {
public:
  StraightLine(const Point &start, const Point &end) : start(start), end(end) {}

  StraightLine(const StraightLine &other)
      : start(other.start), end(other.end) {}

  StraightLine &operator=(const StraightLine &other) {
    start = other.start;
    end = other.end;
    return *this;
  }
  StraightLine &operator=(StraightLine &&other) noexcept {
    if (this != &other) {
      start = other.start;
      end = other.end;
    }
    return *this;
  }

  Point at(double t) const override { return start + (end - start) * t; }

  ~StraightLine() {}

private:
  Point start;
  Point end;
};

class EquationCurve : public Curve {
public:
  static constexpr int MAX_ITER = 100;
  virtual ~EquationCurve() = default;
  Point at(double t) const override {

    CACHE_CHECK(t);

    double s1 = arcLength(1.0);

    double ti = t;
    double sTi = arcLength(ti);

    int iter = 0;
    while (iter < MAX_ITER &&
           std::abs(sTi - t * s1) > std::numeric_limits<double>::epsilon()) {
      double step = (sTi - t * s1) / gammaprime(ti).norm();

      ti = std::clamp(ti - step, 0.0, 1.0);
      sTi = arcLength(ti);
      iter++;
    }
    auto result = gamma(ti);

    CACHE_STORE(t, result);

    return result;
  }

protected:
  CACHE_DECL

private:
  double arcLength(double t) const {
    auto result = boost::math::quadrature::trapezoidal(
        [&](double tau) {
          auto gp = gammaprime(tau);
          return gp.norm();
        },
        0., t);

    return result;
  }
  virtual Point gamma(double t) const = 0;
  virtual Point gammaprime(double t) const = 0;
};

class BottomCurve : public EquationCurve {
public:
  ~BottomCurve() {  }

private:
  double g(double x) const {
    if (x >= -10 && x < -3) {
      return 1 + std::exp(-3 * (x + 6.0));
    }
    if (x >= -3 && x <= 5) {
      return 1 + std::exp(3.0 * x);
    }

    std::string errorMsg = "x { " + std::to_string(x) + " } not in domain";
    throw std::domain_error(errorMsg);
  }

  double gPrime(double x) const {
    if (x >= -10 && x < -3) {
      return -3 * std::exp(-3 * (x + 6.0));
    } else if (x >= -3 && x <= 5) {
      return 3 * std::exp(3.0 * x);
    } else {
      std::string errorMsg =
          "x { " + std::to_string(x) + " } not in prime domain";
      throw std::domain_error(errorMsg);
    }
  }

  Point gamma(double t) const override {
    double x = (1 - t) * -10 + t * 5;
    return {x, 1. / (2 * g(x))};
  }

  Point gammaprime(double t) const override {
    double x = (1 - t) * -10 + t * 5;
    double y_derivative = -1.0 / (2.0 * std::pow(g(x), 2)) * gPrime(x) * 15;

    return {15, y_derivative};
  };
};

class Domain {

public:
  Domain(int n, std::unique_ptr<Curve> lb, std::unique_ptr<Curve> rb,
         std::unique_ptr<Curve> tb, std::unique_ptr<Curve> bb) {

    leftBound = std::move(lb);
    rightBound = std::move(rb);
    topBound = std::move(tb);
    bottomBound = std::move(bb);

    grid = std::vector<std::vector<Point>>(n, std::vector<Point>(n));
    double h = 1.0 / (n - 1);

    auto r00 = bottomBound->at(0);
    auto r01 = topBound->at(0);
    auto r10 = bottomBound->at(1);
    auto r11 = topBound->at(1);

    for (size_t i = 0; i < n; i++) {
      auto xi = static_cast<double>(i) * h;
      auto top = topBound->at(xi);
      auto bottom = bottomBound->at(xi);
      for (size_t j = 0; j < n; j++) {
        auto eta = static_cast<double>(j) * h;
        auto left = leftBound->at(eta);
        auto right = rightBound->at(eta);
        grid[i][j] = (1 - xi) * left + xi * right + (1 - eta) * bottom +
                     eta * top -
                     ((1 - xi) * (1 - eta) * r00 + (1 - xi) * eta * r01 +
                      xi * (1 - eta) * r10 + xi * eta * r11);
      }
    }
  }

  std::vector<std::vector<Point>> grid;

private:
  std::unique_ptr<Curve> leftBound;
  std::unique_ptr<Curve> rightBound;
  std::unique_ptr<Curve> topBound;
  std::unique_ptr<Curve> bottomBound;
};

int main(int argc, char **argv) {

  int n = 5;

  if (argc > 1) {
    n = std::stoi(argv[1]);
  }

  auto bottomBound = std::make_unique<BottomCurve>();

  auto testDomain =
      Domain(n, std::make_unique<StraightLine>(Point{-10, 0}, Point{-10, 3}),
             std::make_unique<StraightLine>(Point{5, 0}, Point{5, 3}),
             std::make_unique<StraightLine>(Point{-10, 3}, Point{5, 3}),
             std::move(bottomBound));

  std::ofstream x_file("x_coords.txt");
  std::ofstream y_file("y_coords.txt");

  for (auto &row : testDomain.grid) {
    auto first = 1;
    for (auto &point : row) {
      if (first) {
        x_file << point.x;
        first = 0;
      }
      x_file << " " << point.x;
    }
    x_file << "\n";
  }

  for (auto &row : testDomain.grid) {
    auto first = 1;
    for (auto &point : row) {
      if (first) {
        y_file << point.y;
        first = 0;
      }

      y_file << " " << point.y;
    }
    y_file << "\n";
  }
}
