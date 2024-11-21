#pragma once

#include <cmath>

class Point {
public:
  double x, y;

  Point() : x(0), y(0) {}
  Point(double x, double y) : x(x), y(y) {}

  Point(const Point &other) : x(other.x), y(other.y) {}

  Point operator+(const Point &other) const {
    return Point(x + other.x, y + other.y);
  }

  Point operator-(const Point &other) const {
    return Point(x - other.x, y - other.y);
  }

  Point operator*(double factor) const { return Point(x * factor, y * factor); }

  Point(Point &&other) noexcept : x(other.x), y(other.y) {
    other.x = 0;
    other.y = 0;
  }

  Point &operator=(const Point &other) {
    if (this != &other) {
      x = other.x;
      y = other.y;
    }
    return *this;
  }

  Point &operator=(Point &&other) noexcept {
    if (this != &other) {
      x = other.x;
      y = other.y;
      other.x = 0;
      other.y = 0;
    }
    return *this;
  }

  double norm() const { return std::sqrt(x * x + y * y); }
};

class Rectangle {
public:
  Point bottomLeft;
  Point topRight;

  Point topLeft() const { return Point(bottomLeft.x, topRight.y); }
  Point bottomRight() const { return Point(topRight.x, bottomLeft.y); }

  Point center() const {
    return Point((bottomLeft.x + topRight.x) / 2,
                 (bottomLeft.y + topRight.y) / 2);
  }

  bool contains(const Point &p) const {
    return bottomLeft.x <= p.x && p.x <= topRight.x && bottomLeft.y <= p.y &&
           p.y <= topRight.y;
  }

  bool overlaps(const Rectangle &r) const {
    return bottomLeft.x < r.topRight.x && topRight.x > r.bottomLeft.x &&
           bottomLeft.y < r.topRight.y && topRight.y > r.bottomLeft.y;
  }

  Rectangle() : bottomLeft(Point(0, 0)), topRight(Point(0, 0)) {}
  Rectangle(Point bl, Point tr) : bottomLeft(bl), topRight(tr) {}
  Rectangle(Rectangle &r) : bottomLeft(r.bottomLeft), topRight(r.topRight) {}
};
