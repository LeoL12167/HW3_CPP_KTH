// HW3_CPP_KTH.cpp: Definiert den Einstiegspunkt für die Anwendung.
//
#include <iostream>
#include <functional>
#include <cmath>
#include <stdexcept>
#include <fstream>
#include <vector>
#include <iomanip>


class Point {
public:
    double x;
    double y;
    Point() : x(0), y(0) {}
    Point(double x_coordinate, double y_coordinate) : x(x_coordinate), y(y_coordinate) {}
};

class Curve {
public:
    virtual ~Curve() = default;
    virtual Point at(double t) const = 0;
    virtual Point at_xi(double xi) const {
        return at(xi);
    }
};

class StraightLine : public Curve {
public:
    Point start, end;
    StraightLine(const Point& start, const Point& end) : start(start), end(end) {}
    Point at(double t) const override {
        double x = (1 - t) * start.x + t * end.x;
        double y = (1 - t) * start.y + t * end.y;
        return Point(x, y);
    }
};


class EquationCurve : public Curve
{
public:


    virtual ~EquationCurve() = default;
    Point at(double t) const override {
        // Step 1: Map t to arclength s (0 to total arclength)
        double s = t * full_arclength();

        // Step 2: Use arclength to find corresponding parameter t'
        double t_prime = find_t_for_arclenght(s);

        // Step 3: Evaluate gamma at the reparameterized t
        return gamma(t_prime);
    };

protected:
    virtual Point gamma(double t) const = 0;
    virtual double full_arclength() const = 0;
    virtual double find_t_for_arclenght(double s, double tol = 10e-6) const = 0;
    Point gammaprime(double t, double epsilon = 10e-5) const {
        
        try {
            double dx = (gamma(t + epsilon).x - gamma(t - epsilon).x) / (2 * epsilon);
            double dy = (gamma(t + epsilon).y - gamma(t - epsilon).y) / (2 * epsilon);
            return Point(dx, dy);
        }
        catch (const std::runtime_error&) {
            try {
                double dx = (gamma(t + epsilon).x - gamma(t).x) / (epsilon);
                double dy = (gamma(t + epsilon).y - gamma(t).y) / (epsilon);
                return Point(dx, dy);
            }
            catch (const std::runtime_error&) {
                double dx = (-gamma(t - epsilon).x + gamma(t).x) / (epsilon);
                double dy = (-gamma(t - epsilon).y + gamma(t).y) / (epsilon);
                return Point(dx, dy);
            }
        }
    }

};


class BottomCurve : public EquationCurve {
public:

    Point start, end;
    BottomCurve(const Point& start, const Point& end)
        : start(start), end(end) {
        full_arclength();
    }


    double full_arclength() const override {
        return ASI(0, 0, 1, 10e-4);
    }

    double arclength(double t) const {
        return ASI(0, 0, t, 10e-4);
    };

    double find_t_for_arclenght(double s, double tol = 10e-6) const override {
        double low = 0.0, high = 1.0, mid;
        while (high - low > tol) {
            mid = (low + high) / 2;
            double arc_length = arclength(mid);
            if (arc_length < s) {
                low = mid;
            }
            else {
                high = mid;
            }
        }
        return mid;
    }
        
        
        
    double iterand(double t)const {
        Point derivative = gammaprime(t);
        return std::sqrt(derivative.x * derivative.x + derivative.y * derivative.y);
    };

    double Simpson_rule(double min, double max) const{
        double Mid = (min + max) / 2;
        return (max - min) / 6 * (iterand(min) + 4 * iterand(Mid) + iterand(max));
    }

    double ASI(double I2, const double& intervalmin, const double& intervalmax, const double& tolerance)const {
        if (tolerance <= 0) {
            throw std::invalid_argument("Negative value is not allowed");
        }

        double midpoint = (intervalmin + intervalmax) / 2;

        // Simpson's rule for the current interval
        double I = Simpson_rule(intervalmin, intervalmax);

        // Simpson's rule for two subintervals
        I2 = Simpson_rule(intervalmin, midpoint) + Simpson_rule(midpoint, intervalmax);

        double errorest = (I2 - I) / 15;

        // Check if error is within tolerance
        if (std::fabs(errorest) < tolerance) {
            return I2;
        }
        else {
            // Recursively apply ASI on two subintervals
            return ASI(0, intervalmin, midpoint, tolerance / 2) + ASI(0, midpoint, intervalmax, tolerance / 2);
        }
    }
  
    
    Point gamma(double t) const override {
        double x = (1 - t) * (start.x) + t * end.x;
        double y;
        if (x >= start.x && x <= -3) {
            y = 1 / (2 *(1 + std::exp(-3 * (x + 6))));
        }
        else if (x > -3 && x <= end.x) {
            y = 1 / (2 * (1 + std::exp(3 * x)));
        }
        else if (x < start.x || x > end.x) {
            throw std::runtime_error("Value of x is out of range.");
        }
        return Point(x, y);
    }
};


class Domain {
public:
    const Curve& left;
    const Curve& right;
    const Curve& top;
    const Curve& bottom;

    Domain(const Curve& left, const Curve& right,
        const Curve& top, const Curve& bottom)
        : left(left), right(right), top(top), bottom(bottom) {}
};


class FTI {
public:
    Point coordinate(double xi, double eta,
        const Curve& bottom, const Curve& top,
        const Curve& left, const Curve& right,
        const Point& bottom_left, const Point& bottom_right,
        const Point& top_left, const Point& top_right) {
        // Compute the blended coordinates
        Point bottom_point = bottom.at_xi(xi);
        Point top_point = top.at_xi(xi);
        Point left_point = left.at_xi(eta);
        Point right_point = right.at_xi(eta);

        double x = (1 - eta) * bottom_point.x + eta * top_point.x +
            (1 - xi) * left_point.x + xi * right_point.x -
            ((1 - xi) * (1 - eta) * bottom_left.x +
                xi * (1 - eta) * bottom_right.x +
                (1 - xi) * eta * top_left.x +
                xi * eta * top_right.x);

        double y = (1 - eta) * bottom_point.y + eta * top_point.y +
            (1 - xi) * left_point.y + xi * right_point.y -
            ((1 - xi) * (1 - eta) * bottom_left.y +
                xi * (1 - eta) * bottom_right.y +
                (1 - xi) * eta * top_left.y +
                xi * eta * top_right.y);

        return Point(x, y);
    }
};



void writeMatrixToFile(const std::string& filename, const std::vector<std::vector<double>>& matrix) {
    std::ofstream file(filename);
    if (file.is_open()) {
        for (const auto& row : matrix) {
            for (size_t j = 0; j < row.size(); ++j) {
                file << std::fixed << std::setprecision(4) << row[j];
                if (j < row.size() - 1) file << " ";
            }
            file << "\n";
        }
        file.close();
        std::cout << "Matrix written to " << filename << "\n";
    }
    else {
        std::cerr << "Unable to open file: " << filename << "\n";
    }
}
int main() {
    // Define grid dimensions and step size
    int grid_size = 50;  // Example grid size
    double h = 1.0 / (grid_size - 1);

    // Define corner points for the domain
    Point bottom_left(-10.0, 0);
    Point bottom_right(5.0, 0);
    Point top_left(-10.0, 3.0);
    Point top_right(5.0, 3.0);

    // Define boundary curves
    //StraightLine bottom_curve(bottom_left, bottom_right);
    BottomCurve bottom_curve(bottom_left, bottom_right);
    StraightLine top_curve(top_left, top_right);
    StraightLine left_curve(bottom_left, top_left);
    StraightLine right_curve(bottom_right, top_right);

    // Create matrices to store x and y grid points
    std::vector<std::vector<double>> x_matrix(grid_size, std::vector<double>(grid_size));
    std::vector<std::vector<double>> y_matrix(grid_size, std::vector<double>(grid_size));

    // Create an instance of FTI
    FTI fti;

    // Populate x_matrix and y_matrix with interpolated values
    for (int i = 0; i < grid_size; ++i) {
        double xi = i * h;
        for (int j = 0; j < grid_size; ++j) {
            double eta = j * h;
            Point p = fti.coordinate(xi, eta, bottom_curve, top_curve, left_curve, right_curve,
                bottom_left, bottom_right, top_left, top_right);
            x_matrix[i][j] = p.x;
            y_matrix[i][j] = p.y;
        }
    }

    // Write matrices to files
    writeMatrixToFile("x_matrix.txt", x_matrix);
    writeMatrixToFile("y_matrix.txt", y_matrix);
    return 0;
}

