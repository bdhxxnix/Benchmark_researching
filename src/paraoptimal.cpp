#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <cassert>

namespace PGM_P::internal {

struct Point {
    double a, b;
    Point(double a_ = 0.0, double b_ = 0.0) : a(a_), b(b_) {}
};

// Evaluate line at x
inline double eval(const Point& p, double x) {
    return p.a * x + p.b;
}

// Line intersection for two lines represented as: b = -x * a + y +/- delta
Point intersect(double x, double y, double delta, const Point& p1, const Point& p2) {
    // Line1: b = -x * a + y + delta
    // Line2: p1 and p2 -> line in (a,b) space
    double A1 = x, B1 = 1.0, C1 = y + delta; // x * a + b = y + delta
    double A2 = p2.a - p1.a;
    double B2 = p2.b - p1.b;
    double C2 = A2 * p1.a + B2 * p1.b;
    double det = A1 * B2 - A2 * B1;
    if (std::abs(det) < 1e-9) return Point(); // Parallel lines
    double a = (C1 * B2 - C2 * B1) / det;
    double b = (A1 * C2 - A2 * C1) / det;
    return Point(a, b);
}

// Feasible region defined by convex hull points
class FeasibleRegion {
public:
    std::vector<Point> poly;

    FeasibleRegion(const Point& pl, const Point& p1, const Point& pr, const Point& p1_dash) {
        poly = {pl, p1, pr, p1_dash};
    }

    // Clip by upper constraint: b = -x*a + y + delta
    void clipUpper(double x, double y, double delta) {
        std::vector<Point> new_poly;
        for (size_t i = 0; i < poly.size(); ++i) {
            Point cur = poly[i];
            Point nxt = poly[(i + 1) % poly.size()];
            bool cur_inside = eval(cur, x) <= y + delta;
            bool nxt_inside = eval(nxt, x) <= y + delta;

            if (cur_inside && nxt_inside) {
                new_poly.push_back(nxt);
            } else if (cur_inside && !nxt_inside) {
                new_poly.push_back(intersect(x, y, delta, cur, nxt));
            } else if (!cur_inside && nxt_inside) {
                new_poly.push_back(intersect(x, y, delta, cur, nxt));
                new_poly.push_back(nxt);
            }
        }
        poly = new_poly;
    }

    // Clip by lower constraint: b = -x*a + y - delta
    void clipLower(double x, double y, double delta) {
        std::vector<Point> new_poly;
        for (size_t i = 0; i < poly.size(); ++i) {
            Point cur = poly[i];
            Point nxt = poly[(i + 1) % poly.size()];
            bool cur_inside = eval(cur, x) >= y - delta;
            bool nxt_inside = eval(nxt, x) >= y - delta;

            if (cur_inside && nxt_inside) {
                new_poly.push_back(nxt);
            } else if (cur_inside && !nxt_inside) {
                new_poly.push_back(intersect(x, y, -delta, cur, nxt));
            } else if (!cur_inside && nxt_inside) {
                new_poly.push_back(intersect(x, y, -delta, cur, nxt));
                new_poly.push_back(nxt);
            }
        }
        poly = new_poly;
    }

    Point getAnyFeasible() const {
        assert(!poly.empty());
        return poly[0];
    }
};

// ParaOptimal algorithm
Point ParaOptimal(const std::vector<std::pair<double, double>>& S, double delta) {
    assert(S.size() >= 2);
    auto [x1, y1] = S[0];
    auto [x2, y2] = S[1];
    double denom = x2 - x1;
    assert(std::abs(denom) > 1e-9);

    Point pl((y2 - y1 - 2 * delta) / denom, (x2 * y1 - x1 * y2 + x2 * delta + x1 * delta) / denom);
    Point pr((y2 - y1 + 2 * delta) / denom, (x2 * y1 - x1 * y2 - x2 * delta - x1 * delta) / denom);
    Point p1((y2 - y1) / denom, (x2 * y1 - x1 * y2 + x2 * delta - x1 * delta) / denom);
    Point p1_dash((y2 - y1) / denom, (x2 * y1 - x1 * y2 - x2 * delta + x1 * delta) / denom);

    FeasibleRegion region(pl, p1, pr, p1_dash);

    for (size_t i = 2; i < S.size(); ++i) {
        double x = S[i].first;
        double y = S[i].second;

        if (eval(pl, x) <= y + delta) {
            if (eval(pr, x) >= y + delta) {
                region.clipUpper(x, y, delta);
                region.clipLower(x, y, delta);
            }
        } else {
            break; // cannot extend further
        }
    }

    return region.getAnyFeasible();
}

}