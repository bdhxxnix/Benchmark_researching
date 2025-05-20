
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#else
#pragma message("Compilation with -fopenmp is optional but recommended")
#define omp_get_num_procs() 8
#define omp_get_max_threads() 8
#endif


template <typename T>
using LargeSigned = typename std::conditional_t<std::is_floating_point_v<T>,
                                                long double,
                                                std::conditional_t<(sizeof(T) < 8), int64_t, __int128>>;

namespace PGM_S::internal{
template <typename X, typename Y>
class SwingPiecewiseLinearModel
{
private:
    using SX = LargeSigned<X>;
    using SY = LargeSigned<Y>;

    struct Point
    {
        X x{};
        Y y{};

        // Get the slope between self and p
        double operator-(const Point &p) const
        {
            if (p.x == x)
            {
                throw std::invalid_argument("can not handle two same points");
            }
            double dx = static_cast<double>(x) - static_cast<double>(p.x);
            double dy = static_cast<double>(y) - static_cast<double>(p.y);
            
            double slope = dy / dx;
            return slope;
        };
    };

    X first_x = 0;
    X last_x = 0;
    Point end[2];
    int points = 0;

    const Y epsilon;
    long double rho_max, rho_min;
    Point so;

    auto cross(const Point &O, const Point &A, const Point &B) const
    {
        auto OA = A - O;
        auto OB = B - O;
        return OA.dx * OB.dy - OA.dy * OB.dx;
    }

public:
    class CanonicalSegment;

    explicit SwingPiecewiseLinearModel(Y epsilon) : epsilon(epsilon)
    {
        if (epsilon < 0)
            throw std::invalid_argument("epsilon cannot be negative");
    }

    bool add_point(const X &x, const Y &y)
    {
        
        last_x = x;
        Point p1{x, y + epsilon};
        Point p2{x, y - epsilon};

        if (points == 0)
        {
            // It is the first point of the segment
            first_x = x;
            so = {x, y};
            points++;
            return true;
        }

        if (points == 1)
        {
            end[0] = p1;
            end[1] = p2;
            points++;
            rho_max = so - end[0];
            rho_min = so - end[1];
            return true;
        }

        double slope1 = so - p2; // lower bound
        double slope2 = so - p1; // upper bound
        bool outside_line1 = slope1 > rho_max;
        bool outside_line2 = slope2 < rho_min;

        if (outside_line1 || outside_line2)
        {
            points = 0;
            return false;
        }

        if (slope2 < rho_max)
        {
            end[0] = p1;
            rho_max = slope2;
        }

        if (slope1 > rho_min)
        {
            end[1] = p2;
            rho_min = slope1;
        }
        return true;
    }

    CanonicalSegment get_segment()
    {
        if (points == 1)
        {
            return CanonicalSegment(so, first_x, last_x);
        }
        return CanonicalSegment(so, end[0], end[1], first_x, last_x);
    }

    void reset()
    {
        points = 0;
    }
};

template <typename X, typename Y>
class SwingPiecewiseLinearModel<X, Y>::CanonicalSegment
{
public:
    CanonicalSegment() = default;

    friend class SwingPiecewiseLinearModel;

    Point so;
    Point end[2];
    X first_x;
    X last_x;

    struct Segment{
        Segment(const long double &slope, const long double &intercept, const X &first_x, const X &last_x) : slope(slope), intercept(intercept), first_x(first_x), last_x(last_x) {};
        Segment() = default;
        public:
        long double slope;
        long double intercept;
        X first_x;
        X last_x;
    };

    CanonicalSegment(Point so, X first_x, X last_x) : so(so), first_x(first_x), last_x(last_x) {
        end[0] = so;
        end[1] = so;
    };

    CanonicalSegment(Point so, const Point &p1, const Point &p2, X first_x, X last_x)
        : so(so), end{p1, p2}, first_x(first_x), last_x(last_x) {};

    bool one_point() const
    {
        return so.x == end[0].x;
    }

    Segment get_Canonicalsegment() const
    {
        Point p1 = this->end[0];
        Point p2 = this->end[1];
        

        if (one_point())
            return {0, static_cast<long double>(end[0].y + end[1].y) / 2, first_x, last_x};

        auto [i_x, i_y] = so;

        
        auto rho_max = (p1.y - i_y) / (static_cast<double>(p1.x) - i_x);
        auto rho_min = (p2.y - i_y) / (static_cast<double>(p2.x) - i_x);

        long double slope = (rho_max + rho_min) / 2.;
        long double intercept = i_y - static_cast<double>(i_x) * slope;
        Segment result =Segment(slope, intercept, first_x, last_x);
        return result;
    }

		
};

template <typename Fin, typename Fout>
size_t make_segmentation(size_t n, size_t start, size_t end, double epsilon, Fin in, Fout out)
{
    using K = typename std::invoke_result_t<Fin, size_t>;
    size_t c = 0;
    SwingPiecewiseLinearModel<K, int> opt(epsilon);
    auto add_point = [&](K x, int y)
    {
        
        if (!opt.add_point(x, y))
        {
            out(opt.get_segment());
            opt.add_point(x, y);
            ++c;
        }
    };

    add_point(in(start),start);
    for (size_t i = start + 1; i < end - 1; ++i)
    {
        if (in(i) == in(i - 1)) {
            // Here there is an adjustment for inputs with duplicate keys: at the end of a run of duplicate keys equal
            // to x=in(i) such that x+1!=in(i+1), we map the values x+1,...,in(i+1)-1 to their correct rank i.
            // For floating-point keys, the value x+1 above is replaced by the next representable value following x.
            if constexpr (std::is_floating_point_v<K>) {
                K next;
                if ((next = std::nextafter(in(i), std::numeric_limits<K>::infinity())) < in(i + 1))
                    add_point(next, i);
            } else {
                if (in(i) + 1 < in(i + 1))
                    add_point(in(i) + 1, i);
            }
        } else {
            add_point(in(i), i);
        }
    }
    if (end >= start + 2 && in(end - 1) != in(end - 2))
        add_point(in(end - 1),(end-1));

    if (end == n)
    {
        add_point(in(n - 1) + 1 ,n);
    }

    out(opt.get_segment());
    return ++c;
}

template <typename Fin, typename Fout>
size_t make_segmentation(size_t n, double epsilon, Fin in, Fout out)
{
    return make_segmentation(n, 0, n, epsilon, in, out);
}

template <typename Fin, typename Fout>
size_t make_segmentation_par(size_t n, size_t epsilon, Fin in, Fout out,int parallelism = 16)
{
    
    // printf("The number of threads is %d\n", parallelism);
    auto chunk_size = n / parallelism;
    auto c = 0ull;

    if (parallelism == 1 || n < 1ull << 15)
        return make_segmentation(n, epsilon, in, out);

    using K = typename std::invoke_result_t<Fin, size_t>;
    using canonical_segment = typename SwingPiecewiseLinearModel<K, int>::CanonicalSegment;
    std::vector<std::vector<canonical_segment>> results(parallelism);

#pragma omp parallel for reduction(+ : c) num_threads(parallelism)
    for (auto i = 0; i < parallelism; ++i)
    {
        auto first = i * chunk_size;
        auto last = i == parallelism - 1 ? n : first + chunk_size;
        if (first > 0)
        {
            for (; first < last; ++first)
                if (in(first) != in(first - 1))
                    break;
            if (first == last)
                continue;
        }

        auto in_fun = [in](auto j){ return in(j); };
        auto out_fun = [&results, i](const auto &cs){ results[i].emplace_back(cs); };
        results[i].reserve(chunk_size / (epsilon > 0 ? epsilon * epsilon : 16));
        c += make_segmentation(n, first, last, epsilon, in_fun, out_fun);
    }

    for (auto &v : results)
        for (auto &cs : v)
            out(cs);

    return c;
}

/**
 * @brief A method to check the rightness of the segmentation of greedyPLR 
 * 
 * @tparam SegmentType 
 * @param data 
 * @param segments 
 */
template <typename Fin>
void checkForEpsilon(Fin data, auto segments,int begin = 0,int ed = 1, double epsilon = 0.1){
    
    auto start = segments[begin].first_x;
    auto end = segments[ed].last_x;
    
    int segmnt_idx = begin;
    auto seg = segments[segmnt_idx];
    auto slope = seg.slope;
    auto intercept = seg.intercept;
    auto first_x = seg.first_x;
    auto last_x = seg.last_x;

    int i = 0;
    auto max_residual = std::numeric_limits<double>::min();
    while (data(i) <= end)
    {
        if (data(i) == last_x)
        {
            // printout the max_residual
            if( i <= 1000){
                printf("The slope is %Lf and the intercept is %Lf\n",slope,intercept);
                printf("The first_x is %f and the last_x is %f\n",first_x,last_x);
            }
            if(max_residual > epsilon + 1){
                printf("The max_residual for Segment %d is %f\n",segmnt_idx,max_residual);
            }
            max_residual = std::numeric_limits<double>::min();
           
            segmnt_idx++;
            if(segmnt_idx>=segments.size()){
                break;
            }
            seg = segments[segmnt_idx];
            slope = seg.slope;
            intercept = seg.intercept;
            first_x = seg.first_x;
            last_x = seg.last_x;
        }
        auto residual = std::abs(i - (slope*data(i) + intercept));
        if(residual>max_residual){
            max_residual = residual;
        }
        i++;
    }
    printf("\n");
    }
}
