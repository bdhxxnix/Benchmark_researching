#pragma once
#include <iostream>
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


namespace Greedy::internal {
template <typename T>
using LargeSigned = typename std::conditional_t<std::is_floating_point_v<T>,
                                                long double,
                                                std::conditional_t<(sizeof(T) < 8), int64_t, __int128>>;

template <typename X, typename Y>
class GreedyPiecewiseLinearModel
{
private:
    using SX = LargeSigned<X>;
    using SY = LargeSigned<Y>;

    struct Point
    {
        X x{};
        Y y{};

        // Get the slope between self and p
        long double operator-(const Point &p) const
        {
            if (p.x == x)
            {
                throw std::invalid_argument("can not handle two same points");
            }
            long double dx = static_cast<long double>(x) - static_cast<long double>(p.x);
            long double dy = static_cast<long double>(y) - static_cast<long double>(p.y);
            
            long double slope = dy / dx;
            return slope;
        };
    };

    X first_x = 0;      // First x value of the segment
    X last_x = 0;       // Last x value of the segment
    Point initial[2];   // Initial points of the segment
    Point end[2];       // End points of the segment
    int points = 0;     // Number of points in the segment, used to check if this is the first point or not

    const Y epsilon;    // Error bound for the PLA method
    long double rho_max, rho_min;   // Maximum and minimum slopes of the segment  
    std::pair<long double, long double> so; // Pivot point of the segment, used to calculate the intercept

    

public:
    class CanonicalSegment;

    explicit GreedyPiecewiseLinearModel(Y epsilon) : epsilon(epsilon)
    {
        if (epsilon < 0)
            throw std::invalid_argument("epsilon cannot be negative");
    }

    bool add_point(const X &x, const Y &y)
    {
        
        
        auto max_y = std::numeric_limits<Y>::max();
        auto min_y = std::numeric_limits<Y>::lowest();
        Point p1{x, y >= max_y - epsilon ? max_y : y + epsilon};
        Point p2{x, y <= min_y + epsilon ? min_y : y - epsilon};


        if (points == 0)
        {
            // It is the first point of the segment
            first_x = x;
            initial[0] = p1;
            initial[1] = p2;
            points++;
            last_x = x;
            return true;
        }

        if (points == 1)
        {   
            // It is the second point of the segment
            end[0] = p1;
            end[1] = p2;
            points++;
            long double so_x = (static_cast<long double>(end[0].x) + initial[0].x)/2.0;
            long double so_y = (static_cast<long double>(end[0].y)+initial[0].y+end[1].y+initial[1].y)/4.0;
            rho_min = initial[0] - end[1];
            rho_max = initial[1] - end[0];
            so = {so_x, so_y};
            last_x = x;
            return true;
        }

        double slope1 = (static_cast<long double>(y -epsilon) - so.second) / (static_cast<long double>(x) - so.first);
        double slope2 = (static_cast<long double>(y +epsilon) - so.second) / (static_cast<long double>(x) - so.first);
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
        last_x = x;
        return true;
    }

    CanonicalSegment get_segment()
    {
        if (points == 1)
        {
            return CanonicalSegment(initial[0], initial[1], first_x, last_x);
        }
        return CanonicalSegment(so, end[0], end[1], first_x, last_x);
    }

    void reset()
    {
        points = 0;
    }
};

template <typename X, typename Y>
class GreedyPiecewiseLinearModel<X, Y>::CanonicalSegment
{
public:
    CanonicalSegment() = default;

    friend class GreedyPiecewiseLinearModel;

    std::pair<long double, long double> so;
    Point end[2];
    X first_x;
    X last_x;

    struct Segment{
        Segment(const long double &slope, const long double &intercept, const X &first_x, const X &last_x) : slope(slope), intercept(intercept), first_x(first_x), last_x(last_x) {};
        Segment() = default;
        public:
        long double slope;
        int64_t intercept;
        X first_x;
        X last_x;
    };

    CanonicalSegment(const Point &p1, const Point &p2, X first_x, X last_x) :end{p1, p2}, first_x(first_x), last_x(last_x) {
        so = { (p1.x + p2.x) / 2, (p1.y + p2.y) / 2 };
    };

    CanonicalSegment(std::pair<long double, long double> so, const Point &p1, const Point &p2, X first_x, X last_x)
        : so(so), end{p1, p2}, first_x(first_x), last_x(last_x) {};

    bool one_point() const
    {
        return so.first == end[0].x;
    }

    Segment get_Canonicalsegment(const X &origin) const
    {
        Point p1 = this->end[0];
        Point p2 = this->end[1];
        

        if (one_point())
            return {0, static_cast<long double>(end[0].y + end[1].y) / 2, first_x, last_x};

        auto [i_x, i_y] = so;
        
        long double rho_max = (static_cast<long double>(p1.y) - i_y) / (static_cast<long double>(p1.x) - i_x);
        long double rho_min = (static_cast<long double>(p2.y) - i_y) / (static_cast<long double>(p2.x) - i_x);

        long double slope =  (rho_max + rho_min) / 2.;
        long double raw_intercept = i_y - slope * (i_x - static_cast<double>(origin));
        int64_t intercept = static_cast<int64_t>(std::llround(raw_intercept));
        Segment result =Segment(slope, intercept, first_x, last_x);
        return result;
    }

    
		
};

/**
 * @brief a method to make a segmentation of the input data using the frs algorithm
 * @param n the number of elements in the data
 * @param start the starting index of the data
 * @param end the ending index of the data
 * @param epsilon the error bound for the segmentation
 * @param in the input function that takes an index and returns the value at that index
 * @param out the output function that puts the created segments in pointed vector
 * @return the number of segments created
 */
template <typename Fin, typename Fout>
size_t make_segmentation(size_t n, size_t start, size_t end, double epsilon, Fin in, Fout out)
{
    using K = typename std::invoke_result_t<Fin, size_t>;
    size_t c = 0;
    GreedyPiecewiseLinearModel<K, int> opt(epsilon);
    auto add_point = [&](K x, size_t y)
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

    out(opt.get_segment());
    return ++c;
}

template <typename Fin, typename Fout>
size_t make_segmentation(size_t n, double epsilon, Fin in, Fout out)
{
    return make_segmentation(n, 0, n, epsilon, in, out);
}


// The parallel version of the segmentation method
template <typename Fin, typename Fout>
size_t make_segmentation_par(size_t n, size_t epsilon, Fin in, Fout out,int parallelism = 16)
{
    
    // printf("The number of threads is %d\n", parallelism);
    auto chunk_size = n / parallelism;
    auto c = 0ull;

    if (parallelism == 1 || n < 1ull << 15)
        return make_segmentation(n, epsilon, in, out);

    using K = typename std::invoke_result_t<Fin, size_t>;
    using canonical_segment = typename GreedyPiecewiseLinearModel<K, int>::CanonicalSegment;
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
 * @brief A method to check the correctness of the segmentation of greedy PLA
 * @tparam Fin the type of the input function
 * @tparam SegmentType the type of the segments
 * @param data the input data function
 * @param segments the list of segments generated
 * @param begin the starting segment index for verification
 * @param ed the ending segment index for verification
 * @param epsilon the error bound
 */
template <typename Fin, typename SegmentType>
void checkForEpsilon(size_t n ,Fin data, std::vector<SegmentType> segments,int begin = 0,int ed = 1, size_t epsilon = 0.1){
    auto start = segments[begin].first_x;
    auto end = segments[ed].last_x;

    size_t segmnt_idx = begin;
    auto seg = segments[segmnt_idx];
    auto slope = seg.slope;
    auto intercept = seg.intercept;
    auto first_x = seg.first_x;
    auto last_x = seg.last_x;

    size_t i = 0;
    while(data(i) < start) i++;

    double max_residual = std::numeric_limits<double>::min();
    while (data(i) <= end && i < n ) {
        auto key = data(i);
        auto predicted = slope * key + intercept;
        auto residual = std::abs(i - predicted);

        if(residual > max_residual && data(i) != last_x) {
            max_residual = residual;
        }

        if (data(i) == last_x) {
            if(max_residual > epsilon + 1) {
                printf("The max_residual for Segment %ld is %f\n", segmnt_idx, max_residual);
                printf("The slope is %Lf and the intercept is %ld\n", slope, intercept);
                printf("The first_x is %ld and the last_x is %ld\n", first_x, last_x);
            }
            max_residual = std::numeric_limits<double>::min();
            segmnt_idx++;
            if(segmnt_idx >= segments.size()) break;
            seg = segments[segmnt_idx];
            slope = seg.slope;
            intercept = seg.intercept;
            first_x = seg.first_x;
            last_x = seg.last_x;
        }
        i++;
    }
}

} // namespace Greedy::internal