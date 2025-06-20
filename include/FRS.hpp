#pragma once

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


namespace FRS::internal {
template <typename T>
using LargeSigned = typename std::conditional_t<std::is_floating_point_v<T>,
                                                long double,
                                                std::conditional_t<(sizeof(T) < 8), int64_t, __int128>>;

template <typename X, typename Y>
class FRSPiecewiseLinearModel
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

    X first_x = 0;  // First x value of the segment
    X last_x = 0;   // Last x value of the segment
    int points = 0; // Number of points in the segment, used to check if this is the first point or not

    const Y epsilon;    // Error bound for the PLA method
    const X upperbound; // Upper bound for the x values, used to calculate the fixed slope
    const Y last_y;     // Last y value, used to calculate the fixed slope
    long double slope;  // Slope of the segment
    long double intercept;// Intercept of the segment

public:
    class CanonicalSegment;

    explicit FRSPiecewiseLinearModel(Y epsilon, X upperbound, Y last_y) 
    : epsilon(epsilon),upperbound(upperbound), last_y(last_y)
    {
        if (epsilon < 0 || upperbound < 0)
            throw std::invalid_argument("epsilon or upperbound cannot be negative");
    }

    bool add_point(const X &x, const Y &y)
    {
        
        last_x = x;
        Point thisP = {x,y};

        if (points == 0)
        {
            // It is the first point of the segment
            first_x = x;
            Point finalPoint = {upperbound,last_y};
            slope = finalPoint - thisP;
            intercept = static_cast<double>(y) - slope * x;
            points++;
            return true;
        }

        long double pre = slope * x + intercept;
        long double residual = std::abs(pre - static_cast<double>(y));

        if (residual > epsilon)
        {
            points = 0;
            return false;
        }
        return true;
    }

    CanonicalSegment get_segment()
    {
        return CanonicalSegment(slope,intercept,first_x, last_x);
    }

    void reset()
    {
        points = 0;
    }
};

template <typename X, typename Y>
class FRSPiecewiseLinearModel<X, Y>::CanonicalSegment
{
public:
    CanonicalSegment() = default;

    friend class FRSPiecewiseLinearModel;

    long double slope ;
    long double intercept;
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

    CanonicalSegment(long double slope, long double intercept, X first_x, X last_x) 
    : slope(slope), intercept(intercept),first_x(first_x), last_x(last_x) {};

    Segment get_Canonicalsegment() const
    {
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
size_t make_segmentation(size_t n, size_t start, size_t end, size_t epsilon, Fin in, Fout out)
{

    using K = typename std::invoke_result_t<Fin, size_t>;
    size_t c = 0;
    FRSPiecewiseLinearModel<K, int> opt(epsilon,in(n-1),n-1);
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

    out(opt.get_segment());
    return ++c;
}

// Overloaded convenience wrapper for the full range of the input data
template <typename Fin, typename Fout>
size_t make_segmentation(size_t n, double epsilon, Fin in, Fout out)
{
    return make_segmentation(n, 0, n, epsilon, in, out);
}


/**
 * @brief A method to make a parallel segmentation of the input data using the FRS algorithm
 * @param n The number of elements in the data
 * @param epsilon The error bound for the segmentation
 * @param in The input function that takes an index and returns the value at that index
 * @param out The output function that puts the created segments in pointed vector
 * @param parallelism The number of threads to use for parallelization
 * @return The number of segments created
 */
template <typename Fin, typename Fout>
size_t make_segmentation_par(size_t n, size_t epsilon, Fin in, Fout out,int parallelism = 16)
{
    auto chunk_size = n / parallelism;
    auto c = 0ull;

    if (parallelism == 1 || n < 1ull << 15)
        return make_segmentation(n, epsilon, in, out);

    using K = typename std::invoke_result_t<Fin, size_t>;
    using canonical_segment = typename FRSPiecewiseLinearModel<K, int>::CanonicalSegment;
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
 * @tparam Fin the type of the input function
 * @tparam SegmentType the type of the segments
 * @param data the input data
 * @param segments the segments created by the segmentation algorithm that to be checked
 * @param begin the index of the first segment to check
 * @param ed the index of the last segment to check
 * @param epsilon the error bound for the segmentation
 * @return void
 */
template <typename Fin, typename SegmentType>
void checkForEpsilon(Fin data, std::vector<SegmentType> segments,int begin = 0,int ed = 1, double epsilon = 0.1){
    
    auto start = segments[begin].first_x;
    auto end = segments[ed].last_x;
    
    int segmnt_idx = begin;
    auto seg = segments[segmnt_idx];
    auto slope = seg.slope;
    auto intercept = seg.intercept;
    auto first_x = seg.first_x;
    auto last_x = seg.last_x;

    size_t i = 0;
    while(data(i)<start) i++;
    auto max_residual = std::numeric_limits<double>::min();
    while (data(i) <= end)
    {
        if (data(i) == last_x)
        {
            // if the max residual is larger than the error bound, printout the max_residual
            if(max_residual > epsilon + 1){
                printf("The max_residual for Segment %d is %f\n",segmnt_idx,max_residual);
                printf("The slope is %Lf and the intercept is %Lf\n",slope,intercept);
                printf("The first_x is %ld and the last_x is %ld\n",first_x,last_x);
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
