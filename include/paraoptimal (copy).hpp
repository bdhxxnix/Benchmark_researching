
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
#include <iostream>
#include <chrono>
#include <fstream>

#ifdef _OPENMP
#include <omp.h>
#else
// #pragma message ("Compilation with -fopenmp is optional but recommended")
#define omp_get_num_procs() 1
#define omp_get_max_threads() 1
#endif

namespace PGM_P::internal {

template<typename T>
using LargeSigned = typename std::conditional_t<std::is_floating_point_v<T>,
                                                long double,
                                                std::conditional_t<(sizeof(T) < 8), int64_t, __int128>>;

template<typename X, typename Y>
class ParaOptimalLinearModel {
private:
    using SX = LargeSigned<X>;
    using SY = LargeSigned<Y>;

    X first_x = 0;
    X last_x = 0;
    
    const Y epsilon;

    struct FeasibleRegion {
        double slope_min, slope_max;
        double intercept_min, intercept_max;

        FeasibleRegion(){
            slope_min = std::numeric_limits<double>::lowest();
            slope_max = std::numeric_limits<double>::max();
            intercept_min = std::numeric_limits<double>::lowest();
            intercept_max = std::numeric_limits<double>::max(); 
        }

        void Reset(){
            slope_min = std::numeric_limits<double>::lowest();
            slope_max = std::numeric_limits<double>::max();
            intercept_min = std::numeric_limits<double>::lowest();
            intercept_max = std::numeric_limits<double>::max(); 
        }
        
    };

    FeasibleRegion feasible_region;

public:

    class CanonicalSegment;

    explicit ParaOptimalLinearModel(Y epsilon) : epsilon(epsilon){
        if (epsilon < 0)
            throw std::invalid_argument("epsilon cannot be negative");
        feasible_region.Reset();
        first_x = 0;
    }

    bool add_point(const X &x, const Y &y) {
        last_x = x;
        if(first_x == 0){
            first_x = x;
        }
        double new_b_min = std::numeric_limits<double>::lowest();
        double new_b_max = std::numeric_limits<double>::max(); 
        // To maintain a convex region, we update the b-boundaries
        if (x != 0) {
            double b1_min = (y - epsilon) - feasible_region.slope_max * x;
            double b1_max = (y + epsilon) - feasible_region.slope_min * x;

            double b2_min = (y - epsilon) - feasible_region.slope_min * x;
            double b2_max = (y + epsilon) - feasible_region.slope_max * x;

            new_b_min = std::max(feasible_region.intercept_min, std::min(b1_min, b2_min));
            new_b_max = std::min(feasible_region.intercept_max, std::max(b1_max, b2_max));
        } else {
            // If xi = 0, then constraint is purely on b
            double b_low = y - epsilon;
            double b_high = y + epsilon;
            new_b_min = std::max(feasible_region.intercept_min, b_low);
            new_b_max = std::min(feasible_region.intercept_max, b_high);
        }

        if (new_b_min > new_b_max){
            first_x = x;
            feasible_region.Reset();
            return false;
        }

        feasible_region.intercept_min = new_b_min;
        feasible_region.intercept_max = new_b_max;
        return true;
    }


    CanonicalSegment get_segment() {
        return CanonicalSegment(feasible_region, first_x);
    }

    
};

template<typename X, typename Y>
class ParaOptimalLinearModel<X, Y>::CanonicalSegment {
    friend class ParaOptimalLinearModel;

    X first;
    double slope ;
    double intercept;
    FeasibleRegion region;

        CanonicalSegment(const FeasibleRegion region, X first)
        : region(region),first(first) {};

public:
    
    struct Segment{
        Segment(const double &slope, const double &intercept, const X &first_x) : slope(slope), intercept(intercept), first_x(first_x){};
        Segment() = default;
        public:
        double slope;
        double intercept;
        X first_x;
    };
    

    CanonicalSegment() = default;

    X get_first_x() const { return first; }

    Segment get_segment() const {
        
        double slope = (region.slope_min + region.slope_max)/2.0;
        double intercept = (region.intercept_min + region.intercept_max)/2.0;

        Segment segment(slope, intercept, first);
        return segment;
    }

    std::pair<double, double> get_slope_range() const {

        return {region.slope_min, region.slope_max};
    }
};

template<typename Fin, typename Fout>
size_t make_segmentation(size_t n, size_t start, size_t end, size_t epsilon, Fin in, Fout out) {
    using K = typename std::invoke_result_t<Fin, size_t>;
    size_t c = 0;
    ParaOptimalLinearModel<K, int> opt(epsilon);
    auto add_point = [&](K x, int y) {
        if (!opt.add_point(x, y)) {
            out(opt.get_segment());
            opt.add_point(x, y);
            ++c;
        }
    };

    add_point(in(start), start);
    for (size_t i = start + 1; i < end - 1; ++i) {
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
        add_point(in(end - 1), end - 1);

    if (end == n) {
        // Ensure values greater than the last one are mapped to n
        if constexpr (std::is_floating_point_v<K>) {
            add_point(std::nextafter(in(n - 1), std::numeric_limits<K>::infinity()), n);
        } else {
            add_point(in(n - 1) + 1, n);
        }
    }

    out(opt.get_segment());
    return ++c;
}

template<typename Fin, typename Fout>
size_t make_segmentation(size_t n, size_t epsilon, Fin in, Fout out) {
    return make_segmentation(n, 0, n, epsilon, in, out);
}

template<typename Fin, typename Fout>
size_t make_segmentation_par(size_t n, size_t epsilon, Fin in, Fout out, int parallelism = 16) {
    
    // printf("The number of threads is %d\n", parallelism);
    auto chunk_size = n / parallelism;
    auto c = 0ull;

    if (parallelism == 1 || n < 1ull << 15)
        return make_segmentation(n, epsilon, in, out);

    using K = typename std::invoke_result_t<Fin, size_t>;
    using canonical_segment = typename ParaOptimalLinearModel<K, int>::CanonicalSegment;
    std::vector<std::vector<canonical_segment>> results(parallelism);

    #pragma omp parallel for reduction(+:c) num_threads(parallelism)
    for (auto i = 0; i < parallelism; ++i) {
        auto first = i * chunk_size;
        auto last = i == parallelism - 1 ? n : first + chunk_size;
        if (first > 0) {
            for (; first < last; ++first)
                if (in(first) != in(first - 1))
                    break;
            if (first == last)
                continue;
        }

        auto in_fun = [in](auto j) { return in(j); };
        auto out_fun = [&results, i](const auto &cs) { results[i].emplace_back(cs); };
        results[i].reserve(chunk_size / (epsilon > 0 ? epsilon * epsilon : 16));
        c += make_segmentation(n, first, last, epsilon, in_fun, out_fun);
    }

    for (auto &v : results)
        for (auto &cs : v)
            out(cs);

    return c;
}


template<typename SegmentType>
void checkForEpsilon(std::vector<double> data, std::vector<SegmentType> segments, int first = 0, int last = 2){
    
    auto start = segments[first].first_x;
    auto end = segments[last].first_x;
    
    int segmnt_idx = first;
    auto seg = segments[segmnt_idx];
    auto slope = seg.slope;
    auto intercept = seg.intercept;
    auto first_x = seg.first_x;
    auto last_x = segments[segmnt_idx+1].first_x;

    int i = 0;
    auto max_residual = std::numeric_limits<double>::min();
    while (data[i] <= end)
    {
        if (data[i] == last_x)
        {
            // printout the max_residual
            
            printf("The %dth segment is (%f,%f) \n",segmnt_idx,first_x,last_x);
            printf("The slope is %f and the intercept is %f\n",slope,intercept);
                
            printf("The max_residual for Segment %d is %f\n",segmnt_idx,max_residual);
            max_residual = std::numeric_limits<double>::min();
            
            segmnt_idx++;
            if(segmnt_idx>=segments.size()){
                break;
            }
            seg = segments[segmnt_idx];
            slope = seg.slope;
            intercept = seg.intercept;
            first_x = seg.first_x;
            last_x = segments[segmnt_idx+1].first_x;
        }
        printf("The real position is %d and the approximate position is %f\n",i,(slope*data[i] + intercept));
        printf("The residual is %f for %dth element\n",std::abs(i - (slope*data[i] + intercept)),i);
        auto residual = std::abs(i - (slope*data[i] + intercept));
        if(residual>max_residual){
            max_residual = residual;
        }
        i++;
    }
}

}