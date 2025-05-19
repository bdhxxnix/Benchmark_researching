
/* This is the v2.0 of optimalPLR foundation, specially designed for PGM-Index */
/**
 * @file OptimalPLR_Foundation.hpp
 * @author Jiayong Qin (2057729401@qq.com)
 * @brief
 * @version 2.0
 * @date 2025-2-28
 *
 * @copyright Copyright (c) 2025
 *
 */
#include <iostream>
#include <vector>
#include <deque>
#include <limits>
#include <algorithm>
#include <fstream>

#ifdef _OPENMP
#include <omp.h>
#else
// #pragma message ("Compilation with -fopenmp is optional but recommended")
#define omp_get_num_procs() 1
#define omp_get_max_threads() 1
#define omp_get_thread_num() 0
#endif

#include <thread>

#include <chrono>


namespace PGM_C2::internal
{

template<typename T>
using LargeSigned = typename std::conditional_t<std::is_floating_point_v<T>,
                                                long double,
                                                std::conditional_t<(sizeof(T) < 8), int64_t, __int128>>;

 
template<typename K> 
struct Segment{
    /**
     * @brief Segments
     * @result the slope and intercept of the line
     * @result the start element of the segment
     */
    double slope, intercept;
    K first_x;

    // Set the start_idx of the segment as key
    K key() const { return first_x; }
    Segment( double slope = 0,  double intercept = 0, K first_x = K(0)) : slope(slope), intercept(intercept), first_x(first_x)  {};
    
    friend inline bool operator<(const Segment &s, const K &k) { return s.key() < k; }
    friend inline bool operator<(const K &k, const Segment &s) { return k < s.key(); }
    friend inline bool operator<(const Segment &s, const Segment &t) { return s.key() < t.key(); }

    operator K() { return key(); };

    /**
     * Returns the approximate position of the specified key.
     * @param k the key whose position must be approximated
     * @return the approximate position of the specified key
     */
    inline double operator()(const K &k) const {
        
        return (slope * k + intercept);
    }
};

template <typename X, typename Y>
class OptimalPLR{
  public:
    /**
     * @brief Construct a new OptimalPLR object
     *
     * @param error error bounds of this implementation
     */
    OptimalPLR(int epsilon) : epsilon(epsilon) {}

    using Segment = PGM_C2::internal::Segment<X>;

    std::vector<Segment> segmentData(const std::vector<X> &data)
    {
        std::vector<Segment> segments;
        if (data.empty())
            return segments;

        // Initialization: the first two points acting as the initial hull and segement
        size_t start = 0;
        Point sa = {data[0], 0 + epsilon};
        Point sc = {data[1], 1 - epsilon};
        Point sb = {data[0], 0 - epsilon};
        Point sd = {data[1], 1 + epsilon};


        double rho_max = computeSlope(sb, sd);
        double rho_min = computeSlope(sa, sc);

        std::deque<Point> upperHull = {sa, sd};
        std::deque<Point> lowerHull = {sb, sc};

        for (int i = 2; i < data.size(); ++i)
        {
            Point s = {data[i], i};

            if (!isWithinBounds(s, rho_max, sc, sd, rho_min))
            {
                // Check if this point is within the bounds , if not , we need to update the segement and initialize again
                std::pair<double,double> so = findIntersection(sa, sc, sb, sd);
                double finalSlope = (rho_min + rho_max) / 2;
                double intercept = computeIntercept_for_intersection(so, finalSlope);
                segments.push_back({finalSlope,intercept,data[start]});

                start = i;
                sa = {data[start], i + epsilon};
                sc = {data[start + 1], i+1 - epsilon};
                sb = {data[start], i - epsilon};
                sd = {data[start + 1], i+1 + epsilon};

                rho_max = computeSlope(sb, sd);
                rho_min = computeSlope(sa, sc);
                upperHull.clear();
                lowerHull.clear();
                upperHull.push_back(sa);
                upperHull.push_back(sd);
                lowerHull.push_back(sb);
                lowerHull.push_back(sc);
                i++;

                continue;
            }
            else
            {
                Point s_upper = {s.x, s.y + epsilon};
                Point s_lower = {s.x, s.y - epsilon};
                bool upper = 0;
                bool lower = 0;
                auto slope1 = computeSlope(s_lower , sa);
                auto slope2 = computeSlope(s_upper , sb);

                if (slope1 > rho_min)
                {
                    // Find the max slope point and delete the points before it so that we can update the min slope
                    lower = 1;
                    sa = findMaxSlopePoint(upperHull, s_lower);
                    sc = s_lower;
                    deletePointsBefore(upperHull, sa);
                    rho_min = computeSlope(sa, sc);
                }

                if (slope2 < rho_max)
                {
                    // Find the min slope point and delete the points before it so that we can update the segement and the hull
                    upper = 1;
                    sb = findMinSlopePoint(lowerHull, s_upper);
                    sd = s_upper;
                    deletePointsBefore(lowerHull, sb);
                    rho_max = computeSlope(sb, sd);
                }
                if (upper)
                {
                    updateConvexHull(upperHull, s_upper, true);
                    sd = upperHull.back();
                }
                if (lower)
                {
                    updateConvexHull(lowerHull, s_lower, false);
                    sc = lowerHull.back();
                }
            }
        }

        std::pair<double,double> so = findIntersection(sa, sc, sb, sd);
        double finalSlope = (rho_min + rho_max) / 2;
        segments.push_back({finalSlope, computeIntercept_for_intersection(so, finalSlope), data[start]});

        return segments;
    }

private:
    int epsilon;

    using SX = LargeSigned<X>;
    using SY = LargeSigned<Y>;

    struct Slope{
        SX dx{};
        SY dy{};

        bool operator<(const Slope &p) const { return dy * p.dx < dx * p.dy; }
        bool operator>(const Slope &p) const { return dy * p.dx > dx * p.dy; }
        bool operator==(const Slope &p) const { return dy * p.dx == dx * p.dy; }
        bool operator!=(const Slope &p) const { return dy * p.dx != dx * p.dy; }
        explicit operator double() const { return dy / (double) dx; }
 
    };

    struct Point{
        SX x{};
        SY y{};
        Slope operator-(const Point &p) const { return {SX(x) - p.x, SY(y) - p.y}; }
 
    };

    /**
     * @brief compute the slope between two points
     *
     * @param p1
     * @param p2
     * @return double
     */
    double computeSlope(const Point &p1, const Point &p2)
    {
        if (p2.x == p1.x)
            return std::numeric_limits<double>::infinity();
       
        double slope = static_cast<double>(p1 - p2);
        return slope;
    }

    /**
     * @brief compute the intercept of a line
     *
     * @param p
     * @param slope
     * @return double
     */
    double computeIntercept(const Point &p, double slope)
    {
        return double(p.y) - slope * p.x;
    }
    // Specialized for intersection point 
    double computeIntercept_for_intersection(std::pair<double,double> &p, double slope)
    {
        return p.second - slope * p.first;
    }

    /**
     * @brief To Check if the point is within the bounds
     *
     * @param s
     * @param rho_max
     * @param sa
     * @param sc
     * @param rho_min
     * @param sb
     * @param sd
     * @return true
     * @return false
     */
    bool isWithinBounds(const Point &s, double rho_max, Point sc, Point sd, double rho_min)
    {
        return (s.y + epsilon >= rho_min * (s.x- sc.x) + sc.y) && (s.y- epsilon <= rho_max * (s.x- sd.x) + sd.y);
    }

    /**
     * @brief update the convex hull
     *
     * @param hull
     * @param s
     * @param upper if the hull is upper or not
     */
    void updateConvexHull(std::deque<Point> &hull, const Point &s, bool upper)
    {
        while (hull.size() >= 2)
        {
            auto &p1 = hull[hull.size() - 2];
            auto &p2 = hull[hull.size() - 1];
            if (isRedundant(p1, p2, s, upper))
            {
                hull.pop_back();
            }
            else
            {
                break;
            }
        }
        hull.push_back(s);
    }

    /**
     * @brief delete the redundant points using triangle area method
     *
     * @param p1
     * @param p2
     * @param p3
     * @param upper
     * @return true
     * @return false
     */
    bool isRedundant(const Point &p1, const Point &p2, const Point &p3, bool upper)
    {
        double slope1 = static_cast<double>(p2 - p1);
        double slope2 = static_cast<double>(p3 - p2);
        return upper ? (slope1 >= slope2) : (slope1 <= slope2);
    }

    /**
     * @brief Find the min slope point
     *
     * @param hull
     * @param s
     * @return Point
     */
    Point findMinSlopePoint(std::deque<Point> &hull, const Point &s)
    {

        return *std::min_element(hull.begin(), hull.end(), [&](const Point &a, const Point &b)
                                 { return computeSlope(a, s) < computeSlope(b, s); });
    }

    /**
     * @brief find the max slope point
     *
     * @param hull
     * @param s
     * @return Point
     */
    Point findMaxSlopePoint(std::deque<Point> &hull, const Point &s)
    {
        return *std::max_element(hull.begin(), hull.end(), [&](const Point &a, const Point &b)
                                 { return computeSlope(a, s) < computeSlope(b, s); });
    }

    /**
     * @brief delete point before the reference point
     *
     * @param hull
     * @param reference
     */
    void deletePointsBefore(std::deque<Point> &hull, const Point &reference)
    {
        while (!hull.empty() && hull[0].x < reference.x)
        {
            hull.pop_front();
        }
    }

    /**
     * @brief Find the intersection point of two lines
     *
     * @param sa
     * @param sc
     * @param sb
     * @param sd
     * @return Point
     */
    std::pair<double,double> findIntersection(Point sa, Point sc, Point sb, Point sd)
    {
        double a1 = computeSlope(sa, sc);
        double b1 = computeIntercept(sa, a1);
        double a2 = computeSlope(sb, sd);
        double b2 = computeIntercept(sb, a2);

        if (a1 == a2)      // Parallel lines
            return {0, 0}; // or handle as needed
        double x_intersect = (double)(b2 - b1) / (double)(a1 - a2);
        double y_intersect = a1 * x_intersect + b1;

        return {x_intersect, y_intersect};
    }

};


/**
 * @brief Make segmentation of the data
 *
 * @tparam Fin the input data type
 * @tparam Fout the output data type
 * @param n the size of the data
 * @param epsilon error boundary
 * @param in input data
 * @param out output data
 * @return size_t the number of the segments
 */
template <typename X, typename Y>
size_t make_segmentation(int n, int epsilon, 
                            std::vector<X> in, std::vector<Segment<X>> &out)
{
    PGM_C2::internal::OptimalPLR<X,Y> opt(epsilon);
    auto segments = opt.segmentData(in);
    out.insert(out.end(), segments.begin(), segments.end());
    return segments.size();
}
/**
 * @brief Translate the segments into points
 * 
 * @param segments 
 * @param first 
 * @param last 
 * @return std::vector<X> 
 */
template <typename X>
auto translate(std::vector<Segment<X>> segments, int first, int last){
    std::vector<X> data;
    for(int i = first; i < last; i++){
        data.push_back(segments[i].first_x);
    }
    return data;
}
/**
 * @brief Make segmentation of the data parallelly
 *
 * @tparam Fin the input data type
 * @tparam Fout the output data type
 * @param n the size of the data
 * @param epsilon error boundary
 * @param in input data
 * @param out output data
 * @return size_t the number of the segments
 */
template <typename X,typename Y>
size_t make_segmentation_par(int n, int epsilon, 
                            std::vector<X> in, std::vector<Segment<X>> &out,
                            int parallelism =16)
{
    
    // printf("The number of threads is %d\n", parallelism);

    int chunk_size = std::max(n / parallelism, 1000);  // Prevent too many threads for small `n`

    if (parallelism == 1 || n < 1ull << 15)
        return make_segmentation<X,Y>(n, epsilon, in, out);
    OptimalPLR<X,Y> opt(epsilon);
    std::vector<std::vector<Segment<X>>> local_segments(parallelism);

    size_t total_segments = 0;

    #pragma omp parallel for num_threads(parallelism)
    for (int i = 0; i < parallelism; i++)
    {
        size_t start = i * chunk_size;
        size_t end = (i == parallelism - 1) ? n : (i + 1) * chunk_size;

        if (start < n) {
            auto segments = opt.segmentData({in.begin() + start, in.begin() + end});
            local_segments[i] = std::move(segments);

            #pragma omp atomic
            total_segments += local_segments[i].size();
        }
    }

    for (size_t i = 0; i < local_segments.size(); i++) {
        out.insert(out.end(), local_segments[i].begin(), local_segments[i].end());
    }

    return total_segments;
}

/**
 * @brief Check the rightness of the segmentation of OptimalPLR
 * 
 * @tparam SegmentType 
 * @param data 
 * @param segments 
 * @param first 
 * @param last 
 */
template <typename K>
void checkForEpsilon(std::vector<K> data, std::vector<Segment<K>> segments, int first = 0, int last = 5,double epsilon = 2){
    
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
            
            if(max_residual > epsilon + 1){
                printf("The max_residual for Segment %d is %f\n",segmnt_idx,max_residual);
                printf("The %dth segment is (%f,%f) \n",segmnt_idx,first_x,last_x);
                printf("The slope is %f and the intercept is %f\n",slope,intercept);
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
            last_x = segments[segmnt_idx+1].first_x;
        }
        // printf("The real position is %d and the approximate position is %f\n",i,(slope*data[i] + intercept));
        // printf("The residual is %f for %dth element\n",std::abs(i - (slope*data[i] + intercept)),i);
        auto residual = std::abs(i - (slope*data[i] + intercept));
        if(residual>max_residual){
            max_residual = residual;
        }
        i++;
    }
}
}