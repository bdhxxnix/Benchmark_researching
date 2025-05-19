/* Write your docstring here*/
/**
 * @file OptimalPLR_Foundation.hpp
 * @author Jiayong Qin (2057729401@qq.com)
 * @brief
 * @version 1.0
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


namespace PGM_C::internal
{
typedef int K;

struct Point
{
    /**
     * @brief Data used to represent a point in 2D space
     *
     */
    double x;
    double y;
    // Initialize the point
    Point(double x = 0, double y = 0) : x(x), y(y) {};
};

struct Segment
{
    /**
     * @brief Segements
     * @result the slope and intercept of the line
     * @result the start element of the segment
     */
    double slope, intercept;
    double first_x;
    int seg_id;

    // Set the start_idx of the segment as key
    int key() const { return first_x; }
    // [start_idx, end_idx]
    Segment(long double slope = 0, long double intercept = 0, K first_x = K(0), int seg_id = 0) : slope(slope), intercept(intercept), first_x(first_x) , seg_id(seg_id) {};
    
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
        double pos = slope *k + intercept;
        
        return (slope * k + intercept);
    }

};

class OptimalPLR
{

public:
    /**
     * @brief Construct a new OptimalPLR object
     *
     * @param error error bounds of this implementation
     */
    OptimalPLR(double error) : epsilon(error) {}

    std::vector<Segment> segmentData(const std::vector<Point> &points)
    {
        std::vector<Segment> segments;
        int seg_id = 0;
        if (points.empty())
            return segments;

        // Initialization: the first two points acting as the initial hull and segement
        size_t start = 0;
        Point sa = {points[0].x, points[0].y + epsilon};
        Point sc = {points[1].x, points[1].y - epsilon};
        Point sb = {points[0].x, points[0].y - epsilon};
        Point sd = {points[1].x, points[1].y + epsilon};

        double rho_max = computeSlope(sb, sd);
        double rho_min = computeSlope(sa, sc);

        std::deque<Point> upperHull = {sa, sd};
        std::deque<Point> lowerHull = {sb, sc};

        for (size_t i = 2; i < points.size(); ++i)
        {
            Point s = {points[i].x, points[i].y};

            if (!isWithinBounds(s, rho_max, sc, sd, rho_min))
            {
                // Check if this point is within the bounds , if not , we need to update the segement and initialize again
                Point so = findIntersection(sa, sc, sb, sd);
                double finalSlope = (rho_min + rho_max) / 2;
                segments.push_back({finalSlope, computeIntercept(so, finalSlope),points[start].x,seg_id++});

                start = i;
                sa = {points[start].x, points[start].y + epsilon};
                sc = {points[i + 1].x, points[i + 1].y - epsilon};
                sb = {points[start].x, points[start].y - epsilon};
                sd = {points[i + 1].x, points[i + 1].y + epsilon};

                rho_max = computeSlope(sb, sd);
                rho_min = computeSlope(sa, sc);
                upperHull.clear();
                lowerHull.clear();
                upperHull = {sa, sd};
                lowerHull = {sb, sc};
                i++;

                continue;
            }
            else
            {
                Point s_upper = {s.x, s.y + epsilon};
                Point s_lower = {s.x, s.y - epsilon};
                bool upper = 0;
                bool lower = 0;

                auto slope1 = computeSlope(s_lower,sa);
                auto slope2 = computeSlope(s_upper,sb);
                
                if ( slope1 > rho_min )
                {
                    // Find the max slope point and delete the points before it so that we can update the min slope
                    lower = 1;
                    sa = findMaxSlopePoint(upperHull, s_lower);
                    sc = s_lower;
                    deletePointsBefore(upperHull, sa);
                    rho_min = computeSlope(sa, sc);
                }

                if ( slope2 < rho_max )
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

        Point so = findIntersection(sa, sc, sb, sd);
        double finalSlope = (rho_min + rho_max) / 2;
        segments.push_back({finalSlope, computeIntercept(so, finalSlope), points[start].x, seg_id++});

        return segments;
    }

private:
    double epsilon;

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
        return (p2.y - p1.y) / (p2.x - p1.x);
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
        return p.y - slope * p.x;
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
        return (s.y + epsilon >= rho_min * (s.x - sc.x) + sc.y) && (s.y - epsilon <= rho_max * (s.x - sd.x) + sd.y);
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
        double slope1 = (p2.y - p1.y) / (p2.x - p1.x);
        double slope2 = (p3.y - p2.y) / (p3.x - p2.x);
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
        while (!hull.empty() && hull.front().x < reference.x)
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
    Point findIntersection(Point sa, Point sc, Point sb, Point sd)
    {
        double a1 = computeSlope(sa, sc);
        double b1 = computeIntercept(sa, a1);

        double a2 = computeSlope(sb, sd);
        double b2 = computeIntercept(sb, a2);

        if (a1 == a2)      // Parallel lines
            return {0, 0}; // or handle as needed
        double x_intersect = (b2 - b1) / (a1 - a2);
        double y_intersect = a1 * x_intersect + b1;

        return {x_intersect, y_intersect};
    }
};

/**
     * @brief Translate the segments in [first,last) into points
     * 
     * @tparam RandomIt the iterator of pointed segment 
     */
    auto translate(std::vector<Segment>segments,int first, int last) -> std::vector<Point>{

        
        std::vector<Point> internal_nodes;    
       
        int n = last - first;
        if(n == 0)
            return internal_nodes;
        
        internal_nodes.clear();
        internal_nodes.reserve(n);
        for(auto it = first; it < last;++it){
            // Get the first_idx and the seg_id of this segment
            Segment seg = segments[it];
            Point p(seg.first_x,seg.seg_id);
            internal_nodes.push_back(p);
        }

        return internal_nodes;
    }

size_t make_segmentation(int n, double epsilon, 
                            std::vector<Point> in, std::vector<Segment> &out)
{
    OptimalPLR opt(epsilon);
    auto segments = opt.segmentData(in);
    out.insert(out.end(), segments.begin(), segments.end());
    return segments.size();
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
size_t make_segmentation_par(int n, double epsilon, 
                            std::vector<Point> in, std::vector<Segment> &out,
                            int parallelism =16)
{
    
    // printf("The number of threads is %d\n", parallelism);

    int chunk_size = std::max(n / parallelism, 1000);  // Prevent too many threads for small `n`

    if (parallelism == 1 || n < 1ull << 15)
        return make_segmentation(n, epsilon, in, out);
    OptimalPLR opt(epsilon);
    std::vector<std::vector<Segment>> local_segments(parallelism);

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
void checkForEpsilon(std::vector<Point> data, std::vector<Segment> segments, int first = 0, int last = 5){
    
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
    while (data[i].x <= end)
    {
        if (data[i].x == last_x)
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
        printf("The real position is %f and the approximate position is %f\n",data[i].y,(slope*data[i].x + intercept));
        printf("The residual is %f for %dth element\n",std::abs(i - (slope*data[i].x + intercept)),i);
        auto residual = std::abs(i - (slope*data[i].x + intercept));
        if(residual>max_residual){
            max_residual = residual;
        }
        i++;
    }
}
}