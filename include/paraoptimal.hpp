
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

namespace PGM_P::internal{
template <typename X, typename Y>
class ParaoptimalPiecewiseLinearModel
{
private:
    using SX = LargeSigned<X>;
    using SY = LargeSigned<Y>;


    const Y epsilon;

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

    struct Line{
        // A line is defined by a point transformed in the feasible region
        double slope;  // slope of the line
        double intercept;  // intercept of the line
        Line(double slope, double intercept) : slope(slope), intercept(intercept) {
            if (slope == 0 && intercept == 0)
                throw std::invalid_argument("Line cannot be zero");
        }

        Line(const Point &p ,const bool lower, const Y epsilon){
            if(!lower){
                slope = static_cast<double> (-p.x);
                intercept = static_cast<double> (p.y + epsilon);
            }
            if(lower){
                slope = static_cast<double> (-p.x);
                intercept = static_cast<double> (p.y - epsilon);
            }
        }
    };



    struct FeasibleRegion
    {
        using RPoint = std::pair<double, double>;   // the first element is the slope, the second element is the intercept
        RPoint pl,pr;   // the point in feasible region that is the leftmost and rightmost point
        std::vector<RPoint> points_upper;   // the points in the feasible region that are placed upper on the polygon
        std::vector<RPoint> points_lower;   // the points in the feasible region that are placed lower on the polygon
        Y epsilon;

        // Default constructor
        FeasibleRegion() : pl{0, 0}, pr{0, 0},epsilon(0) {};

        FeasibleRegion(Point p1, Point p2, Y epsilon): epsilon(epsilon){
            // Initialize the feasible region with two points p1 and p2
            if(p2.x == p1.x){
                throw std::invalid_argument("can not handle two same points");
            }
            double denom = (static_cast<double>(p2.x) - static_cast<double>(p1.x));
            pl.first = (p2.y-p1.y-2*epsilon) / denom;
            pl.second = (p2.x*p1.y-p1.x*p2.y+epsilon*(p2.x+p1.x))/ denom;
            pr.first = (p2.y-p1.y+2*epsilon) / denom;
            pr.second = (p2.x*p1.y-p1.x*p2.y-epsilon*(p2.x+p1.x))/ denom;
           
            RPoint p1_upper, p1_lower;
            p1_upper.first = (p2.y-p1.y)/ denom;
            p1_upper.second = (p2.x*p1.y-p1.x*p2.y+epsilon*(p2.x-p1.x))/ denom;
       
            p1_lower.first = (p2.y-p1.y)/ denom;
            p1_lower.second = (p2.x*p1.y-p1.x*p2.y-epsilon*(p2.x-p1.x))/ denom;
            // Add the points to the feasible region
            points_upper.push_back(p1_upper);
            points_lower.push_back(p1_lower);

        }

        

        // Check the position of a point relative to a certain line
        bool position(const RPoint &p, Line &l)
        {
            // Check if the point is above the line
            double y = l.slope * p.first + l.intercept;
            return p.second > y; // true if the point is above the line, false otherwise
        }

        // Calculate the intersection of two points and a line
        RPoint intersection(const RPoint &p1, const RPoint &p2, Line &l)
        {
            // Calculate the intersection of the line with the segment p1-p2
            double a = (p2.second - p1.second) / (p2.first - p1.first); // slope of the segment
            double b = p1.second - a * p1.first; // intercept of the segment

            double x = (l.intercept - b) / (l.slope - a);
            double y = l.slope * x + l.intercept;
            return {x, y};
        }

        /**
         * @brief Calculate the intersection of a point and a convex hull and update it
         * 
         * @param p the point newly added to the feasible region 
         * @param p_lower the lower line or the upper line defined by this point 
         * @param h_lower the lower convex hull or the upper one 
         * @return RPoint the intersection point of the point and the convex hull 
         */
        RPoint intersect(const Point &p,bool p_lower = 0, bool h_lower = 0)
        {
            RPoint result;
            Line l = Line(p, p_lower, epsilon);     // Define the line that is used to intersect with the convex hull

            bool remove_before = false;
            bool remove_after = false; 
            if(p_lower != h_lower)
                // We need to calculate the leftmost or the rightmost point of the convex hull
                // And remove the point in the convex hull before the intersection point
                remove_before = true;
            if(p_lower == h_lower)
                // We get the new end of the convex hull 
                // remove the point in the convex hull after the intersection point
                remove_after = true; 

            if(h_lower){
                // We are calculating the intersection with the lower convex hull

                RPoint pre = pr;    // Set the first point as the rightmost point
                bool pos_pre = position(pre, l);
                bool pos_cur = pos_pre;    // Initialize the current position as the previous position
                for(size_t i = 0; i < points_lower.size(); ++i){
                    // tranverse the lower points (counter-clockwise)
                    RPoint cur = points_lower[i];
                    pos_cur = position(cur, l);
                    
                    if(pos_cur == pos_pre){
                        // If the current point is on the same side as the previous point,
                        // which means the line does not intersect with the feasible region
                        pos_pre = pos_cur;
                        pre = cur;
                        continue;
                    }
                    else{
                        // If the current point is on the opposite side as the previous point,
                        // which means the line intersects with the feasible region
                        result = intersection(pre, cur, l);
                        
                        // Update the convex hull
                        if(remove_before){
                            // Remove the points before the intersection point and update the pr
                            points_lower.erase(points_lower.begin() , points_lower.begin() + i);
                            pr = result;
                            break;
                        }
                        else if(remove_after){
                            // Remove the point and the points after the intersection point and push the new point
                            points_lower.erase(points_lower.begin() + i  , points_lower.end());
                            points_lower.push_back(result);
                            break;
                        }
                    }
                
                }

                // If we reach here, it means the line will intersect with the last point of the convex hull and pl
                RPoint end = pl;
                bool pos_end = position(end, l);
                if(pos_end != pos_cur){
                    // If the last point is on the opposite side of the line, we can calculate the intersection
                    result = intersection(pre, end, l);
                    if(remove_before){
                        // Remove the points before the intersection point
                        points_lower.erase(points_lower.begin() , points_lower.end());
                        pr = result; // Update the rightmost point
                    }
                    else if(remove_after){
                        // Remove the point and the points after the intersection point
                        points_lower.push_back(result);
                    }
                }
            } 
 
            if(!h_lower){

                // We are calculating the intersection with the upper convex hull

                RPoint pre = pl;    // Set the first point as the leftmost point
                bool pos_pre = position(pre, l);
                bool pos_cur = pos_pre;    // Initialize the current position as the previous position
                for(size_t i = 0; i < points_upper.size(); ++i){
                    // tranverse the lower points (counter-clockwise)
                    RPoint cur = points_upper[i];
                    pos_cur = position(cur, l);
                    
                    if(pos_cur == pos_pre){
                        // If the current point is on the same side as the previous point,
                        // which means the line does not intersect with the feasible region
                        pos_pre = pos_cur;
                        pre = cur;
                        continue;
                    }
                    else{
                        // If the current point is on the opposite side as the previous point,
                        // which means the line intersects with the feasible region
                        result = intersection(pre, cur, l);
                        
                        // Update the convex hull
                        if(remove_before){
                            // Remove the points before the intersection point and update the pl
                            points_upper.erase(points_upper.begin() , points_upper.begin() + i);
                            pl = result;
                        }
                        else if(remove_after){
                            // Remove the point and the points after the intersection point and push the new point
                            points_upper.erase(points_upper.begin() + i  , points_upper.end());
                            points_upper.push_back(result);
                        }
                    }
                
                }

                // If we reach here, it means the line will intersect with the last point of the convex hull and pl
                RPoint end = pr;
                bool pos_end = position(end, l);
                if(pos_end != pos_pre){
                    // If the last point is on the opposite side of the line, we can calculate the intersection
                    result = intersection(pre, end, l);
                    if(remove_before){
                        // Remove the points before the intersection point and update the pl
                        points_upper.erase(points_upper.begin() , points_upper.end());
                        pl = result; // Update the leftmost point
                    }
                    else if(remove_after){
                        // Remove the point and the points after the intersection point
                        points_upper.push_back(result);
                    }
                }
            } 
            
            return result;
        }

        // Add a point to the feasible region
        bool add_point(const Point &p)
        {
        
            // Check if the point is within the feasible region
            bool outside1 = pl.first * p.x + pl.second > p.y + epsilon;
            bool outside2 = pr.first * p.x + pr.second < p.y - epsilon;
            if (outside1 || outside2){
                // The point is outside the feasible region, so we need to update the feasible region
                return false;
            }
            
            // Add the point by clipping the feasible region
            // Update the pr after the upper line intersecting with the p_lower and update the p_lower
            // And update the p_upper after the upper line intersecting with the p_upper
            if(pr.first * p.x + pr.second > p.y + epsilon){
                bool p_lower = 0;
                bool h_lower = 1;
                intersect(p ,p_lower , h_lower = 1); // Update the pr 
                intersect(p ,p_lower , h_lower = 0); // Update the points_upper


            }
            // Update the pl after the lower line intersecting with the p_upper and update the p_upper
            // And update the p_lower after the lower line intersecting with the p_lower
            if(pl.first * p.x + pl.second < p.y - epsilon){
                bool p_lower = 1;
                bool h_lower = 1;
                intersect(p ,p_lower , h_lower = 0); // Update the pl
                intersect(p ,p_lower , h_lower = 1); // Update the points_lower   
            }

            return true;
        }
    
    
        // Get the feasible point in the region
        RPoint get_feasible_point() const
        {
            RPoint p = {(pl.first + pr.first) / 2, (pl.second + pr.second) / 2};
            return p;
        }
    };

    X first_x = 0;
    Y first_y = 0;
    X last_x = 0;
    int points = 0;
    FeasibleRegion region;

public:
    class CanonicalSegment;

    explicit ParaoptimalPiecewiseLinearModel(Y epsilon) : epsilon(epsilon)
    {
        if (epsilon < 0) throw std::invalid_argument("epsilon cannot be negative"); }
    Y get_epsilon(){
        return epsilon;
    }

    bool add_point(const X &x, const Y &y)
    {
        
        last_x = x;
        Point p{x,y};

        if (points == 0)
        {
            // It is the first point of the segment
            first_x = x;
            first_y = y;
            points++;
            region = FeasibleRegion();
            return true;
        }

        if (points == 1)
        {
            points++;
            Point p2{first_x,first_y};
            region = FeasibleRegion(p2, p, epsilon);
            return true;
        }

        if(region.add_point(p)){
            // If the point is added to the feasible region, we can update the first_x and first_y
            points++;
            return true;
        }else{
            points = 0;
            
            return false;
        }
    }

    CanonicalSegment get_segment()
    {
        return CanonicalSegment(region ,first_x, last_x);
    }

    void reset()
    {
        points = 0;
    }
};

template <typename X, typename Y>
class ParaoptimalPiecewiseLinearModel<X, Y>::CanonicalSegment
{
public:
    CanonicalSegment() = default;

    friend class ParaoptimalPiecewiseLinearModel;

    FeasibleRegion region;
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

    CanonicalSegment(FeasibleRegion region, X first_x, X last_x) : region(region), first_x(first_x), last_x(last_x) {};

    bool one_point() const
    {
        return region.pl.first == 0 && region.pl.second == 0 &&
               region.pr.first == 0 && region.pr.second == 0;
    }

    Segment get_Canonicalsegment() const
    {
        
        if (one_point())
            return {0, 0 , first_x, last_x};

        // Calculate the slope and intercept of the segment
        auto p = region.get_feasible_point();
        long double slope = p.first;
        long double intercept = p.second;

        Segment result = Segment(slope, intercept, first_x, last_x);
        return result;
    }

		
};

template <typename Fin, typename Fout>
size_t make_segmentation(size_t n, size_t start, size_t end, double epsilon, Fin in, Fout out)
{
    using K = typename std::invoke_result_t<Fin, size_t>;
    size_t c = 0;
    ParaoptimalPiecewiseLinearModel<K, int> opt(epsilon);
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
    using canonical_segment = typename ParaoptimalPiecewiseLinearModel<K, int>::CanonicalSegment;
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
template <typename Fin, typename SegmentType>
void checkForEpsilon(Fin data, std::vector<SegmentType>segments ,int begin = 0,int ed = 5, double epsilon = 0.1){
    
    auto start = segments[begin].first_x;
    auto end = segments[ed].last_x;
    
    int segmnt_idx = begin;
    auto seg = segments[segmnt_idx];
    long double slope = seg.slope;
    long double intercept = seg.intercept;
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
        printf("The predicted value is %Lf and the actual value is %d\n", slope*data(i) + intercept, i);
        i++;
    }
    printf("\n");
    }
}
