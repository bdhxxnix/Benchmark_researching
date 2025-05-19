#include "optimalPLR_Foundation.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include <cmath>

namespace PGM_C{
    /**
     * A struct that stores the result of a query to a @ref PGMIndex, that is, a range [@ref lo, @ref hi)
     * centered around an approximate position @ref pos of the sought key.
     */
struct ApproxPos {
size_t pos; ///< The approximate position of the key.
size_t lo;  ///< The lower bound of the range.
size_t hi;  ///< The upper bound of the range.
};
/**
 * A space-efficient index that enables fast search operations on a sorted sequence of numbers.
 *
 * A search returns a struct @ref ApproxPos containing an approximate position of the sought key in the sequence and
 * the bounds of a range where the sought key is guaranteed to be found if present.
 * If the key is not present, a @ref std::lower_bound search on the range finds a key that is greater or equal to the
 * sought key, if any.
 * In the case of repeated keys, the index finds the position of the first occurrence of a key.
 *
 * The @p Epsilon template parameter should be set according to the desired space-time trade-off. A smaller value
 * makes the estimation more precise and the range smaller but at the cost of increased space usage.
 *
 * Internally the index uses a succinct piecewise linear mapping from keys to their position in the sorted order.
 * This mapping is represented as a sequence of linear models (segments) which, if @p EpsilonRecursive is not zero, are
 * themselves recursively indexed by other piecewise linear mappings.
 *
 * @tparam K the type of the indexed keys
 * @tparam Epsilon controls the size of the returned search range
 * @tparam EpsilonRecursive controls the size of the search range in the internal structure
 * @tparam Floating the floating-point type to use for slopes
 */
template<typename K, size_t Epsilon = 64, size_t EpsilonRecursive = 4, typename Floating = float>
class PGMIndex {
protected:


    static_assert(Epsilon > 0);
    using Segment = PGM_C::internal::Segment;
    using Point = PGM_C::internal::Point;


    size_t n;                           ///< The number of elements this index was built on.
    K first_key;                        ///< The smallest element.
    std::vector<Segment> segments;      ///< The segments composing the index.
    std::vector<size_t> levels_offsets; ///< The starting position of each level in segments[], in reverse order.


    /// Sentinel value to avoid bounds checking.
    static constexpr K sentinel = std::numeric_limits<K>::has_infinity ? std::numeric_limits<K>::infinity()
                                                                        : std::numeric_limits<K>::max();
    

    
    /**
     * @brief build the whole PGM-Index by recursively using PLR algorithm  
     * 
     * @param first the first iterator of the input data 
     * @param last the last iterator of the input data
     * @param epsilon the episilon value for the first level
     * @param epsilon_recursive the episilon value for the recursive levels
     * @param segments the segments of the whole PGM-Index
     * @param levels_offsets the offsets of the segments in each level in reverse order
     */
    static void build(const std::vector<Point> &data,
                        size_t epsilon, size_t epsilon_recursive,
                        std::vector<Segment> &segments,
                        std::vector<size_t> &levels_offsets) {
        auto n = data.size();
        if (n == 0)
            return;

        levels_offsets.push_back(0);
        segments.reserve(n / (epsilon * epsilon));

        if (data.end()->x == sentinel)
            throw std::invalid_argument("The value " + std::to_string(sentinel) + " is reserved as a sentinel.");

        auto build_level = [&](auto epsilon, std::vector<Point> input, std::vector<Segment> &output, size_t last_n) {
            auto n_segments = PGM_C::internal::make_segmentation_par(last_n, epsilon,input,output);
            if (segments.back() == sentinel)
                --n_segments;
            return n_segments;
        };

        // Build first level
        auto last_n = build_level(epsilon,data,segments, n);
        levels_offsets.push_back(segments.size());
        // PGM_C::internal::checkForEpsilon(data,segments,0,10); // Check the generated segments
        //print out the current segments
        
        // Build upper levels
        while (epsilon_recursive && last_n > 1) {
            auto offset = levels_offsets[levels_offsets.size() - 2];
            auto last = offset + last_n;  
            std::vector<Point> input = PGM_C::internal::translate(segments,offset,last);
            // Check the input data
            
            last_n = build_level(epsilon_recursive, input, segments, last_n);
            // Check the generated segments
            // auto offset1 = levels_offsets[levels_offsets.size() - 1];
            // PGM_C::internal::checkForEpsilon(input,segments,offset1,offset1 + std::min<size_t>(last_n,size_t(10))); // Check the generated segments);
            // printf("\n"); 
            levels_offsets.push_back(segments.size());
            //print out the current segments
        }
    }

    /**
     * Returns the segment responsible for a given key, that is, the rightmost segment having key <= the sought key.
     * @param key the value of the element to search for
     * @return an iterator to the segment responsible for the given key
     */
    auto segment_for_key(const K &key) const {
        if constexpr (EpsilonRecursive == 0) {
            return std::prev(std::upper_bound(segments.begin(), segments.begin() + segments_count(), key));
        }

        auto it = segments.begin() + *(levels_offsets.end() - 2);
        for (auto l = int(height()) - 2; l >= 0; --l) {
            auto level_begin = segments.begin() + levels_offsets[l];
            // auto po = (*it)(key) ;
            // auto po2 = std::next(it)->intercept;
            auto pos = std::min<double>((*it)(key) < 0 ? (*it)(key):0 , std::next(it)->intercept);
            // printf("The position of the key is on level %d , pos :%lf\n",l,pos);
            auto lo = level_begin + PGM_SUB_EPS(pos, EpsilonRecursive + 1);

            static constexpr size_t linear_search_threshold = 8 * 64 / sizeof(Segment);
            if constexpr (EpsilonRecursive <= linear_search_threshold) {
                for (; std::next(lo)->key() <= key; ++lo)
                    continue;
                it = lo;
            } else {
                auto level_size = levels_offsets[l + 1] - levels_offsets[l] - 1;
                auto hi = level_begin + PGM_ADD_EPS(pos, EpsilonRecursive, level_size);
                it = std::prev(std::upper_bound(lo, hi, key));
            }
        }
        return it;
    }

    public:

    static constexpr size_t epsilon_value = Epsilon;

    /**
     * Constructs an empty index.
     */
    PGMIndex() = default;


    /**
     * Constructs the index on the sorted keys in the range [first, last).
     * @param first, last the range containing the sorted keys to be indexed
     */
    PGMIndex(const std::vector<Point> &data)
        : n(data.size()),
            segments(),
            levels_offsets(),
            first_key(data.empty() ? 0 : data.front().x) {
        build(data, Epsilon, EpsilonRecursive, segments, levels_offsets);
    }

    /**
     * Returns the approximate position and the range where @p key can be found.
     * @param key the value of the element to search for
     * @return a struct with the approximate position and bounds of the range
     */
    ApproxPos search(const K &key) const {
        auto k = std::max(first_key, key);
        auto it = segment_for_key(k);
        
        auto ori_pos = std::min<double>((*it)(k) < 0 ? 0 : (*it)(k), std::next(it)->intercept);
        size_t pos = std::max<size_t>(0, std::min<size_t>(n - 1, static_cast<size_t>(ori_pos)));
        auto lo = PGM_SUB_EPS(pos, Epsilon);
        auto hi = PGM_ADD_EPS(pos, Epsilon, n);
        return {pos, lo, hi};
    }

    /**
     * Returns the number of segments in the last level of the index.
     * @return the number of segments
     */
    size_t segments_count() const { return segments.empty() ? 0 : levels_offsets[1] - 1; }

    /**
     * Returns the number of levels of the index.
     * @return the number of levels of the index
     */
    size_t height() const { return levels_offsets.size() - 1; }

    /**
     * Returns the size of the index in bytes.
     * @return the size of the index in bytes
     */
    size_t size_in_bytes() const { return segments.size() * sizeof(Segment) + levels_offsets.size() * sizeof(size_t) ; }

    /**
     * Returns the offsets of the segments in each level of the index.
     * @return a vector with the offsets of the segments in each level
     */
    std::vector<int> get_levels_offsets() const {
        std::vector<int> levels;
        for (size_t i = 0; i < levels_offsets.size(); ++i) {
            levels.push_back(levels_offsets[i]);
        }
        return levels;
    }
    };

    #pragma pack(push, 1)



    #pragma pack(pop)

}
