
#pragma once

#include "swing.hpp"
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
#include <chrono>


namespace Swing {

#define PGM_SUB_EPS(x, epsilon) ((x) <= (epsilon) ? 0 : ((x) - (epsilon)))
#define PGM_ADD_EPS(x, epsilon, size) ((x) + (epsilon) + 2 >= (size) ? (size) : (x) + (epsilon) + 2)

/**
 * A struct that stores the result of a query to a @ref PGMIndex, that is, a range [@ref lo, @ref hi)
 * centered around an approximate position @ref pos of the sought key.
 */
struct ApproxPos {
    size_t pos; ///< The approximate position of the key.
    size_t lo;  ///< The lower bound of the range.
    size_t hi;  ///< The upper bound of the range.
    std::vector<long> residuals;
};

/**
 * @tparam K the type of the indexed keys
 * @tparam Parallelism the number of threads to use for parallelization (default is 16)
 * @tparam Floating the floating-point type to use for slopes
 */
template<typename K, size_t Parallelism = 16, typename Floating = float>
class PGMIndex {
protected:
        size_t Epsilon;
    size_t EpsilonRecursive;
    struct Segment;

    size_t n;                           ///< The number of elements this index was built on.
    K first_key;                        ///< The smallest element.
    std::vector<Segment> segments;      ///< The segments composing the index.
    std::vector<size_t> levels_offsets; ///< The starting position of each level in segments[], in reverse order.

    /// Sentinel value to avoid bounds checking.
    static constexpr K sentinel = std::numeric_limits<K>::has_infinity ? std::numeric_limits<K>::infinity()
                                                                       : std::numeric_limits<K>::max();

    template<typename RandomIt>
    static void build(RandomIt first, RandomIt last,
                      size_t epsilon, size_t epsilon_recursive,
                      std::vector<Segment> &segments,
                      std::vector<size_t> &levels_offsets) {
        auto n = (size_t) std::distance(first, last);
        if (n == 0)
            return;

        levels_offsets.push_back(0);
        segments.reserve(n / (epsilon * epsilon));

        if (*std::prev(last) == sentinel)
            throw std::invalid_argument("The value " + std::to_string(sentinel) + " is reserved as a sentinel.");

        auto build_level = [&](auto epsilon, auto in_fun, auto out_fun, size_t last_n) {
            auto n_segments = internal::make_segmentation_par(last_n, epsilon, in_fun, out_fun,Parallelism);
            if (segments.back() == sentinel)
                --n_segments;
            else {
                if (segments.back()(sentinel - 1) < last_n)
                    segments.emplace_back(*std::prev(last) + 1, 0, last_n); // Ensure keys > last are mapped to last_n
                segments.emplace_back(sentinel, 0, last_n);
            }
            return n_segments;
        };

        // Build first level
        auto in_fun = [&](auto i) { return K(first[i]); };
        auto out_fun = [&](auto cs) { segments.emplace_back(cs); };
        auto last_n = build_level(epsilon, in_fun, out_fun, n);
        levels_offsets.push_back(segments.size());

        // Build upper levels
        while (epsilon_recursive && last_n > 1) {
            auto offset = levels_offsets[levels_offsets.size() - 2];
            auto in_fun_rec = [&](auto i) { return segments[offset + i].key; };
            last_n = build_level(epsilon_recursive, in_fun_rec, out_fun, last_n);
            // Check the correctness of the segments
            
            levels_offsets.push_back(segments.size());
        }
    }

    /**
     * Returns the segment responsible for a given key, that is, the rightmost segment having key <= the sought key.
     * Also captures the residual (prediction error) at each recursive level.
     * @param key the value of the element to search for
     * @return a pair of: 
     *   - an iterator to the segment responsible for the given key
     *   - a vector of residuals (one per recursive level)
     */
    auto segment_for_key(const K &key) const -> std::pair<typename std::vector<typename PGMIndex<K, Parallelism, Floating>::Segment>::const_iterator, std::vector<long>>{
        if (EpsilonRecursive == 0) {
            auto it = std::prev(std::upper_bound(segments.begin(), segments.begin() + segments_count(), key));
            return {it, {}}; // No levels, no residuals
        }

        std::vector<long> residuals_per_level;

        auto it = segments.begin() + *(levels_offsets.end() - 2);
        for (int l = int(height()) - 2; l >= 0; --l) {
            auto level_begin = segments.begin() + levels_offsets[l];
            auto level_size = levels_offsets[l + 1] - levels_offsets[l] - 1;

            // Predict position
            auto pos1 = (*it)(key);
            auto pos = std::min<size_t>((pos1 > 0 ? pos1 : 0),std::next(it)->intercept);

            auto lo = level_begin + PGM_SUB_EPS(pos, EpsilonRecursive + 2);

            static constexpr size_t linear_search_threshold = 8 * 64 / sizeof(Segment);

            if (EpsilonRecursive <= linear_search_threshold) {
                if (lo != level_begin) {
                    for (; std::prev(lo)->key > key; --lo);
                    lo = std::prev(lo);
                }
                for (; std::next(lo)->key <= key; ++lo);
                it = lo;
            } else {
                for (; lo != level_begin + EpsilonRecursive / 2 && lo->key > key; lo -= EpsilonRecursive / 2);
                auto hi = level_begin + PGM_ADD_EPS(pos, EpsilonRecursive, level_size);
                it = std::prev(std::upper_bound(lo, hi, key));
            }

            // Compute residual for this level
            auto actual_pos = std::distance(level_begin, it);
            auto residual = std::abs(static_cast<long>(pos) - actual_pos);
            residuals_per_level.push_back(residual);

            // printf("Level %d | Key %ld | Predicted = %ld | Actual = %ld | Residual = %ld\n",
            //     l, key, pos, actual_pos, residual);
        }

        return {it, residuals_per_level};
    }

public:

    size_t epsilon_value = Epsilon;

    /**
     * Constructs an empty index.
     */
    PGMIndex() = default;

    /**
     * Constructs the index on the given sorted vector.
     * @param data the vector of keys to be indexed, must be sorted
     */
    explicit PGMIndex(const std::vector<K> &data) : PGMIndex(data.begin(), data.end()) {}

    /**
     * Constructs the index on the sorted keys in the range [first, last).
     * @param first, last the range containing the sorted keys to be indexed
     */
    template<typename RandomIt>
    
    PGMIndex(RandomIt first, RandomIt last, size_t Epsilon, size_t EpsilonRecursive)
        : Epsilon(Epsilon),
            EpsilonRecursive(EpsilonRecursive),
         n(std::distance(first, last)),
          first_key(n ? *first : K(0)),
          segments(),
          levels_offsets() {
        build(first, last, Epsilon, EpsilonRecursive, segments, levels_offsets);
    }

    /**
     * Returns the approximate position and the range where @p key can be found.
     * @param key the value of the element to search for
     * @return a struct with the approximate position and bounds of the range
     */
    ApproxPos search(const K &key) const {
        auto k = std::max(first_key, key);
        auto a = segment_for_key(k);
        auto it = a.first;
        auto res = a.second;
        auto pos = std::min<size_t>((*it)(k), std::next(it)->intercept);
        auto lo = PGM_SUB_EPS(pos, Epsilon);
        auto hi = PGM_ADD_EPS(pos, Epsilon, n);
        return {pos, lo, hi,res};
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
    size_t size_in_bytes() const { return segments.size() * sizeof(Segment) + levels_offsets.size() * sizeof(size_t); }

    /**
     * Returns the offsets of the segments in each level of the index.
     * @return a vector with the offsets of the segments in each level
     */
    std::vector<int> get_levels_offsets() const {
        std::vector<int> levels;
        for (size_t i = 1; i < levels_offsets.size(); ++i) {
            levels.push_back(levels_offsets[i] - levels_offsets[i-1]);
        }
        return levels;
    }
};

#pragma pack(push, 1)


template<typename K, size_t Parallelism,typename Floating>
struct PGMIndex<K, Parallelism, Floating>::Segment {
    K key;              ///< The first key that the segment indexes.
    float slope;     ///< The slope of the segment.
    uint32_t intercept; ///< The intercept of the segment.

    Segment() = default;

    Segment(K key, Floating slope, int32_t intercept) : key(key), slope(slope), intercept(intercept) {};

    explicit Segment(const typename internal::SwingPiecewiseLinearModel<K, int>::CanonicalSegment &cs)
    : key(cs.first_x)
     {
        auto segment = cs.get_Canonicalsegment(key);
        auto cs_slope = segment.slope;
        auto cs_intercept = segment.intercept;
        if (cs_intercept > std::numeric_limits<decltype(intercept)>::max())
            throw std::overflow_error("Change the type of Segment::intercept to uint64");
        // if (cs_intercept < 0)
        //     throw std::overflow_error("Unexpected intercept < 0");
        
        slope = cs_slope;
        intercept = cs_intercept;
    }

    friend inline bool operator<(const Segment &s, const K &k) { return s.key < k; }
    friend inline bool operator<(const K &k, const Segment &s) { return k < s.key; }
    friend inline bool operator<(const Segment &s, const Segment &t) { return s.key < t.key; }

    operator K() { return key; };

    /**
     * Returns the approximate position of the specified key.
     * @param k the key whose position must be approximated
     * @return the approximate position of the specified key
     */
    inline double operator()(const K &k) const {
        size_t pos;
        if constexpr (std::is_same_v<K, int64_t> || std::is_same_v<K, int32_t>)
            pos = size_t(slope * double(std::make_unsigned_t<K>(k) - key));
        else
            pos = size_t(slope * double(k - key));
        return pos + intercept;
    }
};

#pragma pack(pop)

}