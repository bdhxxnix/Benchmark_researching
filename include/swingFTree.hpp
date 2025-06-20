#pragma once

#include <vector>
#include <utility>
#include <cstdint>
#include <memory>
#include <stx/btree_map.h>
#include "swing.hpp"

namespace Swing{

template<typename KeyType>
struct LeafNode {
    using Segment = typename Swing::internal::SwingPiecewiseLinearModel<KeyType, int>::CanonicalSegment::Segment; 
    std::vector<Segment> segments;

    size_t bytes_used() const {
        return sizeof(*this) + segments.size() * sizeof(Segment);
    }
};

template<typename KeyType>
class FittingTree {
public:

    FittingTree(double error_bound = 16.0):error_bound_(error_bound) {
        if (error_bound <= 0) {
            throw std::invalid_argument("Error bound must be positive");
        }
    };

    ~FittingTree() = default;

    void build(std::vector<KeyType> data) {

        using Segment = typename Swing::internal::SwingPiecewiseLinearModel<KeyType, int>::CanonicalSegment::Segment; 
        
        // Build the tree index using greedy segmentation
        auto in = [&](int i){return data[i];};
        std::vector<Segment> segments;
        auto out = [&](const auto& cs) { 
            auto first_x = cs.first_x;
            auto result = cs.get_Canonicalsegment(first_x);
            segments.push_back(result); };
        auto n = data.size();
        int parallelism = omp_get_num_procs();
        Swing::internal::make_segmentation_par(n, error_bound_, in, out, parallelism);

        for (const auto& seg : segments) {
            auto leaf = std::make_unique<LeafNode<KeyType>>();
            leaf->segments.push_back(seg);
            tree_index_.insert({seg.first_x, std::move(leaf)});
        }
    };
    
    const int search(KeyType key) const{
        auto it = tree_index_.upper_bound(key);
        if (it == tree_index_.begin()) return 0;
        if (it != tree_index_.begin()) --it;
        using Segment = typename LeafNode<KeyType>::Segment;
        const Segment& seg = it->second->segments[0];
        if (key >= seg.first_x && key <= seg.last_x) {
            size_t idx = seg.slope*(key- seg.first_x)+seg.intercept;
            return idx;
        } 

        return 0; // No segment found for the key
    };

    void clear() {
        tree_index_.clear();
    };

    size_t size_in_bytes() const {
        size_t total_size = sizeof(*this);
        for (const auto& [key, ptr] : tree_index_) {
            total_size += sizeof(key);
            if (ptr) total_size += ptr->bytes_used();
        }
        return total_size;
    } 
    
    

    size_t get_segment_count() const {
        size_t count = tree_index_.size();
        return count;
    };

    size_t get_segment_size_bytes() const {
        size_t total_size = 0;
        total_size += sizeof(LeafNode<KeyType>) * tree_index_.size();
        return total_size;
    };


private:
    using TreeMap = stx::btree_map<KeyType, std::shared_ptr<LeafNode<KeyType>>>;

    TreeMap tree_index_;
    double error_bound_;

};
}