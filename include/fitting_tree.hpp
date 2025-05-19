#pragma once

#include <vector>
#include <memory>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include "greedy.hpp"

namespace fitting {

template<typename K, size_t epsilon = 2>
class FITTING_Tree {
public:
    // 存储Greedy算法生成的线段
    struct Segment {
        K first_key;  // 线段的起始键值
        double slope;
        double intercept;
        size_t start_pos;  // 线段对应的数据起始位置
        size_t end_pos;    // 线段对应的数据结束位置

        Segment(K key, double s, double i, size_t start, size_t end)
            : first_key(key), slope(s), intercept(i), start_pos(start), end_pos(end) {}
    };

    FITTING_Tree(size_t degree = 4) : max_degree(degree), data_size(0) {
        root = std::make_shared<BPlusNode>(true);
        segments.reserve(1000000);  // 预分配空间
    }

    // 构建树
    void build(const std::vector<K>& input_data);

    // 搜索接口
    size_t search(K key) const;

    // 获取树的信息
    size_t height() const;
    size_t size_in_bytes() const;
    std::vector<size_t> get_levels_offsets() const;
    
    // 获取线段信息
    size_t get_segment_count() const { return segments.size(); }
    size_t get_segment_size_bytes() const { return segments.size() * sizeof(Segment); }

private:
    // B+树节点结构
    struct BPlusNode {
        std::vector<K> keys;
        std::vector<std::shared_ptr<BPlusNode>> children;
        bool is_leaf;
        std::shared_ptr<BPlusNode> next;  // 用于叶子节点之间的连接

        BPlusNode(bool leaf = false) : is_leaf(leaf) {
            keys.reserve(4);  // 预分配空间
            children.reserve(4);
        }
    };

    std::shared_ptr<BPlusNode> root;
    size_t max_degree;  // B+树的最大度数
    std::vector<Segment> segments;  // 存储所有线段
    size_t data_size;  // 只存储数据大小

    // 辅助函数
    void build_bplus_tree();
    std::shared_ptr<BPlusNode> find_leaf(K key) const;
};

// 实现部分
template<typename K, size_t epsilon>
void FITTING_Tree<K, epsilon>::build(const std::vector<K>& input_data) {
    data_size = input_data.size();
    segments.clear();  // 清空之前的线段
    
    // 使用Greedy算法构建线段
    auto in_fun = [&input_data](auto i) { return input_data[i]; };
    auto out_fun = [this](auto cs) {
        auto segment = cs.get_Canonicalsegment();
        segments.emplace_back(segment.first_x, segment.slope, segment.intercept,
                            static_cast<size_t>(segment.intercept),
                            static_cast<size_t>(segment.intercept + epsilon));
    };
    
    // 一次性处理所有数据
    Greedy::internal::make_segmentation_par(input_data.size(), epsilon, in_fun, out_fun);

    // 构建B+树索引
    build_bplus_tree();
}

template<typename K, size_t epsilon>
void FITTING_Tree<K, epsilon>::build_bplus_tree() {
    if (segments.empty()) return;

    // 创建叶子节点层
    std::vector<std::shared_ptr<BPlusNode>> current_level;
    current_level.reserve(segments.size());  // 预分配空间
    
    for (const auto& segment : segments) {
        auto node = std::make_shared<BPlusNode>(true);
        node->keys.push_back(segment.first_key);
        current_level.push_back(node);
    }

    // 自底向上构建B+树
    while (current_level.size() > 1) {
        std::vector<std::shared_ptr<BPlusNode>> next_level;
        next_level.reserve((current_level.size() + max_degree - 2) / (max_degree - 1));  // 预分配空间
        
        for (size_t i = 0; i < current_level.size(); i += max_degree - 1) {
            auto parent = std::make_shared<BPlusNode>(false);
            size_t end = std::min(i + max_degree - 1, current_level.size());
            
            for (size_t j = i; j < end; ++j) {
                parent->children.push_back(current_level[j]);
                if (j < end - 1) {
                    parent->keys.push_back(current_level[j + 1]->keys.front());
                }
            }
            next_level.push_back(parent);
        }
        current_level = std::move(next_level);
    }

    root = current_level.front();
}

template<typename K, size_t epsilon>
size_t FITTING_Tree<K, epsilon>::search(K key) const {
    // 首先在B+树中定位到叶子节点
    auto leaf = find_leaf(key);
    if (!leaf) return -1;

    // 找到对应的线段
    auto it = std::upper_bound(segments.begin(), segments.end(), key,
        [](K k, const Segment& seg) { return k < seg.first_key; });
    
    if (it == segments.begin()) return 0;
    --it;

    // 使用线段预测位置
    double pos = it->slope * key + it->intercept;
    size_t predicted_pos = static_cast<size_t>(pos);
    
    // 线性搜索找到精确位置
    predicted_pos = std::max(predicted_pos, it->start_pos);
    predicted_pos = std::min(predicted_pos, it->end_pos);
    
    return predicted_pos;
}

template<typename K, size_t epsilon>
std::shared_ptr<typename FITTING_Tree<K, epsilon>::BPlusNode>
FITTING_Tree<K, epsilon>::find_leaf(K key) const {
    auto current = root;
    while (!current->is_leaf) {
        auto it = std::upper_bound(current->keys.begin(), current->keys.end(), key);
        size_t index = std::distance(current->keys.begin(), it);
        current = current->children[index];
    }
    return current;
}

template<typename K, size_t epsilon>
size_t FITTING_Tree<K, epsilon>::height() const {
    size_t h = 1;
    auto current = root;
    while (!current->is_leaf) {
        h++;
        current = current->children[0];
    }
    return h;
}

template<typename K, size_t epsilon>
size_t FITTING_Tree<K, epsilon>::size_in_bytes() const {
    size_t total_size = 0;
    
    // 计算B+树节点大小
    std::function<void(const std::shared_ptr<BPlusNode>&)> calculate_node_size = 
        [&](const std::shared_ptr<BPlusNode>& node) {
            total_size += sizeof(BPlusNode);
            total_size += node->keys.size() * sizeof(K);
            total_size += node->children.size() * sizeof(std::shared_ptr<BPlusNode>);
            
            if (!node->is_leaf) {
                for (const auto& child : node->children) {
                    calculate_node_size(child);
                }
            }
        };
    
    calculate_node_size(root);
    
    // 计算线段大小
    total_size += segments.size() * sizeof(Segment);
    
    return total_size;
}

template<typename K, size_t epsilon>
std::vector<size_t> FITTING_Tree<K, epsilon>::get_levels_offsets() const {
    std::vector<size_t> offsets;
    size_t current_offset = 0;
    
    // 计算每一层的节点数量
    std::vector<std::shared_ptr<BPlusNode>> current_level = {root};
    while (!current_level.empty()) {
        offsets.push_back(current_offset);
        current_offset += current_level.size();
        
        std::vector<std::shared_ptr<BPlusNode>> next_level;
        for (const auto& node : current_level) {
            if (!node->is_leaf) {
                next_level.insert(next_level.end(), 
                                node->children.begin(), 
                                node->children.end());
            }
        }
        current_level = std::move(next_level);
    }
    
    return offsets;
}

} // namespace fitting 