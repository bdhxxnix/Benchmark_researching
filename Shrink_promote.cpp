#include <vector>
#include <iostream>

struct Segment {
    double slope, intercept;
    int start_idx, end_idx;
};

class DynamicPGM {
public:
    std::vector<Segment> segments;

    void insert(double key) {
        int index = find_segment(key);
        if (!can_fit(segments[index], key)) {
            split_segment(index, key);
        }
    }

    void remove(double key) {
        int index = find_segment(key);
        if (can_merge(index)) {
            merge_segments(index);
        }
    }

private:
    int find_segment(double key) {
        for (size_t i = 0; i < segments.size(); i++) {
            if (key >= segments[i].start_idx && key <= segments[i].end_idx) {
                return i;
            }
        }
        return -1;
    }

    bool can_fit(const Segment& seg, double key) {
        // Check if key fits within segment bounds
        return true;
    }

    void split_segment(int index, double key) {
        Segment new_segment = {segments[index].slope, segments[index].intercept, static_cast<int>(key), segments[index].end_idx};
        segments[index].end_idx = static_cast<int>(key);
        segments.insert(segments.begin() + index + 1, new_segment);
    }

    bool can_merge(int index) {
        // Logic to check if merging is possible
        return true;
    }

    void merge_segments(int index) {
        segments[index].end_idx = segments[index + 1].end_idx;
        segments.erase(segments.begin() + index + 1);
    }
};