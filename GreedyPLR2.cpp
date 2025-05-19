#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <limits>
#include <chrono>

struct Point {
    double x, y;
};

struct Segment {
    double slope, intercept;
    int start_idx, end_idx;
};

// Greedy PGM segmentation with minimal memory usage
std::vector<Segment> segmentPGM(const std::vector<Point>& points, double epsilon) {
    std::vector<Segment> segments;
    if (points.size() < 2) return segments;

    int start = 0;
    double s_upper = std::numeric_limits<double>::infinity();

    for (size_t i = 1; i < points.size(); ++i) {
        double dx = points[i].x - points[start].x;
        if (dx == 0) continue;  // Skip duplicate x-values

        double upper_bound = (points[i].y + epsilon - points[start].y) / dx;

        // If upper slope is decreasing, segment must end
        if (upper_bound < s_upper) {
            int end_idx = i - 1;
            
            // Compute slope and intercept from midpoint
            int mid_idx = (start + end_idx) / 2;
            double final_slope = s_upper;
            double intercept = points[mid_idx].y - final_slope * points[mid_idx].x;

            // Store minimal segment data
            segments.push_back({final_slope, intercept, start, end_idx});

            // Start new segment
            start = end_idx;
            s_upper = std::numeric_limits<double>::infinity();
            i = start;  // Restart processing from new segment
        } else {
            s_upper = upper_bound;  // Update valid upper bound
        }
    }

    // Final segment
    int end_idx = points.size() - 1;
    int mid_idx = (start + end_idx) / 2;
    double final_slope = s_upper;
    double intercept = points[mid_idx].y - final_slope * points[mid_idx].x;
    segments.push_back({final_slope, intercept, start, end_idx});

    return segments;
}

int main() {
    std::ifstream file("data2.txt");
    if (!file) {
        std::cerr << "Error: Cannot open file.\n";
        return 1;
    }

    std::vector<Point> points;
    const size_t CHUNK_SIZE = 1000000;  // Process 1M points at a time
    double x, y;
    double epsilon = 2.0;
    int ser = 0;

    while (ser<=20) {
        ser++;
        points.clear();
        points.reserve(CHUNK_SIZE);

        // Read a chunk of data
        for (size_t i = 0; i < CHUNK_SIZE && file >> x >> y; ++i) {
            points.push_back({x, y});
        }
        if (points.empty()) break;

        auto start1 = std::chrono::high_resolution_clock::now();
        // Process chunk
        auto segments = segmentPGM(points, epsilon);

        auto end1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration1 = end1 - start1;
    std::cout << "The size of the segments is " << segments.size() << std::endl;
    std::cout << "Time taken for cycle"<<ser<<": " << duration1.count() << " seconds\n";
    }

    std::cout << "Segmentation completed.\n";
    return 0;
}
