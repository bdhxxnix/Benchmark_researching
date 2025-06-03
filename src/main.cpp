#include "paraoptimal.hpp"
#include <iostream>
#include <iomanip>
#include <vector>

int main() {
    // Create a ParaOptimal instance with error bound of 0.5
    plr::ParaOptimal plr(0.5);
    
    // Example data points (timestamp, value)
    std::vector<std::pair<double, double>> data = {
        {1.0, 2.0},
        {2.0, 3.0},
        {3.0, 4.0},
        {4.0, 6.0},
        {5.0, 7.0},
        {6.0, 8.0},
        {7.0, 10.0},
        {8.0, 11.0}
    };
    
    // Process each point
    for (const auto& [x, y] : data) {
        plr.processPoint(plr::Point(x, y));
    }
    
    // Print the resulting line segments
    std::cout << "Piecewise Linear Representation:\n";
    std::cout << std::fixed << std::setprecision(2);
    
    for (const auto& segment : plr.getSegments()) {
        std::cout << "Segment: y = " << segment.a << "x + " << segment.b
                  << " (points " << segment.start_idx << " to " << segment.end_idx << ")\n";
    }
    
    return 0;
} 