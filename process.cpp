#include "optimalPLR_Foundation.h"
#include <string>
#include <unordered_map>

#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <algorithm>

/**
 * @brief read the data from the file
 * @param filename the name of the file
 * @param points the container to store the points
 * @return the information about the dataset, such as the size and epsilon
 */
std::unordered_map<std::string, double> BasicInfoAboutFile(std::string filename,std::vector<Point> &points){
    
    std::unordered_map<std::string, double> info;
    std::ifstream file(filename);
    
    int time = 0;
    double maxi = -std::numeric_limits<double>::max();
    double mini = std::numeric_limits<double>::max();
    while (!file.eof())
    {
        time++;
        Point p;
        p.x = time;
        file >> p.y;
        points.push_back(p);
        maxi = std::max(maxi, p.y);
        mini = std::min(mini, p.y);
    }
    
    double range = maxi - mini;
    double epsilon = range * 0.05;
    info["size"] = points.size();

    info["epsilon"] = epsilon;

    return info;
}


/**
 * @brief write the output to a file
 * @param info the information about the dataset
 * @param num_segments the number of segments
 */
void output(std::unordered_map<std::string, double> info,size_t num_segments,std::string filename){
    // Write the output to a file
    std::ofstream outfile(filename);

    outfile  << info["size"] << "\t" << info["epsilon"] << "\t" << num_segments << "\n";
    
    outfile.close();
}

