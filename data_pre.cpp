// preprocess the data
/* Write your docstring here*/
/**
* @file Red-Balck_trees.cpp
* @author Jiayong Qin (2057729401@qq.com)
* @brief 
* @version 0.1
* @date 2024-12-18
* 
* @copyright Copyright (c) 2024
* 
*/
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <fstream>

using namespace std;
struct Point {
    /**
     * @brief Data used to represent a point in 2D space
     * 
     */
    double x, y;
};


/* The code below is often used for testing */
int main() {
    // Your code here
    std::ifstream file("data1.txt");
    std::ofstream out("data2.txt");
    std::string line;
    int i = 0;
    int x;
    while(std::getline(file, line)){
        i++;
        file >> x ;
        out << x << " " << i << std::endl;
    }
    
    return 0;
}