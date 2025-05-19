#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>

int main (){
    std::ifstream file("../data/osmdata.txt");
    int first;
    file >> first;  
    int value;
    std::vector<int> data;
    // Calculate the gap between the data
    long double sum = 0; 
    while(file>> value){
       data.push_back(value-first);
       sum += value-first;
       first = value;
        
    } 
    long double average = sum / data.size();
    long double variance = 0;
    for (int i = 0; i < data.size(); i++){
        variance += (data[i] - average) * (data[i] - average);
    }
    variance /= data.size();
    long double stddev = sqrt(variance);
    std::cout << "Average: " << average << std::endl;
    std::cout << "Standard Deviation: " << stddev << std::endl;
    std::cout << "Variance: " << variance << std::endl;
    std::cout << "Number of data points: " << data.size() << std::endl; 
    std::cout << "First data point: " << data[0] << std::endl;
}