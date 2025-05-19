#include "greedy.hpp"
#include <iostream>
#include <vector>
#include <deque>
#include <limits>
#include <algorithm>
#include <fstream>
#include <chrono>


template <typename SegmentType>
void checkforepi_greplr(std::vector<double> data, std::vector<SegmentType> segments){
    
    auto start = segments[0].first_x;
    auto end = segments[10].last_x;
    
    int segmnt_idx = 0;
    auto seg = segments[segmnt_idx];
    auto slope = seg.slope;
    auto intercept = seg.intercept;
    auto first_x = seg.first_x;
    auto last_x = seg.last_x;

    int i = 0;
    auto max_residual = std::numeric_limits<double>::min();
    while (data[i] <= end)
    {
        if (data[i] == last_x)
        {
            // printout the max_residual
            printf("The max_residual for Segment %d is %f\n",segmnt_idx,max_residual);
            max_residual = std::numeric_limits<double>::min();
           
            segmnt_idx++;
            if(segmnt_idx>=segments.size()){
                break;
            }
            seg = segments[segmnt_idx];
            slope = seg.slope;
            intercept = seg.intercept;
            first_x = seg.first_x;
            last_x = seg.last_x;
        }

        auto residual = std::abs(i - (slope*data[i] + intercept));
        if(residual>max_residual){
            max_residual = residual;
        }
        i++;
    }
}

void printout(std::vector<int>seg_serial,std::vector<int>seg_parallel,std::vector<float>time_serial,std::vector<float>time_parallel){
    for(int i = 0;i<seg_serial.size();i++){
        std::cout<<seg_serial[i] << "," ;
    }
    std::cout<<std::endl;
    for(int i = 0;i<seg_parallel.size();i++){
        std::cout<<seg_parallel[i] << ",";
    }
    std::cout<<std::endl;

    for(int i = 0;i<time_serial.size();i++){
        std::cout<<time_serial[i] << "," ;
    }
    std::cout<<std::endl;
    for(int i = 0;i<time_parallel.size();i++){
        std::cout<<time_parallel[i] << ",";
    }
}


int main(){
    
    // Read the data from data1.txt
    std::ifstream file("../data/data1.txt");
    std::vector<double> data(20000000);
    
    // double maxi = -std::numeric_limits<double>::max();
    // double mini = std::numeric_limits<double>::max();
    for(int i = 0;i<20000000;i++){
        file >> data[i];
    }

    

    // //Check the origin data
    // for(int i = 0;i<=10;i++){
    //     std::cout << "Data " << i << " : " << data[i] << std::endl;
    // }

    std::vector<float> time_serial;
    std::vector<int> seg_serial;
    std::vector<float> time_parallel;
    std::vector<int> seg_parallel;

    using Segment = GreedyPiecewiseLinearModel<double,int>::CanonicalSegment::Segment;
    
    for(int i=5;i<11;i++){
        double epsilon = 1<<i;
        printf("The epsilon is %f\n",epsilon);
        
        auto result_segments_serial = std::vector<Segment>();
        auto result_segments_parallel = std::vector<Segment>();
        // Initialize the epsilon
        
        // Initialize the GreedyPiecewiseLinearModel
        auto in = [&](size_t i) { return data[i]; };
        auto out_serial = [&result_segments_serial](const auto &cs) { 
            auto segment = cs.get_Canonicalsegment();
            result_segments_serial.push_back(segment); 
        };
        auto out_parallel = [&result_segments_parallel](const auto &cs) { 
            auto segment = cs.get_Canonicalsegment();
            result_segments_parallel.push_back(segment); 
        };


        // Compare the time taken for the serial and parallel version
        auto start = std::chrono::high_resolution_clock::now();
        size_t num_segments = make_segmentation(data.size(), epsilon, in, out_serial);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> serial_duration = end - start;
        time_serial.push_back(serial_duration.count());
        seg_serial.push_back(num_segments);
       

        start = std::chrono::high_resolution_clock::now();
        num_segments = make_segmentation_par(data.size(), epsilon, in, out_parallel);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> parallel_duration = end - start;
        time_parallel.push_back(parallel_duration.count());
        seg_parallel.push_back(num_segments);

           
        // checkforepi(data,result_segments_serial);
    }

    // Print the results
    


    printout(seg_serial,seg_parallel,time_serial,time_parallel);    
    
    return 0;
}
