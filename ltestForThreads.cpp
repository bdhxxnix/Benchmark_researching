
#include "greedy.hpp" 
#include "optimal.hpp"
#include "swing.hpp"
#include "FRS.hpp"

#include <iostream>
#include <vector>
#include <deque>
#include <limits>
#include <algorithm>
#include <fstream>
#include <chrono>
#include <string> 


/**
 * @brief Translate the data from the file to the vector
 *
 * @param filename
 * @return std::vector<K>
 */
template<typename K>
std::vector<K> load_data(std::string filename)
{

    /* Open file. */
    std::ifstream in(filename, std::ios::binary);
    if (!in.is_open())
        exit(EXIT_FAILURE);

    /* Read number of keys. */
    K n_keys;
    in.read(reinterpret_cast<char *>(&n_keys), sizeof(K));

    /* Initialize vector. */
    std::vector<K> data;
    data.resize(n_keys);

    /* Read keys. */
    in.read(reinterpret_cast<char *>(data.data()), n_keys * sizeof(K));
    in.close();

    /* Sort the data in increasing order. */
    return data;
}

template<typename K>
/**
 * @brief experiment with the GreedyPLR
 * 
 * @param data 
 * @param epsilon 
 * @param threads_num 
 * @return size_t of the GreedyPLR
 */
std::vector<size_t> experiment_Greedy(std::vector<K> data,size_t epsilon = 32,int threads_num = 16){

    std::vector<size_t> increments;
    // Initialize the results

    // Initialize the GreedyPiecewiseLinearModel 
    using Segment= typename Greedy::internal::GreedyPiecewiseLinearModel<K,int>::CanonicalSegment::Segment;
    auto result_segments_serial = std::vector<Segment>();
    auto result_segments_parallel = std::vector<Segment>();
    
    // Initialize the epsilon
    // Initialize the input and output function
    auto in = [&](size_t i) { return data[i]; };
    auto out_serial = [&result_segments_serial](const auto &cs) { 
        auto segment = cs.get_Canonicalsegment(0);
        result_segments_serial.push_back(segment); 
    };
    auto out_parallel = [&result_segments_parallel](const auto &cs) { 
        auto segment = cs.get_Canonicalsegment(0);
        result_segments_parallel.push_back(segment); 
    };

    size_t num_segments = Greedy::internal::make_segmentation(data.size(), epsilon, in, out_serial);
    // Print out the result of the serial version
    for(int i = 1 ; i< 30 ; i+=2 ){
        size_t c = Greedy::internal::make_segmentation_par(data.size(), epsilon, in, out_parallel,i);
        increments.push_back(c - num_segments);
    }
    
    printf("The result of GreedyPLA: ");
    for(size_t i = 0;i<increments.size();i++){
        printf("%ld,",increments[i]);
    }
    printf("\n");

    return increments; 
}

template<typename K>
/**
 * @brief experiment with the SwingFilter
 * 
 * @param data 
 * @param epsilon 
 * @param threads_num 
 * @return size_t of the SwingFIlter
 */
std::vector<size_t> experiment_Swing(std::vector<K> data,size_t epsilon = 32,int threads_num = 16){

    std::vector<size_t> increments;
    // Initialize the results

    // Initialize the SwingPiecewiseLinearModel 
    using Segment= typename Swing::internal::SwingPiecewiseLinearModel<K,int>::CanonicalSegment::Segment;
    auto result_segments_serial = std::vector<Segment>();
    auto result_segments_parallel = std::vector<Segment>();
    
    // Initialize the epsilon
    // Initialize the input and output function
    auto in = [&](size_t i) { return data[i]; };
    auto out_serial = [&result_segments_serial](const auto &cs) { 
        auto segment = cs.get_Canonicalsegment(0);
        result_segments_serial.push_back(segment); 
    };
    auto out_parallel = [&result_segments_parallel](const auto &cs) { 
        auto segment = cs.get_Canonicalsegment(0);
        result_segments_parallel.push_back(segment); 
    };

    size_t num_segments = Swing::internal::make_segmentation(data.size(), epsilon, in, out_serial);
    // Print out the result of the serial version
    for(int i = 1 ; i<30 ; i+=2 ){
        size_t c = Swing::internal::make_segmentation_par(data.size(), epsilon, in, out_parallel,i);
        increments.push_back(c - num_segments);
    }
    // Swing::internal::checkForEpsilon(data.size(),in,result_segments_serial,0,result.seg_serial-1,epsilon);  

    printf("The result of SwingFilter: ");
    for(size_t i = 0;i<increments.size();i++){
        printf("%ld,",increments[i]);
    }
    printf("\n");

    return increments; 
}



/**
 * @brief Experiment with the OptimalPLR (PGM-Index)
 * 
 * @param data 
 * @param epsilon 
 * @param threads_num 
 * @return size_t of the OptimalPLR
 *
 * 
 */
template<typename K>
std::vector<size_t> experiment_Optimal(std::vector<K> data,size_t epsilon = 32,int threads_num = 16){

    // Initialize the results
    std::vector<size_t> increments;

    // Initialize the OptimalPiecewiseLinearModel  
    std::vector<std::pair<double, double>> serial_segments;
    std::vector<std::pair<double, double>> parallel_segments;
    
    // Initialize the epsilon
    // Initialize the input and output function
    auto in = [&data](size_t i) { return data[i]; };
    auto out_serial = [&serial_segments](const auto &segment) { 
        auto result = segment.get_floating_point_segment(0);
        serial_segments.push_back(result);
    };
    auto out_parallel = [&parallel_segments](const auto &segment) { 
        auto result = segment.get_floating_point_segment(0);
        parallel_segments.push_back(result);
    };

    size_t num_segments = Optimal::internal::make_segmentation(data.size(), epsilon, in, out_serial);
    for(int i = 1;i<30;i+=2){
        size_t c = Optimal::internal::make_segmentation_par(data.size(), epsilon, in, out_parallel,i);
        increments.push_back(c-num_segments);
    } 
    
    printf("The result of OptimalPLA: ");
    for(size_t i = 0;i<increments.size();i++){
        printf("%ld,",increments[i]);
    }
    printf("\n");

    return increments; 
}


int main(int argc ,char* argv[]){

    if (argc < 2) {
        std::cerr << "Usage: ./fitting_tree_test <dataset_path>\n";
        return 1;
    }

    std::string file_path = argv[1];
    // Read the data from data1.txt
    std::vector<uint64_t> data;
    data = load_data<uint64_t>(file_path);
    printf("We are running on linear test with varying threads\n");
    printf("The size of data is %ld\n",data.size());
    
    experiment_Greedy(data, 4);
    experiment_Swing(data,4);
    experiment_Optimal(data,4);
    return 0;
}

