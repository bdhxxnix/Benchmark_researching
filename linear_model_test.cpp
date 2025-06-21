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
 * @brief The result of a certain experiment (augmented with the time taken)
 * 
 */
struct Result
{
    int seg_serial,seg_parallel; double time_serial,time_parallel;

    // Initialize the result
    Result(int seg_serial,int seg_parallel,double time_serial,double time_parallel):
        seg_serial(seg_serial),seg_parallel(seg_parallel),time_serial(time_serial),time_parallel(time_parallel){}
    Result():seg_serial(0),seg_parallel(0),time_serial(0),time_parallel(0){}
    Result(const Result &other):seg_serial(other.seg_serial),seg_parallel(other.seg_parallel),time_serial(other.time_serial),time_parallel(other.time_parallel){}
};



/**
 * @brief Translate the data from the file to the vector
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

/**
 * @brief Calculate the results of the experiment
 * @param results 
 */
void print(std::vector<Result> results){
    
    std::cout << "#Segments (Serial): "<<std::endl;
    for(size_t i = 0;i<results.size();i++){
        std::cout << results[i].seg_serial  << ",";
    }
    std::cout << std::endl;
    std::cout << "#Segments (Parallel): ";
    for(size_t i = 0;i<results.size();i++){
        std::cout << results[i].seg_parallel  << ",";
    }
    std::cout << std::endl;
    std::cout << "Time (Serial): ";
    for(size_t i = 0;i<results.size();i++){
        std::cout << results[i].time_serial  << ",";
    }
    std::cout << std::endl;
    std::cout << "Time (Parallel): ";
    for(size_t i = 0;i<results.size();i++){
        std::cout << results[i].time_parallel  << ",";
    }
}

/**
 * @brief printout the #segments and the time taken for the serial and parallel version 
 *  
 */
void printout(std::vector<Result> results_greedy,
                std::vector<Result>results_swing,
                std::vector<Result>results_linear,
            std::vector<Result> results_FRS){

    
    // Calculateout the result of the GreedyPLR
    std::cout << "OptimalPLA :" << std::endl;
    print(results_linear);
    std::cout<<std::endl;
    std::cout << "SwingFilter:" << std::endl;
    print(results_swing);
    std::cout<<std::endl;
    std::cout << "GreedyPLA:" << std::endl;
    print(results_greedy);    
    std::cout<<std::endl;
    std::cout << "FRS:" << std::endl;
    print(results_FRS);    
    std::cout<<std::endl;

}

template<typename K>
/**
 * @brief experiment with the FRS
 * 
 * @param data 
 * @param epsilon 
 * @param threads_num 
 * @return Result of the GreedyPLR
 */
Result experiment_FRS(std::vector<K> data,size_t epsilon = 32,int threads_num = 16){

    std::vector<size_t> segNumber;
    // Initialize the results
    Result result;

    // Initialize the GreedyPiecewiseLinearModel and OptimalPiecewiseLinearModel and OptimalPLR  
    using Segment= typename FRS::internal::FRSPiecewiseLinearModel<K,int>::CanonicalSegment::Segment;
    auto result_segments_serial = std::vector<Segment>();
    auto result_segments_parallel = std::vector<Segment>();
    
    // Initialize the epsilon
    // Initialize the input and output function
    auto in = [&](size_t i) { return data[i]; };
    auto out_serial = [&result_segments_serial](const auto &cs) { 
        auto segment = cs.get_Canonicalsegment();
        result_segments_serial.push_back(segment); 
    };
    auto out_parallel = [&result_segments_parallel](const auto &cs) { 
        auto segment = cs.get_Canonicalsegment();
        result_segments_parallel.push_back(segment); 
    };
    
    // Calculate the result of the serial version
    auto start = std::chrono::high_resolution_clock::now();
    size_t num_segments = FRS::internal::make_segmentation(data.size(), epsilon,in , out_serial);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> serial_duration = end - start;
    // The time taken for the serial version and the number of the segments
    result.time_serial = serial_duration.count();
    result.seg_serial = num_segments;

    // Calculate the result of the parallel version
    int parallelism = threads_num; 
    start = std::chrono::high_resolution_clock::now();
    num_segments = FRS::internal::make_segmentation_par(data.size(), epsilon, in, out_parallel,parallelism);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> parallel_duration = end - start;
    // The time taken for the parallel version and the number of the segments
    result.time_parallel = parallel_duration.count();
    result.seg_parallel = num_segments;

    // Check the residuals for the segments
    FRS::internal::checkForEpsilon(in,result_segments_parallel,0,result.seg_parallel-1,epsilon);  

    return result; 
}



template<typename K>
/**
 * @brief experiment with the GreedyPLR
 * 
 * @param data 
 * @param epsilon 
 * @param threads_num 
 * @return Result of the GreedyPLR
 */
Result experiment_Greedy(std::vector<K> data,size_t epsilon = 32,int threads_num = 16){

    std::vector<size_t> segNumber;
    // Initialize the results
    Result result;

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


    // Calculate the result of the serial version
    auto start = std::chrono::high_resolution_clock::now();
    size_t num_segments = Greedy::internal::make_segmentation(data.size(), epsilon, in, out_serial);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> serial_duration = end - start;
    // The time taken for the serial version and the number of the segments
    result.time_serial = serial_duration.count();
    result.seg_serial = num_segments;

    // Calculate the result of the serial version
    start = std::chrono::high_resolution_clock::now();
    num_segments = Greedy::internal::make_segmentation_par(data.size(), epsilon, in, out_parallel,threads_num);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> parallel_duration = end - start;
    // The time taken for the parallel version and the number of the segments
    result.time_parallel = parallel_duration.count();
    result.seg_parallel = num_segments;
    
    Greedy::internal::checkForEpsilon(data.size(),in,result_segments_parallel,0,result.seg_parallel-2,epsilon);  
    
    return result; 
}

/**
 * @brief Experiment with the OptimalPLR (PGM-Index)
 * 
 * @param data 
 * @param epsilon 
 * @param threads_num 
 * @return Result of the OptimalPLR
 *
 * 
 */
template<typename K>
Result experiment_Optimal(std::vector<K> data,size_t epsilon = 32,int threads_num = 16){

    // Intialize the results
    Result result;
    std::vector<size_t> segNumber;

    // Initialize the OptimalPiecewiseLinearModel  
    using Segment = typename Optimal::internal::OptimalPiecewiseLinearModel<K,size_t>::CanonicalSegment;
    std::vector<Segment> serial_segments;
    std::vector<Segment> parallel_segments;
    
    // Initialize the epsilon
    // initialize the input and output function
    auto in = [&data](size_t i) { return data[i]; };
    auto out_serial = [&serial_segments](const auto &segment) { 
        auto result = segment;
        serial_segments.push_back(result);
    };
    auto out_parallel = [&parallel_segments](const auto &segment) { 
        auto result = segment;
        parallel_segments.push_back(result);
    };

    auto start = std::chrono::high_resolution_clock::now();
    size_t num_segments = Optimal::internal::make_segmentation(data.size(), epsilon, in, out_serial);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> serial_duration = end - start;
    result.seg_serial = num_segments;
    result.time_serial = serial_duration.count();
    
    // The time taken for the parallel version and the number of the segments
    start = std::chrono::high_resolution_clock::now();
    num_segments = Optimal::internal::make_segmentation_par(data.size(), epsilon, in, out_parallel,threads_num);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> parallel_duration = end - start;
    // The time taken for the parallel version and the number of the segments
    result.time_parallel = parallel_duration.count();
    result.seg_parallel = num_segments;
    
    Optimal::internal::checkForEpsilon(data.size(),in,parallel_segments,0,result.seg_parallel-3,epsilon);  
    return result; 
}
template<typename K>
/**
 * @brief experiment with the SwingFilter
 * 
 * @param data 
 * @param epsilon 
 * @param threads_num 
 * @return Result of the SwingFilter
 */
Result experiment_Swing(std::vector<K> data,size_t epsilon = 32,int threads_num = 16){

    std::vector<size_t> segNumbers; 
    // Initialize the results
    Result result;

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

    // Start the experiment of the SwingFilter
    // Calculate the result of the serial version
    auto start = std::chrono::high_resolution_clock::now();
    size_t num_segments = Swing::internal::make_segmentation(data.size(), epsilon, in, out_serial);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> serial_duration = end - start;
    result.seg_serial = num_segments;
    result.time_serial = serial_duration.count();
    
    // Calculate the result of the parallel version
    start = std::chrono::high_resolution_clock::now();
    num_segments = Swing::internal::make_segmentation_par(data.size(), epsilon, in, out_parallel,threads_num);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> parallel_duration = end - start;
    // The time taken for the parallel version and the number of the segments
    result.time_parallel = parallel_duration.count();
    result.seg_parallel = num_segments;
    segNumbers.push_back(num_segments);
    
    Swing::internal::checkForEpsilon(data.size(),in,result_segments_parallel,0,result.seg_parallel-1,epsilon);  
    
    return result; 
}

template<typename K>
/**
 * 
 * @brief Experiment with the greedyPLR and OptimalPLR (PGM-Index & Customized) 
 * 
 * @param data 
 * @param epsilon 
 * @param threads_num 
 */
void experiment(std::vector<K> data,int start, int end,int threads_num = 32)
{
    // Experiment with the greedyPLR and OptimalPLR (PGM-Index & Customized)
    // Initialize the results
    std::vector<Result> results_greedy;
    std::vector<Result> results_linear;
    std::vector<Result> results_swing;
    std::vector<Result> results_FRS;

    // Run the experiment for the greedyPLR
    for(int i = start;i<end;i++){
        double epsilon = 1<<i;
        printf("The epsilon is %f\n",epsilon);
        
        // Run the experiment for the greedyPLR
        Result result_greedy = experiment_Greedy(data,epsilon,threads_num);
        results_greedy.push_back(result_greedy); 
        
        // Run the experiment for the OptimalPLR (Customized)
        Result result_optimal = experiment_Swing(data,epsilon,threads_num);
        results_swing.push_back(result_optimal);
        
        // Run the experiment for the OptimalPLR (PGM-Index)
        Result result_linear = experiment_Optimal(data,epsilon,threads_num);
        results_linear.push_back(result_linear);

        // Run the experiment for the FRS
        Result result_FRS = experiment_FRS(data,epsilon,threads_num);
        results_FRS.push_back(result_FRS);
    }

    printout(results_greedy,results_swing,results_linear,results_FRS);

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
    
    printf("We are running on linear test\n");
    printf("The size of data is %ld\n",data.size());
    
    experiment(data, 2, 14);
    return 0;
}

