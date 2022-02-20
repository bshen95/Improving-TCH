//
// Created by Bojie Shen on 26/5/21.
//

#include <iostream>
#include <stdlib.h>
#include <string>
#include <iomanip>

#include "io/edge_info.h"
#include "io/tpgr_format.h"
#include "prepro/ordering.h"
#include "util/misc.h"
#include "datastr/graph/search_graph.h"
// it seems like we could set hourly time domain equal to 864000 without losing too
// much performance, since my data is hourly data. Anyway, for larger map we have to change this.


static constexpr double offset = 5000;
//static constexpr double one_TCH_period = 113000.0;
//static constexpr double two_TCH_period = 149000.0;
//static constexpr double four_TCH_period = 221000.0;

//static constexpr double one_TCH_period = 77000.0;
//static constexpr double two_TCH_period = 113000.0;
//static constexpr double four_TCH_period = 185000.0;
static constexpr double STCH_period = 36000.0 * 5 + offset ;
static constexpr double one_MTCH_period = 36000.0 *2 + offset;
static constexpr double two_MTCH_period = 36000.0 *3 + offset;
static constexpr double four_MTCH_period = 36000.0 *5 + offset;

template<const double* h_period>
int construct_Hourly_TCH(const std::string& tpgr_input_file_name,  const int& n_threads){

    using Graph = typename katch::Ordering<h_period>::Graph;
    using EdgeInfo = katch::EdgeInfo<typename katch::Ordering<h_period>::TTF>;
    for(unsigned hour = 0; hour < 24; hour ++) {
        auto t1 = katch::util::time_stamp();
        std::vector<EdgeInfo> edge_list = katch::tpgr_format::read_edges<EdgeInfo>(
                tpgr_input_file_name + "_" + std::to_string(hour) + ".tpgr");

        if (edge_list.empty()) {
            KATCH_ERROR("Empty graph.\n");
            return EXIT_FAILURE;
        }

        Graph graph(std::move(edge_list));
        KATCH_STATUS("Graph has " << graph.get_n_nodes() << " nodes, " << graph.get_n_edges() << " edges.\n");

        katch::Ordering<h_period> ordering(std::move(graph), tpgr_input_file_name + "_" + std::to_string(hour) + ".btch");
        std::cout << "here" << std::endl;
        ordering.order_and_construct(n_threads);
        KATCH_STATUS("Done.\n");
        auto t2 = katch::util::time_stamp();

        std::vector<EdgeInfo> edge_list2 = katch::btch_format::read_edges<EdgeInfo>(
                tpgr_input_file_name + "_" + std::to_string(hour) + ".btch", *h_period);
        katch::SearchGraph<h_period> htch_graph(std::move(edge_list2));

        KATCH_STATUS("HTCH size " + std::to_string((double) htch_graph.get_index_size() / 1000000) + " MB.\n");
        std::string main_name = tpgr_input_file_name.substr(0, tpgr_input_file_name.find_last_of("/"));
        std::ofstream myFile(main_name + "/results/preprocessing/htch_preprocessing_" + std::to_string(hour) + ".csv");
        myFile << "build_time,memory_cost\n";
        myFile << std::fixed << std::setprecision(8) << katch::util::get_duration_in_seconds(t1, t2) / 60
               << "," << (double) htch_graph.get_index_size() / 1000000
               << "\n";
        myFile.close();
    }
    return EXIT_SUCCESS;
}


template<const double* h_period>
int construct_Hourly_TCH(int hour, const std::string& tpgr_input_file_name,  const int& n_threads){

    using Graph = typename katch::Ordering<h_period>::Graph;
    using EdgeInfo = katch::EdgeInfo<typename katch::Ordering<h_period>::TTF>;
        auto t1 = katch::util::time_stamp();
        std::vector<EdgeInfo> edge_list = katch::tpgr_format::read_edges<EdgeInfo>(
                tpgr_input_file_name + "_" + std::to_string(hour) + ".tpgr");

        if (edge_list.empty()) {
            KATCH_ERROR("Empty graph.\n");
            return EXIT_FAILURE;
        }

        Graph graph(std::move(edge_list));
        KATCH_STATUS("Graph has " << graph.get_n_nodes() << " nodes, " << graph.get_n_edges() << " edges.\n");

        katch::Ordering<h_period> ordering(std::move(graph), tpgr_input_file_name + "_" + std::to_string(hour) + ".btch");
        std::cout << "here" << std::endl;
        ordering.order_and_construct(n_threads);
        KATCH_STATUS("Done.\n");
        auto t2 = katch::util::time_stamp();

        std::vector<EdgeInfo> edge_list2 = katch::btch_format::read_edges<EdgeInfo>(
                tpgr_input_file_name + "_" + std::to_string(hour) + ".btch", *h_period);
        katch::SearchGraph<h_period> htch_graph(std::move(edge_list2));

        KATCH_STATUS("HTCH size " + std::to_string((double) htch_graph.get_index_size() / 1000000) + " MB.\n");
        std::string main_name = tpgr_input_file_name.substr(0, tpgr_input_file_name.find_last_of("/"));
        std::ofstream myFile(main_name + "/results/preprocessing/htch_preprocessing_" + std::to_string(hour) + ".csv");
        myFile << "build_time,memory_cost\n";
        myFile << std::fixed << std::setprecision(8) << katch::util::get_duration_in_seconds(t1, t2) / 60
               << "," << (double) htch_graph.get_index_size() / 1000000
               << "\n";
        myFile.close();
    return EXIT_SUCCESS;
}



template<const double* h_period>
int construct_n_Hourly_TCH(int hour_duration, const std::string& tpgr_input_file_name,  const int& n_threads){
    int shifting_hour = *h_period/36000 - hour_duration;
    using Graph = typename katch::Ordering<h_period>::Graph;
    using EdgeInfo = katch::EdgeInfo<typename katch::Ordering<h_period>::TTF>;
    for(unsigned hour = 0; hour < 24; hour += hour_duration) {
        auto t1 = katch::util::time_stamp();
        std::vector<EdgeInfo> edge_list = katch::tpgr_format::read_edges<EdgeInfo>(
                tpgr_input_file_name + "_" + std::to_string(hour_duration) + "_" + std::to_string(shifting_hour) + "_" + std::to_string(hour)+ ".tpgr");

        if (edge_list.empty()) {
            KATCH_ERROR("Empty graph.\n");
            return EXIT_FAILURE;
        }

        Graph graph(std::move(edge_list));
        KATCH_STATUS("Graph has " << graph.get_n_nodes() << " nodes, " << graph.get_n_edges() << " edges.\n");

        katch::Ordering<h_period> ordering(std::move(graph), tpgr_input_file_name + "_" + std::to_string(hour_duration) + "_" + std::to_string(shifting_hour) + "_" + std::to_string(hour)+ ".btch");
        std::cout << "here" << std::endl;
        ordering.order_and_construct(n_threads);
        KATCH_STATUS("Done.\n");
        auto t2 = katch::util::time_stamp();

        std::vector<EdgeInfo> edge_list2 = katch::btch_format::read_edges<EdgeInfo>(
                tpgr_input_file_name + "_" + std::to_string(hour_duration) + "_" + std::to_string(shifting_hour) + "_" + std::to_string(hour)+ ".btch", *h_period);
        katch::SearchGraph<h_period> htch_graph(std::move(edge_list2));

        KATCH_STATUS("HTCH size " + std::to_string((double) htch_graph.get_index_size() / 1000000) + " MB.\n");
        std::string main_name = tpgr_input_file_name.substr(0, tpgr_input_file_name.find_last_of("/"));
        std::ofstream myFile(main_name + "/results/preprocessing/htch_preprocessing_" + std::to_string(hour_duration) + "_" + std::to_string(shifting_hour) + "_" + std::to_string(hour) + ".csv");
        myFile << "build_time,memory_cost\n";
        myFile << std::fixed << std::setprecision(8) << katch::util::get_duration_in_seconds(t1, t2) / 60
               << "," << (double) htch_graph.get_index_size() / 1000000
               << "\n";
        myFile.close();
    }
    return EXIT_SUCCESS;
}























int main(int argc, char** argv)
{

//    using Graph = katch::Ordering<&period>::Graph;
//    using EdgeInfo = katch::EdgeInfo<katch::Ordering<&period>::TTF>;

    const char* binary_name = argv[0];

    if ( argc != 4 )
    {
        std::cerr
                << std::endl
                << "USAGE: " << binary_name
                << " <.tpgr input file> <hour> <n threads>" << std::endl
                << std::endl;

        return EXIT_FAILURE;
    }

    const std::string tpgr_input_file_name(argv[1]);
    const int n_threads = std::stoi(std::string(argv[3]));

    // construct STCH
    construct_Hourly_TCH<&STCH_period>(tpgr_input_file_name,n_threads);

    // construct MTCH
    construct_n_Hourly_TCH<&one_MTCH_period>(1,tpgr_input_file_name,n_threads);
    construct_n_Hourly_TCH<&two_MTCH_period>(1,tpgr_input_file_name,n_threads);
    construct_n_Hourly_TCH<&four_MTCH_period>(1,tpgr_input_file_name,n_threads);

//
//    auto t1 = katch::util::time_stamp();
//    std::vector<EdgeInfo> edge_list = katch::tpgr_format::read_edges<EdgeInfo>(tpgr_input_file_name + "_" + hour+".tpgr");
//
//    if ( edge_list.empty() )
//    {
//        KATCH_ERROR("Empty graph.\n");
//        return EXIT_FAILURE;
//    }
//
//    Graph graph(std::move(edge_list));
//    KATCH_STATUS("Graph has " << graph.get_n_nodes() << " nodes, " << graph.get_n_edges() << " edges.\n");
//
//    katch::Ordering<&period> ordering(std::move(graph), tpgr_input_file_name+"_" + hour+".btch");
//    std::cout<<"here"<<std::endl;
//    ordering.order_and_construct(n_threads);
//    KATCH_STATUS("Done.\n");
//    auto t2 = katch::util::time_stamp();
//
//    std::vector<EdgeInfo> edge_list2 = katch::btch_format::read_edges<EdgeInfo>(tpgr_input_file_name+"_" + hour+".btch",185000.0);
//    katch::SearchGraph<&period> htch_graph(std::move(edge_list2));
//
//    KATCH_STATUS("HTCH size "+ std::to_string((double)htch_graph.get_index_size()/1000000) + " MB.\n");
//    std::string main_name = tpgr_input_file_name.substr(0,tpgr_input_file_name.find_last_of("/"));
//    std::ofstream myFile(main_name+"/results/preprocessing/htch_preprocessing_"+hour+".csv");
//    myFile<<"build_time,memory_cost\n";
//    myFile<<std::fixed<<std::setprecision(8)<<katch::util::get_duration_in_seconds(t1, t2)/60
//          <<","<<(double)htch_graph.get_index_size()/1000000
//          <<"\n";
//    myFile.close();


    return EXIT_SUCCESS;
}