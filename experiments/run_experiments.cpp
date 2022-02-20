//
// Created by Bojie Shen on 16/5/21.
//

#include <iostream>
#include <stdlib.h>
#include <string>
#include <util/my_timer.h>
#include <io/btch_format.h>
#include <iomanip>
#include <cpd_prepro/tch_dijkstra_query.h>


#include "util/misc.h"
#include "io/tpgr_format.h"
#include "io/edge_info.h"
#include "io/vec_io.h"

// original BTCH
#include "query/BTCH.h"
// BTCH + landmark
#include "query/BTCH_L.h"
// FTCH + landmark
#include "query/FTCH_L.h"
// FTCH + TCPD
#include "query/FTCH_TCPD.h"
// FTCH + RTPD
#include "query/FTCH_RTPD.h"
// FTCH + CPD
#include "query/FTCH_CPD.h"
// B-STCH
#include "query/BSTCH.h"
// B-MTCH
#include "query/BMTCH_query.h"
#include "query/BMTCH.h"
// B-STCH + TCPD
#include "query/STCH_TCPD.h"
// B-MTCH + TCPD
#include "query/FMTCH_TCPD_query.h"
#include "query/FMTCH_TCPD.h"


std::vector<int> source;
std::vector<int> target;
std::vector<double>euclidean_distance;
int run_times = 10;
int test_period = 36000;

static constexpr double period = 864000.0;


static constexpr double offset = 5000;
static constexpr double hourly_TCH_period = 36000.0 * 5 + offset;
//static constexpr double one_TCH_period = 113000.0;
//static constexpr double two_TCH_period = 149000.0;
//static constexpr double four_TCH_period = 221000.0;

static constexpr double one_TCH_period = 36000.0 *2 + offset;
static constexpr double two_TCH_period = 36000.0 *3 + offset;
static constexpr double four_TCH_period = 36000.0 *5 + offset;

//static constexpr double hourly_TCH_period = 864000.0;
//static constexpr double one_TCH_period = 864000.0;
//static constexpr double two_TCH_period =  864000.0;
//static constexpr double four_TCH_period =  864000.0;
using S_Graph = katch::TchDijkstraQuery::Graph;
using F_Graph = katch::ForwardSearchGraph<&period>;

template<typename query_type>
void run_search( query_type query, const std::string& result_file,double departure_time){
    std::cout<< "Starting time: " << (double)departure_time/(60*60*10)<<std::endl;


    //record the best, worst, total time;
    std::vector<double> search_time = std::vector<double>(source.size(),0);
    std::vector<double> best_time = std::vector<double>(source.size());
    std::fill(best_time.begin(), best_time.end(), INFINITY);
    std::vector<double> worst_time = std::vector<double>(source.size());

    std::vector<double> path_recover_time = std::vector<double>(source.size(),0);

    std::fill(worst_time.begin(),worst_time.end(), 0);

    std::vector<unsigned> path_length =std::vector<unsigned>(source.size());
    std::vector<double> time_cost = std::vector<double>(source.size(),0);

    std::vector<unsigned long> number_of_nodes_generated = std::vector<unsigned long>(source.size());
    std::vector<unsigned long> number_of_nodes_expanded = std::vector<unsigned long>(source.size());
    std::vector<unsigned long> number_of_first_move = std::vector<unsigned long>(source.size());
    std::vector<unsigned long> number_of_bw_first_move = std::vector<unsigned long>(source.size());
    std::vector<unsigned long> number_of_TD_calculation = std::vector<unsigned long>(source.size());

    //warm up !!!!
    for(int i  = 0; i < source.size(); i ++) {
        int random_source = source[i];
        int random_target = target[i];
        query->load_target_row(departure_time);
        double arr_t = query->one_to_one(random_source,random_target,departure_time);
        const auto& original_path = query->get_path(random_target);
        size_t size = original_path.get_n_edges();
//            original_path.print_nodes();
        if(size == 0){
            //all the queries generated are reachable from s->t.
            std::cout<<"ERROR: path not found"<<std::endl;
        }
        path_length[i]= size;
    }
    my_timer timer1 = my_timer();
    auto t_expand_begin = katch::util::time_stamp();
    // test k times, use global variable run_times!

    for(int runs = 0; runs < run_times; runs ++){
//////        cout<<"Finished: "<< runs <<"/"<<run_times<<endl;
        for(int i  = 0; i < source.size(); i ++) {

            int random_source = source[i];
            int random_target = target[i];

            query->load_target_row(departure_time);
            // This is only used for RTPD, because the RTPD are too large, so load one row and ignore the time
            // if you have more memory available for your RAM, you can change this by storing everything in RAM.

            timer1.start();
            double arr_t = query->one_to_one(random_source,random_target,departure_time);
            timer1.stop();
            double query_time = timer1.elapsed_time_micro();
            double path_retrieval_time = 0;

            path_recover_time[i] += path_retrieval_time;

            search_time[i] += query_time;
            if(best_time[i] > query_time){
                best_time[i]= query_time;
            }
            if(worst_time[i] < query_time){
                worst_time[i]= query_time;
            }
            if(time_cost[i] != 0){
                if(time_cost[i] != arr_t-departure_time){
                    std::cout<<"ERROR: path cost error"<<std::endl;
                }
            }
            time_cost[i] = arr_t-departure_time;
            number_of_nodes_expanded[i] = query->number_of_nodes_expanded;
            number_of_nodes_generated[i] = query->number_of_nodes_generated;
            number_of_first_move[i] = query->number_of_first_move_calls;
            number_of_TD_calculation[i] = query->number_of_TD_calculation;
            number_of_bw_first_move[i] = query-> number_of_bw_first_move_calls;
        }
    }
    auto t_expand_end = katch::util::time_stamp();
    std::cout<< katch::util::get_duration_in_seconds(t_expand_begin, t_expand_end)/run_times/source.size()*1000000<<std::endl;
    double total_search_time = 0;
    double total_path_recover_time = 0;
    double total_path_length = 0 ;
    double total_time_cost= 0 ;
    double total_generated= 0 ;
    double total_expanded= 0 ;
    double total_first_moves= 0;
    double total_TD_calculation= 0;
    double total_bw_first_moves = 0;
    for(int i  = 0; i < source.size(); i ++) {
        // take average;
        search_time[i] = (search_time[i]- best_time[i] - worst_time[i])/(run_times-2);
        path_recover_time[i] = path_recover_time[i]/run_times;
        total_path_recover_time += path_recover_time[i];
//        search_time[i] = (search_time[i]- best_time[i] - worst_time[i])/(run_times-2);
        total_search_time += search_time[i];
        total_path_length += path_length[i];
        total_time_cost += time_cost[i];
        total_generated += number_of_nodes_generated[i];
        total_expanded += number_of_nodes_expanded[i];
        total_first_moves += number_of_first_move[i];
        total_TD_calculation += number_of_TD_calculation[i];
        total_bw_first_moves +=   number_of_bw_first_move[i];
    }
    std::ofstream myFile(result_file);
    myFile<<"run_time,time_cost,euclidean_distance,path_length,nodes_generated,nodes_expanded,first_moves,number_of_TD_calculation,bw_first_moves\n";
    for(int i = 0; i < source.size(); i++){
        myFile<<std::fixed<<std::setprecision(8)<<search_time[i]<<","<<time_cost[i]
              <<","<<euclidean_distance[i]
              <<","<<path_length[i]
              <<","<<number_of_nodes_generated[i]
              <<","<<number_of_nodes_expanded[i]
              <<","<<number_of_first_move[i]
                <<","<<number_of_TD_calculation[i]
                <<","<<  number_of_bw_first_move[i]
              <<"\n";
    }
    myFile.close();
    std::cout<<std::fixed<<std::setprecision(8)<< "Average performance time: " << total_search_time/ source.size()  << std::endl;
    std::cout<<std::fixed<<std::setprecision(8)<< "Average path recover time: " << total_path_recover_time/ source.size()  << std::endl;
    std::cout<<std::fixed<<std::setprecision(8)<< "Average path length: " << total_path_length/ source.size()  << std::endl;
    std::cout<<std::fixed<<std::setprecision(8)<< "Average time cost: " << total_time_cost/ source.size()  << std::endl;
    std::cout<<std::fixed<<std::setprecision(8)<< "Average number of nodes generated: " << total_generated/ source.size()  << std::endl;
    std::cout<<std::fixed<<std::setprecision(8)<< "Average number of nodes expanded: " << total_expanded/ source.size()  << std::endl;
    std::cout<<std::fixed<<std::setprecision(8)<< "Average number of first moves: " << total_first_moves/ source.size()  << std::endl;
    std::cout<<std::fixed<<std::setprecision(8)<< "Average number of bw first moves " << total_bw_first_moves/ source.size()  << std::endl;
    std::cout<<std::fixed<<std::setprecision(8)<< "Average number of TD calculations: " << total_TD_calculation/ source.size()  << std::endl;
    KATCH_STATUS("Save experimental results to "<<result_file<<"\n");
}


void run_experiments( const std::string& algorithm, const std::string& input_dir_name,const std::string& output_dir_name){
    string dir_name = input_dir_name.substr(input_dir_name.find_last_of("/")+1);
    std::string main_name = input_dir_name+"/"+dir_name;
    source =  load_vector<int>(main_name+".source");
    target =  load_vector<int>(main_name+".target");
    euclidean_distance = load_vector<double>(main_name+".Euclidean");
    if(algorithm == "tch"){
        KATCH_STATUS("Benchmark algorithm: TCH\n");
        using Graph = katch::BTCH::Graph;
        using EdgeInfo = katch::EdgeInfo<Graph::TTF>;
        std::vector<EdgeInfo> edge_list = katch::btch_format::read_edges<EdgeInfo>(main_name+".btch",864000.0);
        Graph graph(std::move(edge_list));
        KATCH_STATUS("Graph has " << graph.get_n_nodes() << " nodes, " << graph.get_n_edges() << " edges.\n");
        auto query = new katch::BTCH(std::move(graph));

        for(int departure_time = 0; departure_time <= 864000; departure_time += test_period){

            string result_file = output_dir_name +"/tch_"+ to_string((int)departure_time)+".csv";
            run_search(query,result_file,departure_time);
        }
        std::cout<<std::endl;
    }else if(algorithm == "tch_l_star"){

        KATCH_STATUS("Benchmark algorithm: TCH with landmark heuristic\n");
        using Graph = katch::BTCH_L::Graph;
        using EdgeInfo = katch::EdgeInfo<Graph::TTF>;
        std::vector<EdgeInfo> edge_list = katch::btch_format::read_edges<EdgeInfo>(main_name+".btch",864000.0);
        Graph graph(std::move(edge_list));
        graph.load_coordinate(main_name+".coordinate");
        KATCH_STATUS("Graph has " << graph.get_n_nodes() << " nodes, " << graph.get_n_edges() << " edges.\n");
        auto query = new katch::BTCH_L(std::move(graph));
        std::vector<unsigned> landmarks_options = {12};
        for ( auto num_landmarks : landmarks_options){
            query->load_landmark(main_name+".landmark_"+ to_string(num_landmarks),num_landmarks);
            for(int departure_time = 0; departure_time <= 864000; departure_time += test_period){

                string result_file = output_dir_name +"/tch_l_star_" + to_string(num_landmarks)+"_"+ to_string((int)departure_time)+".csv";
                run_search(query,result_file,departure_time);
            }
            std::cout<<std::endl;
        }


    }else if (algorithm == "tch_cpd"){
        KATCH_STATUS("Benchmark algorithm: TCH CPD\n");
        using Graph = katch::FTCH_TCPD::Graph;
        using EdgeInfo = katch::EdgeInfo<Graph::TTF>;
        std::vector<EdgeInfo> edge_list = katch::btch_format::read_edges<EdgeInfo>(main_name+".btch",864000.0);
        Graph graph(std::move(edge_list));
        auto ordering = katch::btch_format::read_ordering(main_name+".btch",864000.0);
        graph.set_ranking(ordering);

        KATCH_STATUS("Graph has " << graph.get_n_nodes() << " nodes, " << graph.get_n_edges() << " edges.\n");
        auto query = new katch::FTCH_TCPD(std::move(graph));
        query->load_cpd(main_name+"_min.fw_tch_cpd");

//        int departure_time = 0;
        for(int departure_time = 0; departure_time <= 864000; departure_time += test_period){
            string result_file = output_dir_name +"/tch_cpd_"+ to_string((int)departure_time)+".csv";
            run_search(query,result_file,departure_time);
        }
        std::cout<<std::endl;
    }
    else if (algorithm == "rev_tch_cpd"){
        KATCH_STATUS("Benchmark algorithm: Rev TCH CPD\n");
        using Graph = katch::FTCH_RTPD::Graph;
        using EdgeInfo = katch::EdgeInfo<Graph::TTF>;
        std::vector<EdgeInfo> edge_list = katch::btch_format::read_edges<EdgeInfo>(main_name+".btch",864000.0);
        Graph graph(std::move(edge_list));
        auto ordering = katch::btch_format::read_ordering(main_name+".btch",864000.0);
        graph.set_ranking(ordering);

        KATCH_STATUS("Graph has " << graph.get_n_nodes() << " nodes, " << graph.get_n_edges() << " edges.\n");
        auto query = new katch::FTCH_RTPD(std::move(graph));
        query->load_cpd(main_name+"_min.rev_fw_tch_cpd");

//        int departure_time = 0;
        for(int departure_time = 0; departure_time <= 864000; departure_time += test_period){
            string result_file = output_dir_name +"/rev_tch_cpd_"+ to_string((int)departure_time)+".csv";
            run_search(query,result_file,departure_time);
        }
        std::cout<<std::endl;
    }else if (algorithm == "tch_cpd_l"){
        KATCH_STATUS("Benchmark algorithm: TCH CPD with Landmark heuristic\n");
        using Graph = katch::FTCH_L::Graph;
        using EdgeInfo = katch::EdgeInfo<Graph::TTF>;
        std::vector<EdgeInfo> edge_list = katch::btch_format::read_edges<EdgeInfo>(main_name+".btch",864000.0);
        Graph graph(std::move(edge_list));
        graph.load_coordinate(main_name+".coordinate");
//        graph.check_the_graph();
        KATCH_STATUS("Graph has " << graph.get_n_nodes() << " nodes, " << graph.get_n_edges() << " edges.\n");
        auto query = new katch::FTCH_L(std::move(graph));
        query->load_cpd(main_name+"_min.bw_tch_cpd_ch_order");

        std::vector<unsigned> landmarks_options = {12};
        for ( auto num_landmarks : landmarks_options){
            query->load_landmark(main_name+".landmark_"+ to_string(num_landmarks),num_landmarks);
            for(int departure_time = 0; departure_time <= 864000; departure_time += test_period){

                string result_file = output_dir_name +"/tch_cpd_l_" + to_string(num_landmarks)+"_"+ to_string((int)departure_time)+".csv";
                run_search(query,result_file,departure_time);
            }
            std::cout<<std::endl;
        }
    }
    else if (algorithm == "htch"){
        KATCH_STATUS("Benchmark algorithm: Hourly TCH\n");
        using Graph = katch::BSTCH::Graph;
        using EdgeInfo = katch::EdgeInfo<Graph::TTF>;
        std::vector<Graph> graph_vector;
        for(int i  =0 ;i < 24 ; i++){
            std::vector<EdgeInfo> edge_list = katch::btch_format::read_edges<EdgeInfo>(main_name+"_"+std::to_string(i)+".btch",hourly_TCH_period);
            graph_vector.push_back(Graph(std::move(edge_list)));
        }
        graph_vector.shrink_to_fit();
        auto query = new katch::BSTCH(std::move(graph_vector));
        for(int departure_time = 0; departure_time <= 864000; departure_time += test_period){
            string result_file = output_dir_name +"/htch_"+ to_string((int)departure_time)+".csv";
            run_search(query,result_file,departure_time);
        }
        std::cout<<std::endl;
    }    else if(algorithm == "htch_cpd"){
        using Graph = katch::FSTCH_TCPD::Graph;
        using EdgeInfo = katch::EdgeInfo<Graph::TTF>;
        std::vector<Graph> graph_vector;
        for(int i  =0 ;i < 24 ; i++){
            std::vector<EdgeInfo> edge_list = katch::btch_format::read_edges<EdgeInfo>(main_name+"_"+std::to_string(i)+".btch",hourly_TCH_period);
            graph_vector.push_back(Graph(std::move(edge_list)));
        }
        graph_vector.shrink_to_fit();
        auto query = new katch::FSTCH_TCPD(std::move(graph_vector));
        query->load_cpd(main_name+"_min.fw_htch_cpd");
        for(int departure_time = 0; departure_time <= 864000; departure_time += test_period){
            string result_file = output_dir_name +"/htch_cpd_"+ to_string((int)departure_time)+".csv";
            run_search(query,result_file,departure_time);
        }
        std::cout<<std::endl;

    }
    else if (algorithm == "ts_htch"){
        KATCH_STATUS("Benchmark algorithm: Hourly TCH\n");
        using Graph = katch::BMTCHQuery<& one_TCH_period >::Graph;
        using EdgeInfo = katch::EdgeInfo<Graph::TTF>;
        std::vector<Graph> graph_vector;
        for(int i  =0 ;i < 24 ; i+=1){
            std::vector<EdgeInfo> edge_list = katch::btch_format::read_edges<EdgeInfo>(main_name+"_1_1_"+std::to_string(i)+".btch", one_TCH_period );
            graph_vector.push_back(Graph(std::move(edge_list)));
        }
        graph_vector.shrink_to_fit();
        auto hour_query = new katch::BMTCHQuery<&one_TCH_period > (std::move(graph_vector));

        using Graph2 = katch::BMTCHQuery<& two_TCH_period >::Graph;
        using EdgeInfo2 = katch::EdgeInfo<Graph2::TTF>;
        std::vector<Graph2> graph_vector2;
        for(int i  =0 ;i < 24 ; i+=1){
            std::vector<EdgeInfo2> edge_list2 = katch::btch_format::read_edges<EdgeInfo2>(main_name+"_1_2_"+std::to_string(i)+".btch", two_TCH_period );
            graph_vector2.push_back(Graph2(std::move(edge_list2)));
        }
        graph_vector2.shrink_to_fit();
        auto hour_query2 = new katch::BMTCHQuery<&two_TCH_period > (std::move(graph_vector2));

//
        using Graph4 = katch::BMTCHQuery<& four_TCH_period >::Graph;
        using EdgeInfo4 = katch::EdgeInfo<Graph4::TTF>;
        std::vector<Graph4> graph_vector4;
        for(int i  =0 ;i < 24 ; i+=1){
            std::vector<EdgeInfo4> edge_list4 = katch::btch_format::read_edges<EdgeInfo4>(main_name+"_1_4_"+std::to_string(i)+".btch", four_TCH_period );
            graph_vector4.push_back(Graph4(std::move(edge_list4)));
        }
        graph_vector4.shrink_to_fit();
        auto hour_query4 = new katch::BMTCHQuery<&four_TCH_period > (std::move(graph_vector4));
        std::vector<katch::BMTCHQuery_Wrapper*> hour_query_vector;
        hour_query_vector.push_back(hour_query);
        hour_query_vector.push_back(hour_query2);
        hour_query_vector.push_back(hour_query4);

        using ts_tch_cpd_Graph = katch::BMTCH::Graph;
        using ts_tch_cpd_EdgeInfo = katch::EdgeInfo<ts_tch_cpd_Graph::TTF>;
        std::vector<ts_tch_cpd_EdgeInfo> ts_tch_cpd_edge_list = katch::btch_format::read_edges<ts_tch_cpd_EdgeInfo>(main_name+".btch",864000.0);
        ts_tch_cpd_Graph ts_tch_cpd_graph(std::move(ts_tch_cpd_edge_list ));
        auto ordering = katch::btch_format::read_ordering(main_name+".btch",864000.0);
        ts_tch_cpd_graph.set_ranking(ordering);

        KATCH_STATUS("Graph has " << ts_tch_cpd_graph.get_n_nodes() << " nodes, " << ts_tch_cpd_graph.get_n_edges() << " edges.\n");
        auto query = new katch::BMTCH(std::move(ts_tch_cpd_graph));
        query->load_cpd(main_name+"_min.fw_tch_cpd");
        query->set_hourly_queries( hour_query_vector);


        for(int departure_time = 0; departure_time <= 864000; departure_time += test_period){
            string result_file = output_dir_name +"/ts_htch_"+ to_string((int)departure_time)+".csv";
            query->query_answered_by_1_hour_query = 0;
            query->query_answered_by_2_hour_query = 0;
            query->query_answered_by_4_hour_query = 0;
            query->query_answered_by_tch_query = 0;
            run_search(query,result_file,departure_time);
            std::cout<<"Average query answer by 1 hour " << query->query_answered_by_1_hour_query / run_times << std::endl;
            std::cout<<"Average query answer by 2 hour " << query->query_answered_by_2_hour_query / run_times << std::endl;
            std::cout<<"Average query answer by 4 hour " << query->query_answered_by_4_hour_query / run_times << std::endl;
            std::cout<<"Average query answer by TCH " << query->query_answered_by_tch_query / run_times  << std::endl;

        }
        std::cout<<std::endl;
    } else if (algorithm == "ts_htch_cpd") {
        KATCH_STATUS("Benchmark algorithm: Hourly TCH\n");
        using Graph = katch::FMTCH_TCPDQuery<& one_TCH_period >::Graph;
        using EdgeInfo = katch::EdgeInfo<Graph::TTF>;
        std::vector<Graph> graph_vector;
        for(int i  =0 ;i < 24 ; i+=1){
            std::vector<EdgeInfo> edge_list = katch::btch_format::read_edges<EdgeInfo>(main_name+"_1_1_"+std::to_string(i)+".btch", one_TCH_period );
            graph_vector.push_back(Graph(std::move(edge_list)));
        }
        graph_vector.shrink_to_fit();
        auto hour_query = new katch::FMTCH_TCPDQuery<&one_TCH_period > (std::move(graph_vector));
        hour_query->load_cpd(main_name, 1);

        using Graph2 = katch::FMTCH_TCPDQuery<& two_TCH_period >::Graph;
        using EdgeInfo2 = katch::EdgeInfo<Graph2::TTF>;
        std::vector<Graph2> graph_vector2;
        for(int i  =0 ;i < 24 ; i+=1){
            std::vector<EdgeInfo2> edge_list2 = katch::btch_format::read_edges<EdgeInfo2>(main_name+"_1_2_"+std::to_string(i)+".btch", two_TCH_period );
            graph_vector2.push_back(Graph2(std::move(edge_list2)));
        }
        graph_vector2.shrink_to_fit();
        auto hour_query2 = new katch::FMTCH_TCPDQuery<&two_TCH_period > (std::move(graph_vector2));
        hour_query2->load_cpd(main_name, 2);
//
        using Graph4 = katch::FMTCH_TCPDQuery<& four_TCH_period >::Graph;
        using EdgeInfo4 = katch::EdgeInfo<Graph4::TTF>;
        std::vector<Graph4> graph_vector4;
        for(int i  =0 ;i < 24 ; i+=1){
            std::vector<EdgeInfo4> edge_list4 = katch::btch_format::read_edges<EdgeInfo4>(main_name+"_1_4_"+std::to_string(i)+".btch", four_TCH_period );
            graph_vector4.push_back(Graph4(std::move(edge_list4)));
        }
        graph_vector4.shrink_to_fit();
        auto hour_query4 = new katch::FMTCH_TCPDQuery<&four_TCH_period> (std::move(graph_vector4));
        hour_query4->load_cpd(main_name, 4);

        std::vector<katch::FMTCH_TCPDQuery_Wrapper*> hour_query_vector;
        hour_query_vector.push_back(hour_query);
        hour_query_vector.push_back(hour_query2);
        hour_query_vector.push_back(hour_query4);
        hour_query_vector.shrink_to_fit();
        using ts_tch_cpd_Graph = katch::FMTCH_TCPD::Graph;
        using ts_tch_cpd_EdgeInfo = katch::EdgeInfo<ts_tch_cpd_Graph::TTF>;
        std::vector<ts_tch_cpd_EdgeInfo> ts_tch_cpd_edge_list = katch::btch_format::read_edges<ts_tch_cpd_EdgeInfo>(main_name+".btch",864000.0);
        ts_tch_cpd_Graph ts_tch_cpd_graph(std::move(ts_tch_cpd_edge_list ));
        auto ordering = katch::btch_format::read_ordering(main_name+".btch",864000.0);
        ts_tch_cpd_graph.set_ranking(ordering);

        KATCH_STATUS("Graph has " << ts_tch_cpd_graph.get_n_nodes() << " nodes, " << ts_tch_cpd_graph.get_n_edges() << " edges.\n");
        auto query = new katch::FMTCH_TCPD(std::move(ts_tch_cpd_graph));
        query->load_cpd(main_name+"_min.fw_tch_cpd");
        query->set_hourly_queries( hour_query_vector);


        for(int departure_time = 0; departure_time <= 864000; departure_time += test_period){
            string result_file = output_dir_name +"/ts_htch_cpd_"+ to_string((int)departure_time)+".csv";
            query->query_answered_by_1_hour_query = 0;
            query->query_answered_by_2_hour_query = 0;
            query->query_answered_by_4_hour_query = 0;
            query->query_answered_by_tch_query = 0;
            run_search(query,result_file,departure_time);
            std::cout<<"Average query answer by 1 hour " << query->query_answered_by_1_hour_query / run_times << std::endl;
            std::cout<<"Average query answer by 2 hour " << query->query_answered_by_2_hour_query / run_times << std::endl;
            std::cout<<"Average query answer by 4 hour " << query->query_answered_by_4_hour_query / run_times << std::endl;
            std::cout<<"Average query answer by TCH " << query->query_answered_by_tch_query / run_times  << std::endl;

        }
        std::cout<<std::endl;


    }
    else if (algorithm == "tch_ori_cpd"){
        KATCH_STATUS("Benchmark algorithm: TCH with original CPD\n");
        using Graph = katch::FTCH_CPD::Graph;
        using EdgeInfo = katch::EdgeInfo<Graph::TTF>;
        std::vector<EdgeInfo> edge_list = katch::btch_format::read_edges<EdgeInfo>(main_name+".btch",864000.0);
        Graph graph(std::move(edge_list));
        auto ordering = katch::btch_format::read_ordering(main_name+".btch",864000.0);
        graph.set_ranking(ordering);

        KATCH_STATUS("Graph has " << graph.get_n_nodes() << " nodes, " << graph.get_n_edges() << " edges.\n");
        auto query = new katch::FTCH_CPD(std::move(graph));
        query->load_cpd(main_name+"_min.ori_fw_tch_cpd");

//        int departure_time = 0;
        for(int departure_time = 0; departure_time <= 864000; departure_time += test_period){
            string result_file = output_dir_name +"/tch_ori_cpd_"+ to_string((int)departure_time)+".csv";
            run_search(query,result_file,departure_time);
        }
        std::cout<<std::endl;

    }else {
        KATCH_STATUS("Can not find the algorithm\n");
    }
}

int main(int argc, char** argv) {

    const char *binary_name = argv[0];

    if (argc != 4) {
        std::cerr
                << std::endl
                << "USAGE: " << binary_name
                << "<algorithm> <input dir name> <output dir name> <departure time>" << std::endl
                << std::endl;

        return EXIT_FAILURE;
    }
    const std::string algorithm_name(argv[1]);
    const std::string input_dir_name(argv[2]);
    const std::string output_dir_name(argv[3]);
//    const double departure_time(std::stod(argv[4]));

    run_experiments(algorithm_name,input_dir_name,output_dir_name);

}