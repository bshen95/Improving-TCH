//
// Created by Bojie Shen on 29/5/21.
//

#include <iostream>
#include <stdlib.h>
#include <string>
#include <util/my_timer.h>
#include <iomanip>

#include "util/misc.h"
#include "io/btch_format.h"
#include "io/edge_info.h"
#include "io/vec_io.h"
#include "datastr/graph/forward_search_graph.h"
#include "cpd_prepro/tch_downward_dijkstra_query.h"
#include "omp.h"
#include "datastr/graph/dynamic_search_graph.h"
#include "io/tpgr_format.h"
#include "datastr/cpd/bw_cpd.h"
static constexpr double period = 864000.0;
static constexpr double h_period_4 =  149000.0;
static constexpr double STCH_period = 185000.0;
static constexpr double offset = 5000;
static constexpr double hourly_TCH_period = 36000.0 * 5 + offset;
static constexpr double one_MTCH_period = 36000.0 *2 + offset;
static constexpr double two_MTCH_period = 36000.0 *3 + offset;
static constexpr double four_MTCH_period = 36000.0 *5 + offset;

using S_Graph = katch::TchDownardDijkstraQuery::Graph;
using F_Graph = katch::ForwardSearchGraph<&period>;
using EdgeInfo = katch::EdgeInfo<katch::ForwardSearchGraph<&period>::TTF>;

using D_Graph = katch::DynamicSearchGraph;
using DEdgeInfo = katch::EdgeInfo<katch::DynamicSearchGraph::TTF>;

vector<double> construction_time;
vector<double> memory_cost;

//// using double queue, usually slower than int queue, this is confirmed.
//void construct_min_bw_tch_cpd_with_ch_ordering (const std::string& tpgr_input_file_name, const std::string& btch_input_file_name, const string& output_dir_name) {
//    string main_name = btch_input_file_name.substr(0,btch_input_file_name.find_last_of("."));
//
//    // convert forward graph to static graph.
//    std::vector<EdgeInfo> edge_list = katch::btch_format::read_edges<EdgeInfo>(btch_input_file_name,period);
//    F_Graph f_graph(std::move(edge_list));
//    KATCH_STATUS("Graph has " << f_graph.get_n_nodes() << " nodes, " << f_graph.get_n_edges() << " edges.\n");
//    S_Graph s_graph(f_graph,F_Graph::min_graph);
//    auto ordering = katch::btch_format::read_ordering(btch_input_file_name,period);
//    f_graph.test_ordering(ordering);
//    // output the free_flow cost
//    std::vector<double> free_flow = s_graph.get_free_flow_vector();
//    save_vector(main_name+"_min.bw_tch_free_flow_ch_order",free_flow);
//    KATCH_STATUS("Save free flow to " << main_name << "_min.bw_tch_free_flow_ch_order"<<"\n");
//
//
//
//    //generate dfs ordering using the original graph.
//
//
//    my_timer timer1 = my_timer();
//    timer1.start();
//    // resort the graph based on dfs ordering.
//    std::vector<katch::NodeIterator > dfs_ordering = f_graph.generate_ch_based_dfs_ordering(ordering);
//    s_graph.resort_graph(dfs_ordering);
//
//    KATCH_STATUS("Building Min CPD ..." <<"\n");
//    std::vector<katch::TchDownardDijkstraQuery*> dijkstra;
//    for (int i = 0; i < omp_get_max_threads(); i ++){
//        dijkstra.push_back(new katch::TchDownardDijkstraQuery(&s_graph));
//    }
//    unsigned number_of_nodes = s_graph.get_n_nodes();
//    auto* cpd = new katch::CPD();
//    {
//        {
//            katch::TchDownardDijkstraQuery dij(&s_graph);
//            my_timer t;
//            t.start();
//            dij.get_first_move(0);
//            t.stop();
//            double tots = t.elapsed_time_micro()*number_of_nodes / 1000000;
//            printf("Estimated sequential running time : %fmin\n", tots / 60.0);
//
//        }
//
//        printf("Using %d threads\n", omp_get_max_threads());
//        std::vector<katch::CPD>thread_cpd(omp_get_max_threads());
//
//        int progress = 0;
//
//#pragma omp parallel
//        {
//            const int thread_count = omp_get_num_threads();
//            const int thread_id = omp_get_thread_num();
////            const int begin_int  = 0;
//            const int node_count = number_of_nodes;
//
//            int node_begin = (node_count*thread_id) / thread_count ;
//            int node_end = (node_count*(thread_id+1)) / thread_count ;
//            katch::TchDownardDijkstraQuery& thread_dij = *dijkstra[thread_id];
//
//            for(int source_node=node_begin; source_node < node_end; ++source_node){
//                const std::vector<unsigned short>& result = thread_dij.get_first_move(source_node);
//                thread_cpd[thread_id].append_row(source_node,result);
//#pragma omp critical
//                {
//                    ++progress;
//                    if(progress % 100 == 0) {
//                        double ratio = (double)progress / number_of_nodes * 100.0;
//                        std::cout << "Progress: [" << progress << "/" << number_of_nodes << "] "
//                                  << std::setprecision(3) << ratio << "% \r";
//                        std::cout.flush();
//                    }
//                }
//            }
//        }
//
//        for(auto&x:thread_cpd)
//            cpd->append_rows(x);
//    }
//    timer1.stop();
//    for (auto i: dijkstra) delete i;
//
//    printf("Saving data to cpd.txt \n");
//    printf("begin size: %d, entry size: %d\n", cpd->entry_count(), cpd->get_entry_size());
//
//
//    string mapper_file = main_name +".bw_tch_mapper_ch_order";
//    vector<katch::NodeIterator > mapper = invert_permutation(dfs_ordering);
//    save_vector(mapper_file,mapper);
//    KATCH_STATUS("Save mapper to " <<mapper_file<<"\n");
//
//    std::string fname = main_name+"_min.bw_tch_cpd_ch_order";
//    FILE*f = fopen(fname.c_str(), "wb");
//    cpd->save(f);
//    fclose(f);
//    cout << "done" << endl;
//    KATCH_STATUS("Save CPD to " <<fname<<"\n");
//
//    std::ofstream myFile(output_dir_name+"/bw_tch_cpd_preprocessing_min_ch_order.csv");
//    myFile<<"build_time,memory_cost\n";
//    myFile<<std::fixed<<setprecision(8)<<timer1.elapsed_time_mins()<<","<<(double)cpd->get_entry_size()*4/1000000
//          <<"\n";
//    myFile.close();
//
//    std::cout<< "Finished in: " << timer1.elapsed_time_mins()<<" Mins"<< std::endl;
//    std::cout<< "Memory cost: " << (double)cpd->get_entry_size()*4/1000000<<" MB"<< std::endl;
//}


void construct_min_bw_tcpd_with_ch_ordering (const std::string& tpgr_input_file_name, const std::string& btch_input_file_name, const string& output_dir_name) {
    string main_name = btch_input_file_name.substr(0,btch_input_file_name.find_last_of("."));

    // convert forward graph to static graph.
    std::vector<EdgeInfo> edge_list = katch::btch_format::read_edges<EdgeInfo>(btch_input_file_name,period);
    F_Graph f_graph(std::move(edge_list));
    KATCH_STATUS("Graph has " << f_graph.get_n_nodes() << " nodes, " << f_graph.get_n_edges() << " edges.\n");
    S_Graph s_graph(f_graph,F_Graph::min_graph);
    auto ordering = katch::btch_format::read_ordering(btch_input_file_name,period);
    f_graph.test_ordering(ordering);
    std::vector<bool> lowest_level_nodes =  f_graph.get_lowest_level_mapper();
    // output the free_flow cost
    std::vector<double> free_flow = s_graph.get_free_flow_vector();
    save_vector(main_name+"_min.bw_tch_free_flow_ch_order",free_flow);
    KATCH_STATUS("Save free flow to " << main_name << "_min.bw_tch_free_flow_ch_order"<<"\n");



    //generate dfs ordering using the original graph.


    my_timer timer1 = my_timer();
    timer1.start();
    // resort the graph based on dfs ordering.
    std::vector<katch::NodeIterator > dfs_ordering = f_graph.generate_ch_based_dfs_ordering(ordering);
    s_graph.resort_graph(dfs_ordering);
    // dfs_ordering: cpd_node -> v_id;

    string mapper_file = main_name +".bw_tch_mapper_ch_order";
    vector<katch::NodeIterator > mapper = invert_permutation(dfs_ordering);
    save_vector(mapper_file,mapper);
    KATCH_STATUS("Save mapper to " <<mapper_file<<"\n");


    KATCH_STATUS("Building Min CPD ..." <<"\n");
    std::vector<katch::TchDownardDijkstraQuery*> dijkstra;
    for (int i = 0; i < omp_get_max_threads(); i ++){
        dijkstra.push_back(new katch::TchDownardDijkstraQuery(&s_graph));
    }
    unsigned number_of_nodes = s_graph.get_n_nodes();
    auto* bw_cpd = new katch::BW_CPD();
    {
        {
            katch::TchDownardDijkstraQuery dij(&s_graph);
            my_timer t;
            t.start();
            dij.get_first_move(0);
            t.stop();
            double tots = t.elapsed_time_micro()*number_of_nodes / 1000000;
            printf("Estimated sequential running time : %fmin\n", tots / 60.0);

        }

        printf("Using %d threads\n", omp_get_max_threads());
        std::vector<katch::BW_CPD>thread_cpd(omp_get_max_threads());

        int progress = 0;

#pragma omp parallel
        {
            const int thread_count = omp_get_num_threads();
            const int thread_id = omp_get_thread_num();
//            const int begin_int  = 0;
            const int node_count = number_of_nodes;

            int node_begin = (node_count*thread_id) / thread_count ;
            int node_end = (node_count*(thread_id+1)) / thread_count ;
            katch::TchDownardDijkstraQuery& thread_dij = *dijkstra[thread_id];

            for(int source_node=node_begin; source_node < node_end; ++source_node){
//                const std::vector<unsigned short>& result = thread_dij.get_reachability_with_wild_card(mapper[source_node]);
//                thread_cpd[thread_id].append_row(mapper[source_node],result);
                const std::vector<bool>& result = thread_dij.get_reachability(mapper[source_node]);
                thread_cpd[thread_id].append_row(mapper[source_node],result);
#pragma omp critical
                {
                    ++progress;
                    if(progress % 100 == 0) {
                        double ratio = (double)progress / number_of_nodes * 100.0;
                        std::cout << "Progress: [" << progress << "/" << number_of_nodes << "] "
                                  << std::setprecision(3) << ratio << "% \r";
                        std::cout.flush();
                    }
                }
            }
        }

        for(auto&x:thread_cpd)
            bw_cpd->append_rows(x);
    }
    timer1.stop();
    for (auto i: dijkstra) delete i;

    printf("Saving data to cpd.txt \n");
    printf("begin size: %d, entry size: %d\n", bw_cpd->entry_count(), bw_cpd->get_entry_size());



    std::string fname = main_name+"_min.bw_tch_cpd_ch_order";
    FILE*f = fopen(fname.c_str(), "wb");
    bw_cpd->save(f);
    fclose(f);
    cout << "done" << endl;
    KATCH_STATUS("Save CPD to " <<fname<<"\n");

    std::ofstream myFile(output_dir_name+"/bw_tch_cpd_preprocessing_min_ch_order.csv");
    myFile<<"build_time,memory_cost\n";
    myFile<<std::fixed<<setprecision(8)<<timer1.elapsed_time_mins()<<","<<(double)bw_cpd->get_entry_size()*4/1000000
          <<"\n";
    myFile.close();

    std::cout<< "Finished in: " << timer1.elapsed_time_mins()<<" Mins"<< std::endl;
    std::cout<< "Memory cost: " << (double)bw_cpd->get_entry_size()*4/1000000<<" MB"<< std::endl;
}

void construct_min_bw_tch_cpd(const std::string& tpgr_input_file_name, const std::string& btch_input_file_name, const string& output_dir_name) {
    string main_name = btch_input_file_name.substr(0,btch_input_file_name.find_last_of("."));

    // convert forward graph to static graph.
    std::vector<EdgeInfo> edge_list = katch::btch_format::read_edges<EdgeInfo>(btch_input_file_name,period);
    F_Graph f_graph(std::move(edge_list));
    KATCH_STATUS("Graph has " << f_graph.get_n_nodes() << " nodes, " << f_graph.get_n_edges() << " edges.\n");
    S_Graph s_graph(f_graph,F_Graph::min_graph);
    // output the free_flow cost
    std::vector<double> free_flow = s_graph.get_free_flow_vector();
    save_vector(main_name+"_min.bw_tch_free_flow",free_flow);
    KATCH_STATUS("Save free flow to " << main_name << "_min.bw_tch_free_flow"<<"\n");



    //generate dfs ordering using the original graph.
    std::vector<DEdgeInfo> d_edge_list = katch::tpgr_format::read_edges<DEdgeInfo>(tpgr_input_file_name);
    D_Graph d_graph(std::move(d_edge_list));

    my_timer timer1 = my_timer();
    timer1.start();
    // resort the graph based on dfs ordering.
    std::vector<katch::NodeIterator > dfs_ordering = d_graph.generate_DFS_ordering();
    s_graph.resort_graph(dfs_ordering);

    KATCH_STATUS("Building Min CPD ..." <<"\n");
    std::vector<katch::TchDownardDijkstraQuery*> dijkstra;
    for (int i = 0; i < omp_get_max_threads(); i ++){
        dijkstra.push_back(new katch::TchDownardDijkstraQuery(&s_graph));
    }
    unsigned number_of_nodes = s_graph.get_n_nodes();
    auto* cpd = new katch::BW_CPD();
    {
        {
            katch::TchDownardDijkstraQuery dij(&s_graph);
            my_timer t;
            t.start();
            dij.get_first_move(0);
            t.stop();
            double tots = t.elapsed_time_micro()*number_of_nodes / 1000000;
            printf("Estimated sequential running time : %fmin\n", tots / 60.0);

        }

        printf("Using %d threads\n", omp_get_max_threads());
        std::vector<katch::BW_CPD>thread_cpd(omp_get_max_threads());

        int progress = 0;

#pragma omp parallel
        {
            const int thread_count = omp_get_num_threads();
            const int thread_id = omp_get_thread_num();
//            const int begin_int  = 0;
            const int node_count = number_of_nodes;

            int node_begin = (node_count*thread_id) / thread_count ;
            int node_end = (node_count*(thread_id+1)) / thread_count ;
            katch::TchDownardDijkstraQuery& thread_dij = *dijkstra[thread_id];

            for(int source_node=node_begin; source_node < node_end; ++source_node){
                const std::vector<unsigned short>& result = thread_dij.get_first_move(source_node);
                thread_cpd[thread_id].append_row(source_node,result);
#pragma omp critical
                {
                    ++progress;
                    if(progress % 100 == 0) {
                        double ratio = (double)progress / number_of_nodes * 100.0;
                        std::cout << "Progress: [" << progress << "/" << number_of_nodes << "] "
                                  << std::setprecision(3) << ratio << "% \r";
                        std::cout.flush();
                    }
                }
            }
        }

        for(auto&x:thread_cpd)
            cpd->append_rows(x);
    }
    timer1.stop();
    for (auto i: dijkstra) delete i;

    printf("Saving data to cpd.txt \n");
    printf("begin size: %d, entry size: %d\n", cpd->entry_count(), cpd->get_entry_size());


    string mapper_file = main_name +".bw_tch_mapper";
    vector<katch::NodeIterator > mapper = invert_permutation(dfs_ordering);
    save_vector(mapper_file,mapper);
    KATCH_STATUS("Save mapper to " <<mapper_file<<"\n");

    std::string fname = main_name+"_min.bw_tch_cpd";
    FILE*f = fopen(fname.c_str(), "wb");
    cpd->save(f);
    fclose(f);
    cout << "done" << endl;
    KATCH_STATUS("Save CPD to " <<fname<<"\n");

    std::ofstream myFile(output_dir_name+"/bw_tch_cpd_preprocessing_min.csv");
    myFile<<"build_time,memory_cost\n";
    myFile<<std::fixed<<setprecision(8)<<timer1.elapsed_time_mins()<<","<<(double)cpd->get_entry_size()*4/1000000
          <<"\n";
    myFile.close();

    std::cout<< "Finished in: " << timer1.elapsed_time_mins()<<" Mins"<< std::endl;
    std::cout<< "Memory cost: " << (double)cpd->get_entry_size()*4/1000000<<" MB"<< std::endl;
}

template<const double* h_period>
void construct_hourly_min_htch_cpd(const std::string& btch_input_file_name, const string& output_dir_name,
        const std::vector<katch::NodeIterator>& dfs_ordering, int hour) {
    using H_Graph = katch::ForwardSearchGraph<h_period>;
    using HEdgeInfo = katch::EdgeInfo<typename H_Graph::TTF>;
    string main_name = btch_input_file_name.substr(0,btch_input_file_name.find_last_of("."));
    std::vector<HEdgeInfo> edge_list = katch::btch_format::read_edges<HEdgeInfo>(btch_input_file_name,*h_period);
    H_Graph h_graph(std::move(edge_list));
    KATCH_STATUS("Graph has " << h_graph.get_n_nodes() << " nodes, " << h_graph.get_n_edges() << " edges.\n");
    S_Graph s_graph(h_graph, H_Graph::min_graph);

    my_timer timer1 = my_timer();
    timer1.start();
    // resort the graph based on dfs ordering.
    s_graph.resort_graph(dfs_ordering);
    KATCH_STATUS("Building Min CPD ..." <<"\n");
    std::vector<katch::TchDownardDijkstraQuery*> dijkstra;
    for (int i = 0; i < omp_get_max_threads(); i ++){
        dijkstra.push_back(new katch::TchDownardDijkstraQuery(&s_graph));
    }
    unsigned number_of_nodes = s_graph.get_n_nodes();
    auto* cpd = new katch::BW_CPD();
    {
        {
            katch::TchDownardDijkstraQuery dij(&s_graph);
            my_timer t;
            t.start();
            dij.get_first_move(0);
            t.stop();
            double tots = t.elapsed_time_micro()*number_of_nodes / 1000000;
            printf("Estimated sequential running time : %fmin\n", tots / 60.0);

        }

        printf("Using %d threads\n", omp_get_max_threads());
        std::vector<katch::BW_CPD>thread_cpd(omp_get_max_threads());

        int progress = 0;

#pragma omp parallel
        {
            const int thread_count = omp_get_num_threads();
            const int thread_id = omp_get_thread_num();
//            const int begin_int  = 0;
            const int node_count = number_of_nodes;

            int node_begin = (node_count*thread_id) / thread_count ;
            int node_end = (node_count*(thread_id+1)) / thread_count ;
            katch::TchDownardDijkstraQuery& thread_dij = *dijkstra[thread_id];

            for(int source_node=node_begin; source_node < node_end; ++source_node){
                const std::vector<unsigned short>& result = thread_dij.get_first_move(source_node);
                thread_cpd[thread_id].append_row(source_node,result);
#pragma omp critical
                {
                    ++progress;
                    if(progress % 100 == 0) {
                        double ratio = (double)progress / number_of_nodes * 100.0;
                        std::cout << "Progress: [" << progress << "/" << number_of_nodes << "] "
                                  << std::setprecision(3) << ratio << "% \r";
                        std::cout.flush();
                    }
                }
            }
        }

        for(auto&x:thread_cpd)
            cpd->append_rows(x);
    }
    timer1.stop();
    for (auto i: dijkstra) delete i;

    printf("Saving data to cpd.txt \n");
    printf("begin size: %d, entry size: %d\n", cpd->entry_count(), cpd->get_entry_size());


    std::string fname = main_name+"_min.bw_htch_cpd";
    FILE*f = fopen(fname.c_str(), "wb");
    cpd->save(f);
    fclose(f);
    cout << "done" << endl;
    KATCH_STATUS("Save CPD to " <<fname<<"\n");

    std::ofstream myFile(output_dir_name+"/bw_htch_cpd_preprocessing_min_"+to_string(hour)+".csv");
    myFile<<"build_time,memory_cost\n";
    myFile<<std::fixed<<setprecision(8)<<timer1.elapsed_time_mins()<<","<<(double)cpd->get_entry_size()*4/1000000
          <<"\n";
    myFile.close();

    std::cout<< "Finished in: " << timer1.elapsed_time_mins()<<" Mins"<< std::endl;
    std::cout<< "Memory cost: " << (double)cpd->get_entry_size()*4/1000000<<" MB"<< std::endl;
}


template<const double* h_period>
void construct_min_stch_bw_tcpd(const std::string& btch_input_file_name, const string& output_dir_name, int hour) {
    using H_Graph = katch::ForwardSearchGraph<h_period>;
    using HEdgeInfo = katch::EdgeInfo<typename H_Graph::TTF>;
    string main_name = btch_input_file_name.substr(0,btch_input_file_name.find_last_of("."));
    std::vector<HEdgeInfo> edge_list = katch::btch_format::read_edges<HEdgeInfo>(btch_input_file_name,*h_period);
    H_Graph h_graph(std::move(edge_list));
    KATCH_STATUS("Graph has " << h_graph.get_n_nodes() << " nodes, " << h_graph.get_n_edges() << " edges.\n");
    S_Graph s_graph(h_graph, H_Graph::min_graph);


    auto ordering = katch::btch_format::read_ordering(btch_input_file_name,*h_period);
    h_graph.test_ordering(ordering);



    my_timer timer1 = my_timer();
    timer1.start();
    // resort the graph based on dfs ordering.
    std::vector<katch::NodeIterator > dfs_ordering = h_graph.generate_ch_based_dfs_ordering(ordering);
    s_graph.resort_graph(dfs_ordering);

    string mapper_file = main_name +"_min.bw_htch_mapper";
    vector<katch::NodeIterator > mapper = invert_permutation(dfs_ordering);
    save_vector(mapper_file,mapper);

    KATCH_STATUS("Building Min CPD ..." <<"\n");
    std::vector<katch::TchDownardDijkstraQuery*> dijkstra;
    for (int i = 0; i < omp_get_max_threads(); i ++){
        dijkstra.push_back(new katch::TchDownardDijkstraQuery(&s_graph));
    }
    unsigned number_of_nodes = s_graph.get_n_nodes();
    auto* cpd = new katch::BW_CPD();
    {
        {
            katch::TchDownardDijkstraQuery dij(&s_graph);
            my_timer t;
            t.start();
            dij.get_first_move(0);
            t.stop();
            double tots = t.elapsed_time_micro()*number_of_nodes / 1000000;
            printf("Estimated sequential running time : %fmin\n", tots / 60.0);

        }

        printf("Using %d threads\n", omp_get_max_threads());
        std::vector<katch::BW_CPD>thread_cpd(omp_get_max_threads());

        int progress = 0;

#pragma omp parallel
        {
            const int thread_count = omp_get_num_threads();
            const int thread_id = omp_get_thread_num();
//            const int begin_int  = 0;
            const int node_count = number_of_nodes;

            int node_begin = (node_count*thread_id) / thread_count ;
            int node_end = (node_count*(thread_id+1)) / thread_count ;
            katch::TchDownardDijkstraQuery& thread_dij = *dijkstra[thread_id];

            for(int source_node=node_begin; source_node < node_end; ++source_node){
                const std::vector<bool>& result = thread_dij.get_reachability( mapper[source_node]);
                thread_cpd[thread_id].append_row( mapper[source_node],result);
#pragma omp critical
                {
                    ++progress;
                    if(progress % 100 == 0) {
                        double ratio = (double)progress / number_of_nodes * 100.0;
                        std::cout << "Progress: [" << progress << "/" << number_of_nodes << "] "
                                  << std::setprecision(3) << ratio << "% \r";
                        std::cout.flush();
                    }
                }
            }
        }

        for(auto&x:thread_cpd)
            cpd->append_rows(x);
    }
    timer1.stop();
    for (auto i: dijkstra) delete i;

    printf("Saving data to cpd.txt \n");
    printf("begin size: %d, entry size: %d\n", cpd->entry_count(), cpd->get_entry_size());


    std::string fname = main_name+"_min.bw_htch_cpd";
    FILE*f = fopen(fname.c_str(), "wb");
    cpd->save(f);
    fclose(f);
    cout << "done" << endl;
    KATCH_STATUS("Save CPD to " <<fname<<"\n");

    std::ofstream myFile(output_dir_name+"/bw_htch_cpd_preprocessing_min_"+to_string(hour)+".csv");
    myFile<<"build_time,memory_cost\n";
    myFile<<std::fixed<<setprecision(8)<<timer1.elapsed_time_mins()<<","<<(double)cpd->get_entry_size()*4/1000000
          <<"\n";
    myFile.close();

    std::cout<< "Finished in: " << timer1.elapsed_time_mins()<<" Mins"<< std::endl;
    std::cout<< "Memory cost: " << (double)cpd->get_entry_size()*4/1000000<<" MB"<< std::endl;
}


template<const double* h_period>
void construct_stch_bw_tcpd(const std::string& tpgr_input_file_name,const std::string& btch_input_file_name, const string& output_dir_name){
    string main_name = btch_input_file_name.substr(0,btch_input_file_name.find_last_of("."));
//    //generate dfs ordering using the original graph.
//    std::vector<DEdgeInfo> d_edge_list = katch::tpgr_format::read_edges<DEdgeInfo>(tpgr_input_file_name);
//    D_Graph d_graph(std::move(d_edge_list));
//    std::vector<katch::NodeIterator > dfs_ordering = d_graph.generate_DFS_ordering();
    //generate ch based ordering using the original graph.
//    std::vector<EdgeInfo> edge_list = katch::btch_format::read_edges<EdgeInfo>(btch_input_file_name,period);
//    F_Graph f_graph(std::move(edge_list));
//    auto ordering = katch::btch_format::read_ordering(btch_input_file_name,period);
//    f_graph.test_ordering(ordering);
//    std::vector<katch::NodeIterator > dfs_ordering = f_graph.generate_ch_based_dfs_ordering(ordering);
//
//    string mapper_file = main_name +".bw_htch_mapper";
//    vector<katch::NodeIterator > mapper = invert_permutation(dfs_ordering);
//    save_vector(mapper_file,mapper);
//    KATCH_STATUS("Save mapper to " <<mapper_file<<"\n");
    for(int i = 0; i < 24; i ++){
        // use one mapper only should be fine.
        construct_min_stch_bw_tcpd<h_period>(main_name+"_"+to_string(i)+".btch", output_dir_name,i);
    }
}






template<const double* h_period>
void construct_min_mtch_bw_tcpd(const std::string& btch_input_file_name, const string& output_dir_name, int hour_duration, int shifting_hour,
                                     int hour) {
    using H_Graph = katch::ForwardSearchGraph<h_period>;
    using HEdgeInfo = katch::EdgeInfo<typename H_Graph::TTF>;
    string main_name = btch_input_file_name.substr(0,btch_input_file_name.find_last_of("."));
    std::vector<HEdgeInfo> edge_list = katch::btch_format::read_edges<HEdgeInfo>(btch_input_file_name,*h_period);
    H_Graph h_graph(std::move(edge_list));
    KATCH_STATUS("Graph has " << h_graph.get_n_nodes() << " nodes, " << h_graph.get_n_edges() << " edges.\n");
    S_Graph s_graph(h_graph, H_Graph::min_graph);


    auto ordering = katch::btch_format::read_ordering(btch_input_file_name,*h_period);
    h_graph.test_ordering(ordering);



    my_timer timer1 = my_timer();
    timer1.start();
    // resort the graph based on dfs ordering.
    std::vector<katch::NodeIterator > dfs_ordering = h_graph.generate_ch_based_dfs_ordering(ordering);
    s_graph.resort_graph(dfs_ordering);

    string mapper_file = main_name +"_min.bw_ts_tch_mapper";
    vector<katch::NodeIterator > mapper = invert_permutation(dfs_ordering);
    save_vector(mapper_file,mapper);

    KATCH_STATUS("Building Min CPD ..." <<"\n");
    std::vector<katch::TchDownardDijkstraQuery*> dijkstra;
    for (int i = 0; i < omp_get_max_threads(); i ++){
        dijkstra.push_back(new katch::TchDownardDijkstraQuery(&s_graph));
    }
    unsigned number_of_nodes = s_graph.get_n_nodes();
    auto* cpd = new katch::BW_CPD();
    {
        {
            katch::TchDownardDijkstraQuery dij(&s_graph);
            my_timer t;
            t.start();
            dij.get_first_move(0);
            t.stop();
            double tots = t.elapsed_time_micro()*number_of_nodes / 1000000;
            printf("Estimated sequential running time : %fmin\n", tots / 60.0);

        }

        printf("Using %d threads\n", omp_get_max_threads());
        std::vector<katch::BW_CPD>thread_cpd(omp_get_max_threads());

        int progress = 0;

#pragma omp parallel
        {
            const int thread_count = omp_get_num_threads();
            const int thread_id = omp_get_thread_num();
//            const int begin_int  = 0;
            const int node_count = number_of_nodes;

            int node_begin = (node_count*thread_id) / thread_count ;
            int node_end = (node_count*(thread_id+1)) / thread_count ;
            katch::TchDownardDijkstraQuery& thread_dij = *dijkstra[thread_id];

            for(int source_node=node_begin; source_node < node_end; ++source_node){
                const std::vector<bool>& result = thread_dij.get_reachability( mapper[source_node]);
                thread_cpd[thread_id].append_row( mapper[source_node],result);
#pragma omp critical
                {
                    ++progress;
                    if(progress % 100 == 0) {
                        double ratio = (double)progress / number_of_nodes * 100.0;
                        std::cout << "Progress: [" << progress << "/" << number_of_nodes << "] "
                                  << std::setprecision(3) << ratio << "% \r";
                        std::cout.flush();
                    }
                }
            }
        }

        for(auto&x:thread_cpd)
            cpd->append_rows(x);
    }
    timer1.stop();
    for (auto i: dijkstra) delete i;

    printf("Saving data to cpd.txt \n");
    printf("begin size: %d, entry size: %d\n", cpd->entry_count(), cpd->get_entry_size());


    std::string fname = main_name+"_min.bw_ts_tch_cpd";
    FILE*f = fopen(fname.c_str(), "wb");
    cpd->save(f);
    fclose(f);
    cout << "done" << endl;
    KATCH_STATUS("Save CPD to " <<fname<<"\n");

    std::ofstream myFile(output_dir_name+"/bw_ts_tch_cpd_preprocessing_min_"+ std::to_string(hour_duration) + "_" + std::to_string(shifting_hour) + "_" + std::to_string(hour)+".csv");
    myFile<<"build_time,memory_cost\n";
    myFile<<std::fixed<<setprecision(8)<<timer1.elapsed_time_mins()<<","<<(double)cpd->get_entry_size()*4/1000000
          <<"\n";
    myFile.close();

    std::cout<< "Finished in: " << timer1.elapsed_time_mins()<<" Mins"<< std::endl;
    std::cout<< "Memory cost: " << (double)cpd->get_entry_size()*4/1000000<<" MB"<< std::endl;
}
template<const double* h_period>
void construct_mtch_bw_tcpd(int hour_duration, int shifting_hour, const std::string& tpgr_input_file_name, const std::string& btch_input_file_name, const string& output_dir_name){
    string main_name = btch_input_file_name.substr(0,btch_input_file_name.find_last_of("."));
    for(unsigned hour = 0; hour < 24; hour += hour_duration) {
        string btch_input_file = main_name+ "_" + std::to_string(hour_duration) + "_" + std::to_string(shifting_hour) + "_" + std::to_string(hour)+".btch";
        construct_min_mtch_bw_tcpd<h_period>(btch_input_file, output_dir_name,hour_duration, shifting_hour, hour);
    }
}



int main(int argc, char** argv) {

    const char *binary_name = argv[0];

    if (argc != 4) {
        std::cerr
                << std::endl
                << "USAGE: " << binary_name
                << "<.tpgr input file> <.btch input file> <output file>" << std::endl
                << std::endl;

        return EXIT_FAILURE;
    }
    const std::string tpgr_input_file_name(argv[1]);
    const std::string btch_input_file_name(argv[2]);
    const std::string output_dir_name(argv[3]);
    // build BW CPD for TCH
    construct_min_bw_tcpd_with_ch_ordering(tpgr_input_file_name,btch_input_file_name, output_dir_name);

    // build BW CPD for STCH
    construct_stch_bw_tcpd<&STCH_period>(tpgr_input_file_name,btch_input_file_name, output_dir_name);

    // build BW CPD for MTCH
    construct_mtch_bw_tcpd<&one_MTCH_period>(1,1,tpgr_input_file_name,btch_input_file_name, output_dir_name);
    construct_mtch_bw_tcpd<&two_MTCH_period>(1,2,tpgr_input_file_name,btch_input_file_name, output_dir_name);
    construct_mtch_bw_tcpd<&four_MTCH_period>(1,4,tpgr_input_file_name,btch_input_file_name, output_dir_name);


}