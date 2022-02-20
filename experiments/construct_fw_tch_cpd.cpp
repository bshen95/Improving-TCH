//
// Created by Bojie Shen on 23/5/21.
//

#include <iostream>
#include <stdlib.h>
#include <string>
#include <util/my_timer.h>
#include <unordered_set>
#include <iomanip>

#include "util/misc.h"
#include "io/btch_format.h"
#include "io/tpgr_format.h"
#include "io/edge_info.h"
#include "io/vec_io.h"
#include "datastr/graph/forward_search_graph.h"
#include "datastr/graph/dynamic_search_graph.h"
#include "cpd_prepro/tch_dijkstra_query.h"
#include "omp.h"
#include "datastr/cpd/cpd.h"
#include <map>
#include "cpd_prepro/tch_updown_query.h"
static constexpr double period = 864000.0;
static constexpr double STCH_period = 185000.0;
static constexpr double offset = 5000;
static constexpr double hourly_TCH_period = 36000.0 * 5 + offset;
static constexpr double one_MTCH_period = 36000.0 *2 + offset;
static constexpr double two_MTCH_period = 36000.0 *3 + offset;
static constexpr double four_MTCH_period = 36000.0 *5 + offset;

using S_Graph = katch::StaticForwardSearchGraph;
using F_Graph = katch::ForwardSearchGraph<&period>;
using EdgeInfo = katch::EdgeInfo<katch::ForwardSearchGraph<&period>::TTF>;


using D_Graph = katch::DynamicSearchGraph;
using DEdgeInfo = katch::EdgeInfo<katch::DynamicSearchGraph::TTF>;

katch::CPD* _global_cpd;
vector<double> construction_time;
vector<double> memory_cost;

// using double queue, usually slower than int queue, this is confirmed.

void construct_min_cpd(const std::string& tpgr_input_file_name, const std::string& btch_input_file_name, const string& output_dir_name) {
    string main_name = btch_input_file_name.substr(0,btch_input_file_name.find_last_of("."));
    std::vector<EdgeInfo> edge_list = katch::btch_format::read_edges<EdgeInfo>(btch_input_file_name,period);
    F_Graph f_graph(std::move(edge_list));
    KATCH_STATUS("Graph has " << f_graph.get_n_nodes() << " nodes, " << f_graph.get_n_edges() << " edges.\n");
    S_Graph s_graph(f_graph,F_Graph::min_graph);

    // output the free_flow cost
    std::vector<double> free_flow = s_graph.get_free_flow_vector();
    save_vector(main_name+"_min.ori_fw_tch_free_flow",free_flow);
    KATCH_STATUS("Save free flow to " << main_name << "_min.ori_fw_tch_free_flow"<<"\n");


    //take original ordering.
    std::vector<DEdgeInfo> d_edge_list = katch::tpgr_format::read_edges<DEdgeInfo>(tpgr_input_file_name);
    D_Graph d_graph(std::move(d_edge_list));


    my_timer timer1 = my_timer();
    timer1.start();
    std::vector<katch::NodeIterator > dfs_ordering = d_graph.generate_DFS_ordering();
    // resort the graph based on dfs ordering.;
    s_graph.resort_graph(dfs_ordering);

    string mapper_file = main_name +".ori_fw_tch_mapper";
    vector<katch::NodeIterator > mapper = invert_permutation(dfs_ordering);
    save_vector(mapper_file,mapper);
    KATCH_STATUS("Save mapper to " <<mapper_file<<"\n");



    KATCH_STATUS("Building Min CPD ..." <<"\n");
    std::vector<katch::TchDijkstraQuery*> dijkstra;
    for (int i = 0; i < omp_get_max_threads(); i ++){
        dijkstra.push_back(new katch::TchDijkstraQuery(&s_graph));
    }

    unsigned number_of_nodes = s_graph.get_n_nodes();
    auto* cpd = new katch::CPD();
    {
        {
            katch::TchDijkstraQuery dij(&s_graph);
            my_timer t;
            t.start();
            dij.get_first_move(0);
            t.stop();
            double tots = t.elapsed_time_micro()*number_of_nodes / 1000000/ omp_get_max_threads();
            printf("Estimated sequential running time : %fmin\n", tots / 60.0);

        }

        printf("Using %d threads\n", omp_get_max_threads());
        std::vector<katch::CPD>thread_cpd(omp_get_max_threads());

        int progress = 0;

#pragma omp parallel
        {
            const int thread_count = omp_get_num_threads();
            const int thread_id = omp_get_thread_num();
//            const int begin_int  = 0;
            const int node_count = number_of_nodes;

            int node_begin = (node_count*thread_id) / thread_count ;
            int node_end = (node_count*(thread_id+1)) / thread_count ;
            katch::TchDijkstraQuery& thread_dij = *dijkstra[thread_id];

            for(int source_node=node_begin; source_node < node_end; ++source_node){
                const std::vector<unsigned short>& result = thread_dij.get_first_move(mapper[source_node]);
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
            cpd->append_rows(x);
    }
    timer1.stop();
    for (auto i: dijkstra) delete i;

    printf("Saving data to cpd.txt \n");
    printf("begin size: %d, entry size: %d\n", cpd->entry_count(), cpd->get_entry_size());


    std::string fname = main_name+"_min.ori_fw_tch_cpd";
    FILE*f = fopen(fname.c_str(), "wb");
    cpd->save(f);
    fclose(f);
    cout << "done" << endl;
    KATCH_STATUS("Save CPD to " <<fname<<"\n");



    std::ofstream myFile(output_dir_name+"/ori_fw_tch_cpd_preprocessing_min.csv");
    myFile<<"build_time,memory_cost\n";
    myFile<<std::fixed<<setprecision(8)<<timer1.elapsed_time_mins()<<","<<(double)cpd->get_entry_size()*4/1000000
          <<"\n";
    myFile.close();


    std::cout<< "Finished in: " << timer1.elapsed_time_mins()<<" Mins"<< std::endl;
    std::cout<< "Memory cost: " << (double)cpd->get_entry_size()*4/1000000<<" MB"<< std::endl;
}


void construct_min_tcpd(const std::string& tpgr_input_file_name, const std::string& btch_input_file_name, const string& output_dir_name) {


    string main_name = btch_input_file_name.substr(0,btch_input_file_name.find_last_of("."));
    std::vector<EdgeInfo> edge_list = katch::btch_format::read_edges<EdgeInfo>(btch_input_file_name,period);
    F_Graph f_graph(std::move(edge_list));
    KATCH_STATUS("Graph has " << f_graph.get_n_nodes() << " nodes, " << f_graph.get_n_edges() << " edges.\n");
    S_Graph s_graph(f_graph,F_Graph::min_graph);

    // output the free_flow cost
    std::vector<double> free_flow = s_graph.get_free_flow_vector();
    save_vector(main_name+"_min.fw_tch_free_flow",free_flow);
    KATCH_STATUS("Save free flow to " << main_name << "_min.fw_tch_free_flow"<<"\n");


    //take original ordering.
    std::vector<DEdgeInfo> d_edge_list = katch::tpgr_format::read_edges<DEdgeInfo>(tpgr_input_file_name);
    D_Graph d_graph(std::move(d_edge_list));


    my_timer timer1 = my_timer();
    timer1.start();
    std::vector<katch::NodeIterator > dfs_ordering = d_graph.generate_DFS_ordering();
    // resort the graph based on dfs ordering.;
    s_graph.resort_graph(dfs_ordering);

    string mapper_file = main_name +".fw_tch_mapper";
    vector<katch::NodeIterator > mapper = invert_permutation(dfs_ordering);
    save_vector(mapper_file,mapper);
    KATCH_STATUS("Save mapper to " <<mapper_file<<"\n");



    KATCH_STATUS("Building Min CPD ..." <<"\n");
    std::vector<katch::TchUpDownQuery*> dijkstra;
    for (int i = 0; i < omp_get_max_threads(); i ++){
        dijkstra.push_back(new katch::TchUpDownQuery(&s_graph));
    }

    unsigned number_of_nodes = s_graph.get_n_nodes();
    auto* cpd = new katch::CPD();
    {
        {
            katch::TchUpDownQuery dij(&s_graph);
            my_timer t;
            t.start();
            dij.get_first_move(0);
            t.stop();
            double tots = t.elapsed_time_micro()*number_of_nodes / 1000000/ omp_get_max_threads();
            printf("Estimated sequential running time : %fmin\n", tots / 60.0);

        }

        printf("Using %d threads\n", omp_get_max_threads());
        std::vector<katch::CPD>thread_cpd(omp_get_max_threads());

        int progress = 0;

#pragma omp parallel
        {
            const int thread_count = omp_get_num_threads();
            const int thread_id = omp_get_thread_num();
//            const int begin_int  = 0;
            const int node_count = number_of_nodes;

            int node_begin = (node_count*thread_id) / thread_count ;
            int node_end = (node_count*(thread_id+1)) / thread_count ;
            katch::TchUpDownQuery& thread_dij = *dijkstra[thread_id];

            for(int source_node=node_begin; source_node < node_end; ++source_node){
                const std::vector<unsigned short>& result = thread_dij.get_first_move(mapper[source_node]);
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
            cpd->append_rows(x);
    }
    timer1.stop();
    for (auto i: dijkstra) delete i;

    printf("Saving data to cpd.txt \n");
    printf("begin size: %d, entry size: %d\n", cpd->entry_count(), cpd->get_entry_size());


    std::string fname = main_name+"_min.fw_tch_cpd";
    FILE*f = fopen(fname.c_str(), "wb");
    cpd->save(f);
    fclose(f);
    cout << "done" << endl;
    KATCH_STATUS("Save CPD to " <<fname<<"\n");


    std::ofstream disFile(output_dir_name+"/fw_tch_cpd_distribution_min"+".csv");
    disFile<<"row_id,number_of_entries\n";
    for(int i = 0; i < number_of_nodes; i ++){
        disFile<<i <<","<<cpd->get_row_size(i)<<"\n";
    }
    disFile.close();

    std::ofstream myFile(output_dir_name+"/fw_tch_cpd_preprocessing_min.csv");
    myFile<<"build_time,memory_cost\n";
    myFile<<std::fixed<<setprecision(8)<<timer1.elapsed_time_mins()<<","<<(double)cpd->get_entry_size()*4/1000000
          <<"\n";
    myFile.close();


    std::cout<< "Finished in: " << timer1.elapsed_time_mins()<<" Mins"<< std::endl;
    std::cout<< "Memory cost: " << (double)cpd->get_entry_size()*4/1000000<<" MB"<< std::endl;
}




void construct_min_tch_cpd_down_only(const std::string& tpgr_input_file_name, const std::string& btch_input_file_name, const string& output_dir_name) {


    string main_name = btch_input_file_name.substr(0,btch_input_file_name.find_last_of("."));
    std::vector<EdgeInfo> edge_list = katch::btch_format::read_edges<EdgeInfo>(btch_input_file_name,period);
    F_Graph f_graph(std::move(edge_list));
    KATCH_STATUS("Graph has " << f_graph.get_n_nodes() << " nodes, " << f_graph.get_n_edges() << " edges.\n");
    S_Graph s_graph(f_graph,F_Graph::min_graph);
    auto ordering = katch::btch_format::read_ordering(btch_input_file_name,period);
    f_graph.test_ordering(ordering);

    // output the free_flow cost
    std::vector<double> free_flow = s_graph.get_free_flow_vector();
    save_vector(main_name+"_min.fw_tch_free_flow_down_only",free_flow);
    KATCH_STATUS("Save free flow to " << main_name << "_min.fw_tch_free_flow_down_only"<<"\n");


    //take original ordering.
    std::vector<DEdgeInfo> d_edge_list = katch::tpgr_format::read_edges<DEdgeInfo>(tpgr_input_file_name);
    D_Graph d_graph(std::move(d_edge_list));


    my_timer timer1 = my_timer();
    timer1.start();
    std::vector<katch::NodeIterator > dfs_ordering = d_graph.generate_DFS_ordering();
    // resort the graph based on dfs ordering.;
    s_graph.resort_graph(dfs_ordering);

    string mapper_file = main_name +".fw_tch_mapper_down_only";
    vector<katch::NodeIterator > mapper = invert_permutation(dfs_ordering);
    save_vector(mapper_file,mapper);
    KATCH_STATUS("Save mapper to " <<mapper_file<<"\n");


    vector<unsigned> cpd_2_CH_ranking_mapper = vector<unsigned>(ordering.size());
    for(unsigned i = 0; i < cpd_2_CH_ranking_mapper.size(); i ++){
        cpd_2_CH_ranking_mapper[i] = ordering[dfs_ordering[i]];
    }


    KATCH_STATUS("Building Min CPD Down Only..." <<"\n");
    std::vector<katch::TchDijkstraQuery*> dijkstra;
    for (int i = 0; i < omp_get_max_threads(); i ++){
        dijkstra.push_back(new katch::TchDijkstraQuery(&s_graph, &cpd_2_CH_ranking_mapper));
    }

    unsigned number_of_nodes = s_graph.get_n_nodes();
    auto* cpd = new katch::CPD();
    {
        {
            katch::TchDijkstraQuery dij(&s_graph, &cpd_2_CH_ranking_mapper);
            my_timer t;
            t.start();
            dij.get_first_move_down_only(0);
            t.stop();
            double tots = t.elapsed_time_micro()*number_of_nodes / 1000000/ omp_get_max_threads();
            printf("Estimated sequential running time : %fmin\n", tots / 60.0);

        }

        printf("Using %d threads\n", omp_get_max_threads());
        std::vector<katch::CPD>thread_cpd(omp_get_max_threads());

        int progress = 0;

#pragma omp parallel
        {
            const int thread_count = omp_get_num_threads();
            const int thread_id = omp_get_thread_num();
//            const int begin_int  = 0;
            const int node_count = number_of_nodes;

            int node_begin = (node_count*thread_id) / thread_count ;
            int node_end = (node_count*(thread_id+1)) / thread_count ;
            katch::TchDijkstraQuery& thread_dij = *dijkstra[thread_id];

            for(int source_node=node_begin; source_node < node_end; ++source_node){
                const std::vector<unsigned short>& result = thread_dij.get_first_move_down_only(mapper[source_node]);
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
            cpd->append_rows(x);
    }
    timer1.stop();
    for (auto i: dijkstra) delete i;

    printf("Saving data to cpd.txt \n");
    printf("begin size: %d, entry size: %d\n", cpd->entry_count(), cpd->get_entry_size());


    std::string fname = main_name+"_min.fw_tch_cpd_down_only";
    FILE*f = fopen(fname.c_str(), "wb");
    cpd->save(f);
    fclose(f);
    cout << "done" << endl;
    KATCH_STATUS("Save CPD to " <<fname<<"\n");


    std::ofstream disFile(output_dir_name+"/fw_tch_cpd_distribution_min_down_only.csv");
    disFile<<"row_id,number_of_entries\n";
    for(int i = 0; i < number_of_nodes; i ++){
        disFile<<i <<","<<cpd->get_row_size(i)<<"\n";
    }
    disFile.close();

    std::ofstream myFile(output_dir_name+"/fw_tch_cpd_preprocessing_min_down_only.csv");
    myFile<<"build_time,memory_cost\n";
    myFile<<std::fixed<<setprecision(8)<<timer1.elapsed_time_mins()<<","<<(double)cpd->get_entry_size()*4/1000000
          <<"\n";
    myFile.close();


    std::cout<< "Finished in: " << timer1.elapsed_time_mins()<<" Mins"<< std::endl;
    std::cout<< "Memory cost: " << (double)cpd->get_entry_size()*4/1000000<<" MB"<< std::endl;
}




void construct_min_tch_cpd_with_ch_ordering(const std::string& tpgr_input_file_name, const std::string& btch_input_file_name, const string& output_dir_name) {


    string main_name = btch_input_file_name.substr(0,btch_input_file_name.find_last_of("."));
    std::vector<EdgeInfo> edge_list = katch::btch_format::read_edges<EdgeInfo>(btch_input_file_name,period);
    F_Graph f_graph(std::move(edge_list));
    auto ordering = katch::btch_format::read_ordering(btch_input_file_name,period);
    f_graph.test_ordering(ordering);


    KATCH_STATUS("Graph has " << f_graph.get_n_nodes() << " nodes, " << f_graph.get_n_edges() << " edges.\n");
    S_Graph s_graph(f_graph,F_Graph::min_graph);

    // output the free_flow cost
    std::vector<double> free_flow = s_graph.get_free_flow_vector();
    save_vector(main_name+"_min.fw_tch_free_flow_ch_order",free_flow);
    KATCH_STATUS("Save free flow to " << main_name << "_min.fw_tch_free_flow_ch_order"<<"\n");



    my_timer timer1 = my_timer();
    timer1.start();
    std::vector<katch::NodeIterator > dfs_ordering = f_graph.generate_ch_based_dfs_ordering(ordering);
    // resort the graph based on dfs ordering.;
    s_graph.resort_graph(dfs_ordering);


    KATCH_STATUS("Building Min CPD ..." <<"\n");
    std::vector<katch::TchDijkstraQuery*> dijkstra;
    for (int i = 0; i < omp_get_max_threads(); i ++){
        dijkstra.push_back(new katch::TchDijkstraQuery(&s_graph));
    }
    unsigned number_of_nodes = s_graph.get_n_nodes();
    auto* cpd = new katch::CPD();
    {
        {
            katch::TchDijkstraQuery dij(&s_graph);
            my_timer t;
            t.start();
            dij.get_first_move(0);
            t.stop();
            double tots = t.elapsed_time_micro()*number_of_nodes / 1000000/ omp_get_max_threads();
            printf("Estimated sequential running time : %fmin\n", tots / 60.0);

        }

        printf("Using %d threads\n", omp_get_max_threads());
        std::vector<katch::CPD>thread_cpd(omp_get_max_threads());

        int progress = 0;

#pragma omp parallel
        {
            const int thread_count = omp_get_num_threads();
            const int thread_id = omp_get_thread_num();
//            const int begin_int  = 0;
            const int node_count = number_of_nodes;

            int node_begin = (node_count*thread_id) / thread_count ;
            int node_end = (node_count*(thread_id+1)) / thread_count ;
            katch::TchDijkstraQuery& thread_dij = *dijkstra[thread_id];

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

    string mapper_file = main_name +".fw_tch_mapper_ch_order";
    vector<katch::NodeIterator > mapper = invert_permutation(dfs_ordering);
    save_vector(mapper_file,mapper);
    KATCH_STATUS("Save mapper to " <<mapper_file<<"\n");

    std::string fname = main_name+"_min.fw_tch_cpd_ch_order";
    FILE*f = fopen(fname.c_str(), "wb");
    cpd->save(f);
    fclose(f);
    cout << "done" << endl;
    KATCH_STATUS("Save CPD to " <<fname<<"\n");


    std::ofstream disFile(output_dir_name+"/fw_tch_cpd_distribution_min_ch_order"+".csv");
    disFile<<"row_id,number_of_entries\n";
    for(int i = 0; i < number_of_nodes; i ++){
        disFile<<i <<","<<cpd->get_row_size(i)<<"\n";
    }
    disFile.close();

    std::ofstream myFile(output_dir_name+"/fw_tch_cpd_preprocessing_min_ch_order.csv");
    myFile<<"build_time,memory_cost\n";
    myFile<<std::fixed<<setprecision(8)<<timer1.elapsed_time_mins()<<","<<(double)cpd->get_entry_size()*4/1000000
          <<"\n";
    myFile.close();


    std::cout<< "Finished in: " << timer1.elapsed_time_mins()<<" Mins"<< std::endl;
    std::cout<< "Memory cost: " << (double)cpd->get_entry_size()*4/1000000<<" MB"<< std::endl;
}


template<const double* h_period>
void construct_min_stch_tcpd(const std::string& btch_input_file_name, const string& output_dir_name,
                                   const std::vector<katch::NodeIterator>& dfs_ordering,  const std::vector<katch::NodeIterator>& mapper,int hour) {

    string main_name = btch_input_file_name.substr(0,btch_input_file_name.find_last_of("."));

    using H_Graph = katch::ForwardSearchGraph<h_period>;
    using HEdgeInfo = katch::EdgeInfo<typename H_Graph::TTF>;
    std::vector<HEdgeInfo> edge_list = katch::btch_format::read_edges<HEdgeInfo>(btch_input_file_name,*h_period);
    H_Graph h_graph(std::move(edge_list));
    KATCH_STATUS("Graph has " << h_graph.get_n_nodes() << " nodes, " << h_graph.get_n_edges() << " edges.\n");
    S_Graph s_graph(h_graph,H_Graph::min_graph);

    my_timer timer1 = my_timer();
    timer1.start();
    // resort the graph based on dfs ordering.
    s_graph.resort_graph(dfs_ordering);

    KATCH_STATUS("Building Min CPD ..." <<"\n");

    std::vector<katch::TchUpDownQuery*> dijkstra;
    for (int i = 0; i < omp_get_max_threads(); i ++){
        dijkstra.push_back(new katch::TchUpDownQuery(&s_graph));
    }
    unsigned number_of_nodes = s_graph.get_n_nodes();
    auto* cpd = new katch::CPD();
    {
        {
            katch::TchUpDownQuery dij(&s_graph);
            my_timer t;
            t.start();
            dij.get_first_move(0);
            t.stop();
            double tots = t.elapsed_time_micro()*number_of_nodes/10 /1000000 / omp_get_max_threads();
            printf("Estimated sequential running time : %fmin\n", tots / 60.0);

        }

        printf("Using %d threads\n", omp_get_max_threads());
        std::vector<katch::CPD>thread_cpd(omp_get_max_threads());

        int progress = 0;

#pragma omp parallel
        {
            const int thread_count = omp_get_num_threads();
            const int thread_id = omp_get_thread_num();
//            const int begin_int  = 0;
            const int node_count = number_of_nodes;

            int node_begin = (node_count*thread_id) / thread_count ;
            int node_end = (node_count*(thread_id+1)) / thread_count ;
            katch::TchUpDownQuery& thread_dij = *dijkstra[thread_id];

            for(int source_node=node_begin; source_node < node_end; ++source_node){
                const std::vector<unsigned short>& result = thread_dij.get_first_move( mapper[source_node]);
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
            cpd->append_rows(x);
    }
    timer1.stop();
    for (auto i: dijkstra) delete i;

    printf("Saving data to cpd.txt \n");
    printf("begin size: %d, entry size: %d\n", cpd->entry_count(), cpd->get_entry_size());


    std::string fname = main_name+"_min.fw_htch_cpd";
    FILE*f = fopen(fname.c_str(), "wb");
    cpd->save(f);
    fclose(f);
    cout << "done" << endl;
    KATCH_STATUS("Save CPD to " <<fname<<"\n");

    std::ofstream myFile(output_dir_name+"/fw_htch_cpd_preprocessing_min_"+to_string(hour)+".csv");
    myFile<<"build_time,memory_cost\n";
    myFile<<std::fixed<<setprecision(8)<<timer1.elapsed_time_mins()<<","<<(double)cpd->get_entry_size()*4/1000000
          <<"\n";
    myFile.close();


    std::cout<< "Finished in: " << timer1.elapsed_time_mins()<<" Mins"<< std::endl;
    std::cout<< "Memory cost: " << (double)cpd->get_entry_size()*4/1000000<<" MB"<< std::endl;
}




template<const double* h_period>
void construct_stch_tcpd(const std::string& tpgr_input_file_name, const std::string& btch_input_file_name, const string& output_dir_name){
    string main_name = btch_input_file_name.substr(0,btch_input_file_name.find_last_of("."));
    std::vector<DEdgeInfo> d_edge_list = katch::tpgr_format::read_edges<DEdgeInfo>(tpgr_input_file_name);
    D_Graph d_graph(std::move(d_edge_list));
    std::vector<katch::NodeIterator > dfs_ordering = d_graph.generate_DFS_ordering();
    string mapper_file = main_name +".fw_htch_mapper";
    vector<katch::NodeIterator > mapper = invert_permutation(dfs_ordering);
    save_vector(mapper_file,mapper);
    KATCH_STATUS("Save mapper to " <<mapper_file<<"\n");
    for(int i = 0; i < 24; i ++){
        // use one mapper only should be fine.
        construct_min_stch_tcpd<h_period>(main_name+"_"+to_string(i)+".btch", output_dir_name,dfs_ordering,mapper,i);
    }
}













template<const double* h_period>
void construct_min_mtch_tcpd(const std::string& btch_input_file_name, const string& output_dir_name,
                                     const std::vector<katch::NodeIterator>& dfs_ordering,  const std::vector<katch::NodeIterator>& mapper,
                                     int hour_duration, int shifting_hour,
                                     int hour) {

    string main_name = btch_input_file_name.substr(0,btch_input_file_name.find_last_of("."));

    using H_Graph = katch::ForwardSearchGraph<h_period>;
    using HEdgeInfo = katch::EdgeInfo<typename H_Graph::TTF>;
    std::vector<HEdgeInfo> edge_list = katch::btch_format::read_edges<HEdgeInfo>(btch_input_file_name,*h_period);
    H_Graph h_graph(std::move(edge_list));
    KATCH_STATUS("Graph has " << h_graph.get_n_nodes() << " nodes, " << h_graph.get_n_edges() << " edges.\n");
    S_Graph s_graph(h_graph,H_Graph::min_graph);

    my_timer timer1 = my_timer();
    timer1.start();
    // resort the graph based on dfs ordering.
    s_graph.resort_graph(dfs_ordering);

    KATCH_STATUS("Building Min CPD ..." <<"\n");

    std::vector<katch::TchUpDownQuery*> dijkstra;
    for (int i = 0; i < omp_get_max_threads(); i ++){
        dijkstra.push_back(new katch::TchUpDownQuery(&s_graph));
    }
    unsigned number_of_nodes = s_graph.get_n_nodes();
    auto* cpd = new katch::CPD();
    {
        {
            katch::TchUpDownQuery dij(&s_graph);
            my_timer t;
            t.start();
            dij.get_first_move(0);
            t.stop();
            double tots = t.elapsed_time_micro()*number_of_nodes/10 /1000000 / omp_get_max_threads();
            printf("Estimated sequential running time : %fmin\n", tots / 60.0);

        }

        printf("Using %d threads\n", omp_get_max_threads());
        std::vector<katch::CPD>thread_cpd(omp_get_max_threads());

        int progress = 0;

#pragma omp parallel
        {
            const int thread_count = omp_get_num_threads();
            const int thread_id = omp_get_thread_num();
//            const int begin_int  = 0;
            const int node_count = number_of_nodes;

            int node_begin = (node_count*thread_id) / thread_count ;
            int node_end = (node_count*(thread_id+1)) / thread_count ;
            katch::TchUpDownQuery& thread_dij = *dijkstra[thread_id];

            for(int source_node=node_begin; source_node < node_end; ++source_node){
                const std::vector<unsigned short>& result = thread_dij.get_first_move( mapper[source_node]);
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
            cpd->append_rows(x);
    }
    timer1.stop();
    for (auto i: dijkstra) delete i;

    printf("Saving data to cpd.txt \n");
    printf("begin size: %d, entry size: %d\n", cpd->entry_count(), cpd->get_entry_size());


    std::string fname = main_name+"_min.fw_ts_tch_cpd";
    FILE*f = fopen(fname.c_str(), "wb");
    cpd->save(f);
    fclose(f);
    cout << "done" << endl;
    KATCH_STATUS("Save CPD to " <<fname<<"\n");

    std::ofstream myFile(output_dir_name+"/fw_ts_tch_cpd_preprocessing_min_"+ std::to_string(hour_duration) + "_" + std::to_string(shifting_hour) + "_" + std::to_string(hour)+".csv");
    myFile<<"build_time,memory_cost\n";
    myFile<<std::fixed<<setprecision(8)<<timer1.elapsed_time_mins()<<","<<(double)cpd->get_entry_size()*4/1000000
          <<"\n";
    myFile.close();


    std::cout<< "Finished in: " << timer1.elapsed_time_mins()<<" Mins"<< std::endl;
    std::cout<< "Memory cost: " << (double)cpd->get_entry_size()*4/1000000<<" MB"<< std::endl;
}



template<const double* h_period>
void construct_mtch_tcpd(int hour_duration, int shifting_hour, const std::string& tpgr_input_file_name, const std::string& btch_input_file_name, const string& output_dir_name){
    string main_name = btch_input_file_name.substr(0,btch_input_file_name.find_last_of("."));
    std::vector<DEdgeInfo> d_edge_list = katch::tpgr_format::read_edges<DEdgeInfo>(tpgr_input_file_name);
    D_Graph d_graph(std::move(d_edge_list));
    std::vector<katch::NodeIterator > dfs_ordering = d_graph.generate_DFS_ordering();
    string mapper_file = main_name + "_" + std::to_string(hour_duration) + "_" + std::to_string(shifting_hour) + ".fw_ts_tch_mapper";
    vector<katch::NodeIterator > mapper = invert_permutation(dfs_ordering);
    save_vector(mapper_file,mapper);
    KATCH_STATUS("Save mapper to " <<mapper_file<<"\n");
    for(unsigned hour = 0; hour < 24; hour += hour_duration) {
        string btch_input_file = main_name+ "_" + std::to_string(hour_duration) + "_" + std::to_string(shifting_hour) + "_" + std::to_string(hour)+".btch";

        // use one mapper only should be fine.
        construct_min_mtch_tcpd<h_period>(btch_input_file, output_dir_name,dfs_ordering,mapper,hour_duration, shifting_hour, hour);
    }
}















void construct_global_gtch_cpd(const std::string& tpgr_input_file_name, const std::string& btch_input_file_name, const string& output_dir_name) {

    string main_name = btch_input_file_name.substr(0,btch_input_file_name.find_last_of("."));
    std::vector<EdgeInfo> edge_list = katch::btch_format::read_edges<EdgeInfo>(btch_input_file_name,period);
    F_Graph f_graph(std::move(edge_list));
    KATCH_STATUS("Graph has " << f_graph.get_n_nodes() << " nodes, " << f_graph.get_n_edges() << " edges.\n");
    S_Graph s_graph(f_graph,F_Graph::min_graph);

    // output the free_flow cost
    std::vector<double> free_flow = s_graph.get_free_flow_vector();
    save_vector(main_name+"_global.fw_gtch_free_flow",free_flow);
    KATCH_STATUS("Save free flow to " << main_name << "_global.fw_gtch_free_flow"<<"\n");


    //take original ordering.
    std::vector<DEdgeInfo> d_edge_list = katch::tpgr_format::read_edges<DEdgeInfo>(tpgr_input_file_name);
    D_Graph d_graph(std::move(d_edge_list));

    my_timer timer1 = my_timer();
    timer1.start();
    std::vector<katch::NodeIterator > dfs_ordering = d_graph.generate_DFS_ordering();
    // resort the graph based on dfs ordering.;
    s_graph.resort_graph(dfs_ordering);


    KATCH_STATUS("Building Min CPD ..." <<"\n");
    std::vector<katch::TchDijkstraQuery*> dijkstra;
    for (int i = 0; i < omp_get_max_threads(); i ++){
        dijkstra.push_back(new katch::TchDijkstraQuery(&s_graph));
    }
    unsigned number_of_nodes = s_graph.get_n_nodes();
    auto* cpd = new katch::CPD();
    {
        {
            katch::TchDijkstraQuery dij(&s_graph);
            my_timer t;
            t.start();
            dij.get_first_move(0);
            t.stop();
            double tots = t.elapsed_time_micro()*number_of_nodes / 1000000/ omp_get_max_threads();
            printf("Estimated sequential running time : %fmin\n", tots / 60.0);

        }

        printf("Using %d threads\n", omp_get_max_threads());
        std::vector<katch::CPD>thread_cpd(omp_get_max_threads());

        int progress = 0;

#pragma omp parallel
        {
            const int thread_count = omp_get_num_threads();
            const int thread_id = omp_get_thread_num();
//            const int begin_int  = 0;
            const int node_count = number_of_nodes;

            int node_begin = (node_count*thread_id) / thread_count ;
            int node_end = (node_count*(thread_id+1)) / thread_count ;
            katch::TchDijkstraQuery& thread_dij = *dijkstra[thread_id];

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

    string mapper_file = main_name +"_global.fw_gtch_mapper";
    vector<katch::NodeIterator > mapper = invert_permutation(dfs_ordering);
    save_vector(mapper_file,mapper);
    KATCH_STATUS("Save mapper to " <<mapper_file<<"\n");

    std::string fname = main_name+"_global.fw_gtch_cpd";
    FILE*f = fopen(fname.c_str(), "wb");
    cpd->save(f);
    fclose(f);
    cout << "done" << endl;
    KATCH_STATUS("Save CPD to " <<fname<<"\n");


    std::ofstream disFile(output_dir_name+"/fw_gtch_cpd_distribution_global"+".csv");
    disFile<<"row_id,number_of_entries\n";
    for(int i = 0; i < number_of_nodes; i ++){
        disFile<<i <<","<<cpd->get_row_size(i)<<"\n";
    }
    disFile.close();

    std::ofstream myFile(output_dir_name+"/fw_gtch_cpd_preprocessing_global.csv");
    myFile<<"build_time,memory_cost\n";
    myFile<<std::fixed<<setprecision(8)<<timer1.elapsed_time_mins()<<","<<(double)cpd->get_entry_size()*4/1000000
          <<"\n";
    myFile.close();


    std::cout<< "Finished in: " << timer1.elapsed_time_mins()<<" Mins"<< std::endl;
    std::cout<< "Memory cost: " << (double)cpd->get_entry_size()*4/1000000<<" MB"<< std::endl;
}


void construct_hourly_gtch_cpd(const std::string& btch_input_file_name, const string& output_dir_name, const std::vector<katch::NodeIterator>& dfs_ordering, int hour) {


    string main_name = btch_input_file_name.substr(0,btch_input_file_name.find_last_of("."));
    std::vector<EdgeInfo> edge_list = katch::btch_format::read_edges<EdgeInfo>(btch_input_file_name,period);
    F_Graph f_graph(std::move(edge_list));
    KATCH_STATUS("Graph has " << f_graph.get_n_nodes() << " nodes, " << f_graph.get_n_edges() << " edges.\n");
    S_Graph s_graph(f_graph,{hour*36000.0,hour*36000.0+3*36000.0});

    // output the free_flow cost
    std::vector<double> free_flow = s_graph.get_free_flow_vector();
    save_vector(main_name+"_" + to_string(hour) + ".fw_gtch_free_flow",free_flow);
    KATCH_STATUS("Save free flow to " << main_name << "_" + to_string(hour) + ".fw_gtch_free_flow"<<"\n");



    my_timer timer1 = my_timer();
    timer1.start();
    // resort the graph based on dfs ordering.
    s_graph.resort_graph(dfs_ordering);


    KATCH_STATUS("Building Min CPD ..." <<"\n");

    std::vector<katch::TchDijkstraQuery*> dijkstra;
    for (int i = 0; i < omp_get_max_threads(); i ++){
        dijkstra.push_back(new katch::TchDijkstraQuery(&s_graph));
    }
    unsigned number_of_nodes = s_graph.get_n_nodes();
    std::vector<unsigned>number_of_same_symbols(number_of_nodes,0);
    auto* cpd = new katch::CPD();
    {
        {
            katch::TchDijkstraQuery dij(&s_graph);
            my_timer t;
            t.start();
            dij.get_first_move(0);
            t.stop();
            double tots = t.elapsed_time_micro()*number_of_nodes/10 /1000000 / omp_get_max_threads();
            printf("Estimated sequential running time : %fmin\n", tots / 60.0);

        }

        printf("Using %d threads\n", omp_get_max_threads());
        std::vector<katch::CPD>thread_cpd(omp_get_max_threads());
        int progress = 0;

#pragma omp parallel
        {
            const int thread_count = omp_get_num_threads();
            const int thread_id = omp_get_thread_num();
//            const int begin_int  = 0;
            const int node_count = number_of_nodes;

            int node_begin = (node_count*thread_id) / thread_count ;
            int node_end = (node_count*(thread_id+1)) / thread_count ;
            katch::TchDijkstraQuery& thread_dij = *dijkstra[thread_id];

            for(int source_node=node_begin; source_node < node_end; ++source_node){
                std::vector<unsigned short> result = thread_dij.get_first_move(source_node);
                const std::vector<unsigned short> decode_results = _global_cpd->decode_row(source_node);
                std::vector<unsigned short> cpd_result (number_of_nodes*2);
                for(int i = 0; i < number_of_nodes;i++){
                    cpd_result[i*2] = result[i];
                    if(result[i] == decode_results[i]){
                        // 0XFD for same
                        cpd_result[i * 2 + 1 ] =0XFD;
                        number_of_same_symbols[source_node]++;
                    }else{
                        // 0XFC for not same
                        cpd_result[i * 2 + 1 ] =0XFC;
                    }
                }
//                cpd_result[source_node] = unordered_set<unsigned short>();
//                thread_cpd[thread_id].append_row(source_node,result);
                thread_cpd[thread_id].append_row_multiple_symbols(source_node,cpd_result);
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

        for(auto&x:thread_cpd){
            cpd->append_rows(x);
        }
    }
    timer1.stop();
    for (auto i: dijkstra) delete i;

    unsigned long total_same = 0;
    for (auto i : number_of_same_symbols ){
        total_same += i;
    }
    std::cout<<"Percentage of same symbol: " << (double)total_same/(dfs_ordering.size()*dfs_ordering.size()) <<std::endl;
    printf("Saving data to cpd.txt \n");
    printf("begin size: %d, entry size: %d\n", cpd->entry_count(), cpd->get_entry_size());


    std::string fname = main_name+"_"+ to_string(hour) +".fw_gtch_cpd";
    FILE*f = fopen(fname.c_str(), "wb");
    cpd->save(f);
    fclose(f);
    cout << "done" << endl;
    KATCH_STATUS("Save CPD to " <<fname<<"\n");

    std::ofstream disFile(output_dir_name+"/fw_gtch_cpd_distribution_"+to_string(hour)+".csv");
    disFile<<"row_id,number_of_entries,number_of_same_symbols\n";
    for(int i = 0; i < number_of_nodes; i ++){
        disFile<<i <<","<<cpd->get_row_size(i)<<","<< number_of_same_symbols[i]<<"\n";
    }
    disFile.close();


    std::ofstream myFile(output_dir_name+"/fw_gtch_cpd_preprocessing_"+to_string(hour)+".csv");
    myFile<<"build_time,memory_cost\n";
    myFile<<std::fixed<<setprecision(8)<<timer1.elapsed_time_mins()<<","<<(double)cpd->get_entry_size()*4/1000000
          <<"\n";
    myFile.close();


    std::cout<< "Finished in: " << timer1.elapsed_time_mins()<<" Mins"<< std::endl;
    std::cout<< "Memory cost: " << (double)cpd->get_entry_size()*4/1000000<<" MB"<< std::endl;
}

void construct_hourly_gtch_cpd(const std::string& tpgr_input_file_name, const std::string& btch_input_file_name, const string& output_dir_name){

//    construct_global_gtch_cpd(tpgr_input_file_name, btch_input_file_name, output_dir_name);

    string main_name = btch_input_file_name.substr(0,btch_input_file_name.find_last_of("."));
    std::vector<katch::NodeIterator > dfs_ordering = load_vector<katch::NodeIterator >(main_name+"_global.fw_gtch_mapper");
    dfs_ordering = invert_permutation(dfs_ordering);
    _global_cpd = new katch::CPD();
    FILE*f = fopen((main_name +"_global.fw_gtch_cpd").c_str(), "r");
    _global_cpd->load(f);

    for(int i = 0; i < 24; i ++){
        // use one mapper only should be fine.
        construct_hourly_gtch_cpd(btch_input_file_name, output_dir_name,dfs_ordering,i);
    }
}








void construct_global_gtch_cpd_average(const std::string& tpgr_input_file_name, const std::string& btch_input_file_name, const string& output_dir_name) {

    string main_name = btch_input_file_name.substr(0,btch_input_file_name.find_last_of("."));
    std::vector<EdgeInfo> edge_list = katch::btch_format::read_edges<EdgeInfo>(btch_input_file_name,period);
    F_Graph f_graph(std::move(edge_list));
    KATCH_STATUS("Graph has " << f_graph.get_n_nodes() << " nodes, " << f_graph.get_n_edges() << " edges.\n");
    S_Graph s_graph(f_graph,F_Graph::average_graph);

    // output the free_flow cost
    std::vector<double> free_flow = s_graph.get_free_flow_vector();
    save_vector(main_name+"_global.fw_gtch_free_flow_avg",free_flow);
    KATCH_STATUS("Save free flow to " << main_name << "_global.fw_gtch_free_flow_avg"<<"\n");


    //take original ordering.
    std::vector<DEdgeInfo> d_edge_list = katch::tpgr_format::read_edges<DEdgeInfo>(tpgr_input_file_name);
    D_Graph d_graph(std::move(d_edge_list));

    my_timer timer1 = my_timer();
    timer1.start();
    std::vector<katch::NodeIterator > dfs_ordering = d_graph.generate_DFS_ordering();
    // resort the graph based on dfs ordering.;
    s_graph.resort_graph(dfs_ordering);


    KATCH_STATUS("Building Min CPD ..." <<"\n");
    std::vector<katch::TchDijkstraQuery*> dijkstra;
    for (int i = 0; i < omp_get_max_threads(); i ++){
        dijkstra.push_back(new katch::TchDijkstraQuery(&s_graph));
    }
    unsigned number_of_nodes = s_graph.get_n_nodes();
    auto* cpd = new katch::CPD();
    {
        {
            katch::TchDijkstraQuery dij(&s_graph);
            my_timer t;
            t.start();
            dij.get_first_move(0);
            t.stop();
            double tots = t.elapsed_time_micro()*number_of_nodes / 1000000/ omp_get_max_threads();
            printf("Estimated sequential running time : %fmin\n", tots / 60.0);

        }

        printf("Using %d threads\n", omp_get_max_threads());
        std::vector<katch::CPD>thread_cpd(omp_get_max_threads());

        int progress = 0;

#pragma omp parallel
        {
            const int thread_count = omp_get_num_threads();
            const int thread_id = omp_get_thread_num();
//            const int begin_int  = 0;
            const int node_count = number_of_nodes;

            int node_begin = (node_count*thread_id) / thread_count ;
            int node_end = (node_count*(thread_id+1)) / thread_count ;
            katch::TchDijkstraQuery& thread_dij = *dijkstra[thread_id];

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

    string mapper_file = main_name +"_global.fw_gtch_mapper_avg";
    vector<katch::NodeIterator > mapper = invert_permutation(dfs_ordering);
    save_vector(mapper_file,mapper);
    KATCH_STATUS("Save mapper to " <<mapper_file<<"\n");

    std::string fname = main_name+"_global.fw_gtch_cpd_avg";
    FILE*f = fopen(fname.c_str(), "wb");
    cpd->save(f);
    fclose(f);
    cout << "done" << endl;
    KATCH_STATUS("Save CPD to " <<fname<<"\n");


    std::ofstream disFile(output_dir_name+"/fw_gtch_cpd_distribution_global_avg.csv");
    disFile<<"row_id,number_of_entries\n";
    for(int i = 0; i < number_of_nodes; i ++){
        disFile<<i <<","<<cpd->get_row_size(i)<<"\n";
    }
    disFile.close();

    std::ofstream myFile(output_dir_name+"/fw_gtch_cpd_preprocessing_global_avg.csv");
    myFile<<"build_time,memory_cost\n";
    myFile<<std::fixed<<setprecision(8)<<timer1.elapsed_time_mins()<<","<<(double)cpd->get_entry_size()*4/1000000
          <<"\n";
    myFile.close();


    std::cout<< "Finished in: " << timer1.elapsed_time_mins()<<" Mins"<< std::endl;
    std::cout<< "Memory cost: " << (double)cpd->get_entry_size()*4/1000000<<" MB"<< std::endl;
}


void construct_hourly_gtch_cpd_average(const std::string& btch_input_file_name, const string& output_dir_name, const std::vector<katch::NodeIterator>& dfs_ordering, int hour) {


    string main_name = btch_input_file_name.substr(0,btch_input_file_name.find_last_of("."));
    std::vector<EdgeInfo> edge_list = katch::btch_format::read_edges<EdgeInfo>(btch_input_file_name,period);
    F_Graph f_graph(std::move(edge_list));
    KATCH_STATUS("Graph has " << f_graph.get_n_nodes() << " nodes, " << f_graph.get_n_edges() << " edges.\n");
    S_Graph s_graph(f_graph,{hour*36000.0,hour*36000.0+3*36000.0});

    // output the free_flow cost
    std::vector<double> free_flow = s_graph.get_free_flow_vector();
    save_vector(main_name+"_" + to_string(hour) + ".fw_gtch_free_flow_avg",free_flow);
    KATCH_STATUS("Save free flow to " << main_name << "_" + to_string(hour) + ".fw_gtch_free_flow_avg"<<"\n");



    my_timer timer1 = my_timer();
    timer1.start();
    // resort the graph based on dfs ordering.
    s_graph.resort_graph(dfs_ordering);


    KATCH_STATUS("Building Min CPD ..." <<"\n");

    std::vector<katch::TchDijkstraQuery*> dijkstra;
    for (int i = 0; i < omp_get_max_threads(); i ++){
        dijkstra.push_back(new katch::TchDijkstraQuery(&s_graph));
    }
    unsigned number_of_nodes = s_graph.get_n_nodes();
    std::vector<unsigned>number_of_same_symbols(number_of_nodes,0);
    auto* cpd = new katch::CPD();
    {
        {
            katch::TchDijkstraQuery dij(&s_graph);
            my_timer t;
            t.start();
            dij.get_first_move(0);
            t.stop();
            double tots = t.elapsed_time_micro()*number_of_nodes/10 /1000000 / omp_get_max_threads();
            printf("Estimated sequential running time : %fmin\n", tots / 60.0);

        }

        printf("Using %d threads\n", omp_get_max_threads());
        std::vector<katch::CPD>thread_cpd(omp_get_max_threads());
        int progress = 0;

#pragma omp parallel
        {
            const int thread_count = omp_get_num_threads();
            const int thread_id = omp_get_thread_num();
//            const int begin_int  = 0;
            const int node_count = number_of_nodes;

            int node_begin = (node_count*thread_id) / thread_count ;
            int node_end = (node_count*(thread_id+1)) / thread_count ;
            katch::TchDijkstraQuery& thread_dij = *dijkstra[thread_id];

            for(int source_node=node_begin; source_node < node_end; ++source_node){
                std::vector<unsigned short> result = thread_dij.get_first_move(source_node);
                const std::vector<unsigned short> decode_results = _global_cpd->decode_row(source_node);
                std::vector<unsigned short> cpd_result (number_of_nodes*2);
                for(int i = 0; i < number_of_nodes;i++){
                    cpd_result[i*2] = result[i];
                    if(result[i] == decode_results[i]){
                        // 0XFD for same
                        cpd_result[i * 2 + 1 ] =0XFD;
                        number_of_same_symbols[source_node]++;
                    }else{
                        // 0XFC for not same
                        cpd_result[i * 2 + 1 ] =0XFC;
                    }
                }
//                cpd_result[source_node] = unordered_set<unsigned short>();
//                thread_cpd[thread_id].append_row(source_node,result);
                thread_cpd[thread_id].append_row_multiple_symbols(source_node,cpd_result);
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

        for(auto&x:thread_cpd){
            cpd->append_rows(x);
        }
    }
    timer1.stop();
    for (auto i: dijkstra) delete i;

    unsigned long total_same = 0;
    for (auto i : number_of_same_symbols ){
        total_same += i;
    }
    std::cout<<"Percentage of same symbol: " << (double)total_same/(dfs_ordering.size()*dfs_ordering.size()) <<std::endl;
    printf("Saving data to cpd.txt \n");
    printf("begin size: %d, entry size: %d\n", cpd->entry_count(), cpd->get_entry_size());


    std::string fname = main_name+"_"+ to_string(hour) +".fw_gtch_cpd_avg";
    FILE*f = fopen(fname.c_str(), "wb");
    cpd->save(f);
    fclose(f);
    cout << "done" << endl;
    KATCH_STATUS("Save CPD to " <<fname<<"\n");

    std::ofstream disFile(output_dir_name+"/fw_gtch_cpd_distribution_"+to_string(hour)+"_avg.csv");
    disFile<<"row_id,number_of_entries,number_of_same_symbols\n";
    for(int i = 0; i < number_of_nodes; i ++){
        disFile<<i <<","<<cpd->get_row_size(i)<<","<< number_of_same_symbols[i]<<"\n";
    }
    disFile.close();


    std::ofstream myFile(output_dir_name+"/fw_gtch_cpd_preprocessing_"+to_string(hour)+"_avg.csv");
    myFile<<"build_time,memory_cost\n";
    myFile<<std::fixed<<setprecision(8)<<timer1.elapsed_time_mins()<<","<<(double)cpd->get_entry_size()*4/1000000
          <<"\n";
    myFile.close();


    std::cout<< "Finished in: " << timer1.elapsed_time_mins()<<" Mins"<< std::endl;
    std::cout<< "Memory cost: " << (double)cpd->get_entry_size()*4/1000000<<" MB"<< std::endl;
}


void construct_hourly_gtch_cpd_average(const std::string& tpgr_input_file_name, const std::string& btch_input_file_name, const string& output_dir_name){

    construct_global_gtch_cpd_average(tpgr_input_file_name, btch_input_file_name, output_dir_name);

    string main_name = btch_input_file_name.substr(0,btch_input_file_name.find_last_of("."));
    std::vector<katch::NodeIterator > dfs_ordering = load_vector<katch::NodeIterator >(main_name+"_global.fw_gtch_mapper_avg");
    dfs_ordering = invert_permutation(dfs_ordering);
    _global_cpd = new katch::CPD();
    FILE*f = fopen((main_name +"_global.fw_gtch_cpd_avg").c_str(), "r");
    _global_cpd->load(f);

    for(int i = 0; i < 24; i ++){
        // use one mapper only should be fine.
        construct_hourly_gtch_cpd_average(btch_input_file_name, output_dir_name,dfs_ordering,i);
    }
}




void construct_global_gtch_cpd_4h(const std::string& btch_input_file_name, const string& output_dir_name, const std::vector<katch::NodeIterator>& dfs_ordering, int hour) {

    string main_name = btch_input_file_name.substr(0,btch_input_file_name.find_last_of("."));
    std::vector<EdgeInfo> edge_list = katch::btch_format::read_edges<EdgeInfo>(btch_input_file_name,period);
    F_Graph f_graph(std::move(edge_list));
    KATCH_STATUS("Graph has " << f_graph.get_n_nodes() << " nodes, " << f_graph.get_n_edges() << " edges.\n");
    S_Graph s_graph(f_graph,{hour*36000.0,hour*36000.0+3*36000.0});

    // output the free_flow cost
    std::vector<double> free_flow = s_graph.get_free_flow_vector();
    save_vector(main_name+"_" + to_string(hour) + ".fw_gtch_free_flow_4h",free_flow);
    KATCH_STATUS("Save free flow to " << main_name << "_" + to_string(hour) + ".fw_gtch_free_flow_4h"<<"\n");



    my_timer timer1 = my_timer();
    timer1.start();
    // resort the graph based on dfs ordering.
    s_graph.resort_graph(dfs_ordering);


    KATCH_STATUS("Building Min CPD ..." <<"\n");

    std::vector<katch::TchDijkstraQuery*> dijkstra;
    for (int i = 0; i < omp_get_max_threads(); i ++){
        dijkstra.push_back(new katch::TchDijkstraQuery(&s_graph));
    }
    unsigned number_of_nodes = s_graph.get_n_nodes();
    std::vector<unsigned>number_of_same_symbols(number_of_nodes,0);
    auto* cpd = new katch::CPD();
    {
        {
            katch::TchDijkstraQuery dij(&s_graph);
            my_timer t;
            t.start();
            dij.get_first_move(0);
            t.stop();
            double tots = t.elapsed_time_micro()*number_of_nodes/10 /1000000 / omp_get_max_threads();
            printf("Estimated sequential running time : %fmin\n", tots / 60.0);

        }

        printf("Using %d threads\n", omp_get_max_threads());
        std::vector<katch::CPD>thread_cpd(omp_get_max_threads());
        int progress = 0;

#pragma omp parallel
        {
            const int thread_count = omp_get_num_threads();
            const int thread_id = omp_get_thread_num();
//            const int begin_int  = 0;
            const int node_count = number_of_nodes;

            int node_begin = (node_count*thread_id) / thread_count ;
            int node_end = (node_count*(thread_id+1)) / thread_count ;
            katch::TchDijkstraQuery& thread_dij = *dijkstra[thread_id];

            for(int source_node=node_begin; source_node < node_end; ++source_node){
                std::vector<unsigned short> result = thread_dij.get_first_move(source_node);
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

        for(auto&x:thread_cpd){
            cpd->append_rows(x);
        }
    }
    timer1.stop();
    for (auto i: dijkstra) delete i;

    unsigned long total_same = 0;
    for (auto i : number_of_same_symbols ){
        total_same += i;
    }
    std::cout<<"Percentage of same symbol: " << (double)total_same/(dfs_ordering.size()*dfs_ordering.size()) <<std::endl;
    printf("Saving data to cpd.txt \n");
    printf("begin size: %d, entry size: %d\n", cpd->entry_count(), cpd->get_entry_size());


    std::string fname = main_name+"_"+ to_string(hour) +".fw_gtch_cpd_4h";
    FILE*f = fopen(fname.c_str(), "wb");
    cpd->save(f);
    fclose(f);
    cout << "done" << endl;
    KATCH_STATUS("Save CPD to " <<fname<<"\n");

    std::ofstream disFile(output_dir_name+"/fw_gtch_cpd_distribution_"+to_string(hour)+"_4h.csv");
    disFile<<"row_id,number_of_entries,number_of_same_symbols\n";
    for(int i = 0; i < number_of_nodes; i ++){
        disFile<<i <<","<<cpd->get_row_size(i)<<","<< number_of_same_symbols[i]<<"\n";
    }
    disFile.close();


    std::ofstream myFile(output_dir_name+"/fw_gtch_cpd_preprocessing_"+to_string(hour)+"_4h.csv");
    myFile<<"build_time,memory_cost\n";
    myFile<<std::fixed<<setprecision(8)<<timer1.elapsed_time_mins()<<","<<(double)cpd->get_entry_size()*4/1000000
          <<"\n";
    myFile.close();


    std::cout<< "Finished in: " << timer1.elapsed_time_mins()<<" Mins"<< std::endl;
    std::cout<< "Memory cost: " << (double)cpd->get_entry_size()*4/1000000<<" MB"<< std::endl;
}



void construct_hourly_gtch_cpd_4h(const std::string& btch_input_file_name, const string& output_dir_name, const std::vector<katch::NodeIterator>& dfs_ordering, int hour) {


    string main_name = btch_input_file_name.substr(0,btch_input_file_name.find_last_of("."));
    std::vector<EdgeInfo> edge_list = katch::btch_format::read_edges<EdgeInfo>(btch_input_file_name,period);
    F_Graph f_graph(std::move(edge_list));
    KATCH_STATUS("Graph has " << f_graph.get_n_nodes() << " nodes, " << f_graph.get_n_edges() << " edges.\n");
    S_Graph s_graph(f_graph,{hour*36000.0,hour*36000.0+3*36000.0});

    // output the free_flow cost
    std::vector<double> free_flow = s_graph.get_free_flow_vector();
    save_vector(main_name+"_" + to_string(hour) + ".fw_gtch_free_flow_4h",free_flow);
    KATCH_STATUS("Save free flow to " << main_name << "_" + to_string(hour) + ".fw_gtch_free_flow_4h"<<"\n");



    my_timer timer1 = my_timer();
    timer1.start();
    // resort the graph based on dfs ordering.
    s_graph.resort_graph(dfs_ordering);


    KATCH_STATUS("Building Min CPD ..." <<"\n");

    std::vector<katch::TchDijkstraQuery*> dijkstra;
    for (int i = 0; i < omp_get_max_threads(); i ++){
        dijkstra.push_back(new katch::TchDijkstraQuery(&s_graph));
    }
    unsigned number_of_nodes = s_graph.get_n_nodes();
    std::vector<unsigned>number_of_same_symbols(number_of_nodes,0);
    auto* cpd = new katch::CPD();
    {
        {
            katch::TchDijkstraQuery dij(&s_graph);
            my_timer t;
            t.start();
            dij.get_first_move(0);
            t.stop();
            double tots = t.elapsed_time_micro()*number_of_nodes/10 /1000000 / omp_get_max_threads();
            printf("Estimated sequential running time : %fmin\n", tots / 60.0);

        }

        printf("Using %d threads\n", omp_get_max_threads());
        std::vector<katch::CPD>thread_cpd(omp_get_max_threads());
        int progress = 0;

#pragma omp parallel
        {
            const int thread_count = omp_get_num_threads();
            const int thread_id = omp_get_thread_num();
//            const int begin_int  = 0;
            const int node_count = number_of_nodes;

            int node_begin = (node_count*thread_id) / thread_count ;
            int node_end = (node_count*(thread_id+1)) / thread_count ;
            katch::TchDijkstraQuery& thread_dij = *dijkstra[thread_id];

            for(int source_node=node_begin; source_node < node_end; ++source_node){
                std::vector<unsigned short> result = thread_dij.get_first_move(source_node);
                const std::vector<unsigned short> decode_results = _global_cpd->decode_row(source_node);
                std::vector<unsigned short> cpd_result (number_of_nodes*2);
                for(int i = 0; i < number_of_nodes;i++){
                    cpd_result[i*2] = result[i];
                    if(result[i] == decode_results[i]){
                        // 0XFD for same
                        cpd_result[i * 2 + 1 ] =0XFD;
                        number_of_same_symbols[source_node]++;
                    }else{
                        // 0XFC for not same
                        cpd_result[i * 2 + 1 ] =0XFC;
                    }
                }
//                cpd_result[source_node] = unordered_set<unsigned short>();
//                thread_cpd[thread_id].append_row(source_node,result);
                thread_cpd[thread_id].append_row_multiple_symbols(source_node,cpd_result);
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

        for(auto&x:thread_cpd){
            cpd->append_rows(x);
        }
    }
    timer1.stop();
    for (auto i: dijkstra) delete i;

    unsigned long total_same = 0;
    for (auto i : number_of_same_symbols ){
        total_same += i;
    }
    std::cout<<"Percentage of same symbol: " << (double)total_same/(dfs_ordering.size()*dfs_ordering.size()) <<std::endl;
    printf("Saving data to cpd.txt \n");
    printf("begin size: %d, entry size: %d\n", cpd->entry_count(), cpd->get_entry_size());


    std::string fname = main_name+"_"+ to_string(hour) +".fw_gtch_cpd_4h";
    FILE*f = fopen(fname.c_str(), "wb");
    cpd->save(f);
    fclose(f);
    cout << "done" << endl;
    KATCH_STATUS("Save CPD to " <<fname<<"\n");

    std::ofstream disFile(output_dir_name+"/fw_gtch_cpd_distribution_"+to_string(hour)+"_4h.csv");
    disFile<<"row_id,number_of_entries,number_of_same_symbols\n";
    for(int i = 0; i < number_of_nodes; i ++){
        disFile<<i <<","<<cpd->get_row_size(i)<<","<< number_of_same_symbols[i]<<"\n";
    }
    disFile.close();


    std::ofstream myFile(output_dir_name+"/fw_gtch_cpd_preprocessing_"+to_string(hour)+"_4h.csv");
    myFile<<"build_time,memory_cost\n";
    myFile<<std::fixed<<setprecision(8)<<timer1.elapsed_time_mins()<<","<<(double)cpd->get_entry_size()*4/1000000
          <<"\n";
    myFile.close();


    std::cout<< "Finished in: " << timer1.elapsed_time_mins()<<" Mins"<< std::endl;
    std::cout<< "Memory cost: " << (double)cpd->get_entry_size()*4/1000000<<" MB"<< std::endl;
}


void construct_hourly_gtch_cpd_4h(const std::string& tpgr_input_file_name, const std::string& btch_input_file_name, const string& output_dir_name){

    //take original ordering.
    string main_name = btch_input_file_name.substr(0,btch_input_file_name.find_last_of("."));
    std::vector<DEdgeInfo> d_edge_list = katch::tpgr_format::read_edges<DEdgeInfo>(tpgr_input_file_name);
    D_Graph d_graph(std::move(d_edge_list));
    std::vector<katch::NodeIterator > dfs_ordering = d_graph.generate_DFS_ordering();
    string mapper_file = main_name +"_global.fw_gtch_mapper_4h";
    vector<katch::NodeIterator > mapper = invert_permutation(dfs_ordering);
    save_vector(mapper_file,mapper);

    for(int i = 0; i < 24; i ++){
        // use one mapper only should be fine.
        if(i%4 ==0){
            construct_global_gtch_cpd_4h(btch_input_file_name, output_dir_name,dfs_ordering,i);
            delete(_global_cpd);
            _global_cpd = new katch::CPD();
            FILE*f = fopen((main_name+"_"+ to_string(i) +".fw_gtch_cpd_4h").c_str(), "r");
            _global_cpd->load(f);
        }else{
            construct_hourly_gtch_cpd_4h(btch_input_file_name, output_dir_name,dfs_ordering,i);
        }
    }
}














void construct_hourly_tch_cpd(const std::string& btch_input_file_name, const string& output_dir_name, const std::vector<katch::NodeIterator>& dfs_ordering, int hour) {

    string main_name = btch_input_file_name.substr(0,btch_input_file_name.find_last_of("."));
    std::vector<EdgeInfo> edge_list = katch::btch_format::read_edges<EdgeInfo>(btch_input_file_name,period);
    F_Graph f_graph(std::move(edge_list));
    KATCH_STATUS("Graph has " << f_graph.get_n_nodes() << " nodes, " << f_graph.get_n_edges() << " edges.\n");
    S_Graph s_graph(f_graph,{hour*36000.0,hour*36000.0+3*36000.0});

    // output the free_flow cost
    std::vector<double> free_flow = s_graph.get_free_flow_vector();
    save_vector(main_name+"_" + to_string(hour) + ".fw_tch_free_flow",free_flow);
    KATCH_STATUS("Save free flow to " << main_name << "_" + to_string(hour) + ".fw_tch_free_flow"<<"\n");



    my_timer timer1 = my_timer();
    timer1.start();
    // resort the graph based on dfs ordering.
    s_graph.resort_graph(dfs_ordering);


    KATCH_STATUS("Building Min CPD ..." <<"\n");

    std::vector<katch::TchDijkstraQuery*> dijkstra;
    for (int i = 0; i < omp_get_max_threads(); i ++){
        dijkstra.push_back(new katch::TchDijkstraQuery(&s_graph));
    }
    unsigned number_of_nodes = s_graph.get_n_nodes();
    std::vector<unsigned>number_of_same_symbols(number_of_nodes,0);
    auto* cpd = new katch::CPD();
    {
        {
            katch::TchDijkstraQuery dij(&s_graph);
            my_timer t;
            t.start();
            dij.get_first_move(0);
            t.stop();
            double tots = t.elapsed_time_micro()*number_of_nodes/10 /1000000 / omp_get_max_threads();
            printf("Estimated sequential running time : %fmin\n", tots / 60.0);

        }

        printf("Using %d threads\n", omp_get_max_threads());
        std::vector<katch::CPD>thread_cpd(omp_get_max_threads());
        int progress = 0;

#pragma omp parallel
        {
            const int thread_count = omp_get_num_threads();
            const int thread_id = omp_get_thread_num();
//            const int begin_int  = 0;
            const int node_count = number_of_nodes;

            int node_begin = (node_count*thread_id) / thread_count ;
            int node_end = (node_count*(thread_id+1)) / thread_count ;
            katch::TchDijkstraQuery& thread_dij = *dijkstra[thread_id];

            for(int source_node=node_begin; source_node < node_end; ++source_node){
                std::vector<unsigned short> result = thread_dij.get_first_move(source_node);
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

        for(auto&x:thread_cpd){
            cpd->append_rows(x);
        }
    }
    timer1.stop();
    for (auto i: dijkstra) delete i;

    unsigned long total_same = 0;
    for (auto i : number_of_same_symbols ){
        total_same += i;
    }
    std::cout<<"Percentage of same symbol: " << (double)total_same/(dfs_ordering.size()*dfs_ordering.size()) <<std::endl;
    printf("Saving data to cpd.txt \n");
    printf("begin size: %d, entry size: %d\n", cpd->entry_count(), cpd->get_entry_size());


    std::string fname = main_name+"_"+ to_string(hour) +".fw_tch_cpd";
    FILE*f = fopen(fname.c_str(), "wb");
    cpd->save(f);
    fclose(f);
    cout << "done" << endl;
    KATCH_STATUS("Save CPD to " <<fname<<"\n");

    std::ofstream disFile(output_dir_name+"/fw_tch_cpd_distribution_"+to_string(hour)+".csv");
    disFile<<"row_id,number_of_entries,number_of_same_symbols\n";
    for(int i = 0; i < number_of_nodes; i ++){
        disFile<<i <<","<<cpd->get_row_size(i)<<","<< number_of_same_symbols[i]<<"\n";
    }
    disFile.close();


    std::ofstream myFile(output_dir_name+"/fw_tch_cpd_preprocessing_"+to_string(hour)+".csv");
    myFile<<"build_time,memory_cost\n";
    myFile<<std::fixed<<setprecision(8)<<timer1.elapsed_time_mins()<<","<<(double)cpd->get_entry_size()*4/1000000
          <<"\n";
    myFile.close();


    std::cout<< "Finished in: " << timer1.elapsed_time_mins()<<" Mins"<< std::endl;
    std::cout<< "Memory cost: " << (double)cpd->get_entry_size()*4/1000000<<" MB"<< std::endl;
}


void construct_hourly_tch_cpd(const std::string& tpgr_input_file_name, const std::string& btch_input_file_name, const string& output_dir_name){
    //take original ordering.
    string main_name = btch_input_file_name.substr(0,btch_input_file_name.find_last_of("."));
    std::vector<DEdgeInfo> d_edge_list = katch::tpgr_format::read_edges<DEdgeInfo>(tpgr_input_file_name);
    D_Graph d_graph(std::move(d_edge_list));
    std::vector<katch::NodeIterator > dfs_ordering = d_graph.generate_DFS_ordering();
    string mapper_file = main_name +".fw_tch_mapper";
    vector<katch::NodeIterator > mapper = invert_permutation(dfs_ordering);
    save_vector(mapper_file,mapper);
    KATCH_STATUS("Save mapper to " <<mapper_file<<"\n");

    for(int i = 0; i < 24; i ++){
        // use one mapper only should be fine.
        construct_hourly_tch_cpd(btch_input_file_name, output_dir_name,dfs_ordering,i);
    }
}



void check_graph_difference (const std::string& tpgr_input_file_name, const std::string& btch_input_file_name, const string& output_dir_name){

    string main_name = btch_input_file_name.substr(0,btch_input_file_name.find_last_of("."));
    std::vector<EdgeInfo> edge_list = katch::btch_format::read_edges<EdgeInfo>(btch_input_file_name,period);
    F_Graph f_graph(std::move(edge_list));
    KATCH_STATUS("Graph has " << f_graph.get_n_nodes() << " nodes, " << f_graph.get_n_edges() << " edges.\n");
    S_Graph s_graph(f_graph,F_Graph::min_graph);

    // output the free_flow cost
    std::vector<double> global_free_flow = s_graph.get_free_flow_vector();

    unsigned num_td_edges = 0;
    for(int i =0 ; i < global_free_flow.size(); i ++){
        if(!f_graph.get_ttf((katch::EdgeIterator) i).is_constant()){
            num_td_edges ++;
        }
    }
    std::cout<<" Number_of_TD_edges: "<<(double)num_td_edges/ global_free_flow.size()<<std::endl;
    for(int hour  = 0 ; hour < 24;  hour++){

        S_Graph s_graph_cur(f_graph,{hour*36000.0,hour*36000.0+3*36000.0});
        std::vector<double> cur_free_flow = s_graph_cur.get_free_flow_vector();
        unsigned edges_have_different_cost = 0;
        for(int i =0 ; i < cur_free_flow.size(); i ++){
            if(cur_free_flow[i] != global_free_flow[i]){
                edges_have_different_cost ++;
            }
        }
        std::cout<<"Hour: "<< hour<< " Number of edges has different cost: "<<(double)edges_have_different_cost / cur_free_flow.size()<<std::endl;
    }



}

void check_raw_first_move (const std::string& tpgr_input_file_name, const std::string& btch_input_file_name, const string& output_dir_name) {
    _global_cpd = new katch::CPD();

    string main_name = btch_input_file_name.substr(0,btch_input_file_name.find_last_of("."));
    FILE*f = fopen((main_name +"_global.fw_gtch_cpd").c_str(), "r");
    _global_cpd->load(f);
    std::vector<katch::CPD* >gtch_cpd_vec =  std::vector<katch::CPD* >(24);
    std::vector<katch::CPD* >ori_cpd_vec =  std::vector<katch::CPD* >(24);
    for(int i = 0;  i < 24; i ++){
        katch::CPD* _cur_gtch_cpd = new katch::CPD();
        katch::CPD* _cur_ori_cpd = new katch::CPD();
        string gtchcpd_name =  main_name+"_"+ to_string(i) +".fw_gtch_cpd";
        string ori_cpd_name =  main_name+"_"+ to_string(i) +".fw_tch_cpd";
        FILE* gtchcpd_f = fopen( gtchcpd_name.c_str(), "r");
        FILE* ori_cpd_f = fopen(ori_cpd_name .c_str(), "r");
        _cur_gtch_cpd ->load(gtchcpd_f);
        _cur_ori_cpd->load(ori_cpd_f);
        gtch_cpd_vec[i] = _cur_gtch_cpd;
        ori_cpd_vec[i] = _cur_ori_cpd;
        std::cout<<"Finished loading "<<i<<" CPD"<<std::endl;
    }

    for(int i = 0 ; i < _global_cpd->node_count(); i++){
        const std::vector<unsigned short> decode_results = _global_cpd->decode_row(i);

        for(int j = 0;  j < 24; j ++){
            std::vector<unsigned short> gtch_cpd_decode_results = gtch_cpd_vec[j]->decode_row(i);
            std::vector<unsigned short> ori_cpd_decode_results = ori_cpd_vec[j]->decode_row(i);
            for(int k = 0 ; k < _global_cpd->node_count(); k++){
                if(gtch_cpd_decode_results[k] == 0XFD){
                    gtch_cpd_decode_results[k] = decode_results[k];
                }
                if(gtch_cpd_decode_results[k] != ori_cpd_decode_results[k]){
                    std::cout<<"Symbol not matched" << std::endl;
                    return;
                }
            }

        }
        std::cout<<"Finished checking raw "<<i<<std::endl;
    }
}




void check_compression (const std::string& tpgr_input_file_name, const std::string& btch_input_file_name, const string& output_dir_name) {
    string main_name = btch_input_file_name.substr(0,btch_input_file_name.find_last_of("."));
    std::vector<katch::CPD* >gtch_cpd_vec =  std::vector<katch::CPD* >(24);
    for(int i = 0;  i < 24; i ++){
        katch::CPD* _cur_gtch_cpd = new katch::CPD();
        string gtchcpd_name =  main_name+"_"+ to_string(i) +".fw_gtch_cpd";
        FILE* gtchcpd_f = fopen( gtchcpd_name.c_str(), "r");
        _cur_gtch_cpd ->load(gtchcpd_f);
        gtch_cpd_vec[i] = _cur_gtch_cpd;
        std::cout<<"Finished loading "<<i<<" CPD"<<std::endl;
    }

    for(int i = 0 ; i < gtch_cpd_vec[0]->node_count(); i++){
        for(int j = 0;  j < 24; j ++){
            std::vector<unsigned short> gtch_cpd_decode_results = gtch_cpd_vec[j]->get_compressed_symbol(i);
            if(gtch_cpd_decode_results.size() != 1){
                for(int k = 0; k < gtch_cpd_decode_results.size()-2; k++){
                    if(gtch_cpd_decode_results[k] == gtch_cpd_decode_results[k+1]){
                        std::cout<<"Compression errors"<<std::endl;
                    }
                }
            }



        }

    }
}


void print_raw_first_move (const std::string& tpgr_input_file_name, const std::string& btch_input_file_name, const string& output_dir_name, int row_id) {
    _global_cpd = new katch::CPD();

    string main_name = btch_input_file_name.substr(0,btch_input_file_name.find_last_of("."));
    FILE*f = fopen((main_name +"_global.fw_gtch_cpd").c_str(), "r");
    _global_cpd->load(f);
    std::vector<katch::CPD* >gtch_cpd_vec =  std::vector<katch::CPD* >(24);
    std::vector<katch::CPD* >ori_cpd_vec =  std::vector<katch::CPD* >(24);
    for(int i = 0;  i < 24; i ++){
        katch::CPD* _cur_gtch_cpd = new katch::CPD();
        katch::CPD* _cur_ori_cpd = new katch::CPD();
        string gtchcpd_name =  main_name+"_"+ to_string(i) +".fw_gtch_cpd";
        string ori_cpd_name =  main_name+"_"+ to_string(i) +".fw_tch_cpd";
        FILE* gtchcpd_f = fopen( gtchcpd_name.c_str(), "r");
        FILE* ori_cpd_f = fopen(ori_cpd_name .c_str(), "r");
        _cur_gtch_cpd ->load(gtchcpd_f);
        _cur_ori_cpd->load(ori_cpd_f);
        gtch_cpd_vec[i] = _cur_gtch_cpd;
        ori_cpd_vec[i] = _cur_ori_cpd;
        std::cout<<"Finished loading "<<i<<" CPD"<<std::endl;
    }

    std::vector<std::vector<unsigned short>> raw_first_moves =     std::vector<std::vector<unsigned short>>(73);
    raw_first_moves[0] =  _global_cpd->decode_row(row_id);
    for(int j = 0;  j < 24; j ++){
        raw_first_moves[j+1] = ori_cpd_vec[j]->decode_row(row_id);
    }

    for(int j = 0;  j < 24; j ++){
        raw_first_moves[j+25] = gtch_cpd_vec[j]->decode_row(row_id);
    }

    for(int j = 0;  j < 24; j ++){
        std::vector<unsigned short> global = _global_cpd->decode_row(row_id);
        std::vector<unsigned short> cur = ori_cpd_vec[j]->decode_row(row_id);
        for(int k = 0; k < _global_cpd->node_count(); k++ ){
            if(global[k] == cur[k]){
                cur[k] = 0XFD;
            }
        }

        raw_first_moves[j+49] = cur;
    }


    std::ofstream myFile(output_dir_name+"/raw_first_move"+to_string(row_id)+".csv");
    myFile<<"0";
    for(int i = 1; i <73; i ++){
        myFile<<","+to_string(i);
    }
    myFile<<"\n";

    for(int i = 0; i <_global_cpd->node_count(); i ++){
        myFile<<to_string(raw_first_moves[0][i]);
        for(int j= 1; j <73; j ++){
            myFile<<","+to_string(raw_first_moves[j][i]);
        }
        myFile<<"\n";
    }
    myFile.close();

}


void print_raw_first_move_with_wkt (const std::string& tpgr_input_file_name, const std::string& btch_input_file_name, const string& output_dir_name, int row_id) {

    string main_name = btch_input_file_name.substr(0,btch_input_file_name.find_last_of("."));
    std::vector<EdgeInfo> edge_list = katch::btch_format::read_edges<EdgeInfo>(btch_input_file_name,period);
    F_Graph f_graph(std::move(edge_list));
    KATCH_STATUS("Graph has " << f_graph.get_n_nodes() << " nodes, " << f_graph.get_n_edges() << " edges.\n");

    f_graph.load_coordinate(main_name +".coordinate");
    _global_cpd = new katch::CPD();
    FILE*f = fopen((main_name +"_global.fw_gtch_cpd").c_str(), "r");
    _global_cpd->load(f);
    std::vector<katch::CPD* >gtch_cpd_vec =  std::vector<katch::CPD* >(24);
    std::vector<katch::CPD* >ori_cpd_vec =  std::vector<katch::CPD* >(24);
    for(int i = 0;  i < 24; i ++){
        katch::CPD* _cur_gtch_cpd = new katch::CPD();
        katch::CPD* _cur_ori_cpd = new katch::CPD();
        string gtchcpd_name =  main_name+"_"+ to_string(i) +".fw_gtch_cpd";
        string ori_cpd_name =  main_name+"_"+ to_string(i) +".fw_tch_cpd";
        FILE* gtchcpd_f = fopen( gtchcpd_name.c_str(), "r");
        FILE* ori_cpd_f = fopen(ori_cpd_name .c_str(), "r");
        _cur_gtch_cpd ->load(gtchcpd_f);
        _cur_ori_cpd->load(ori_cpd_f);
        gtch_cpd_vec[i] = _cur_gtch_cpd;
        ori_cpd_vec[i] = _cur_ori_cpd;
        std::cout<<"Finished loading "<<i<<" CPD"<<std::endl;
    }
    //convert to inverse mapper
    std::vector<katch::NodeIterator > mapper = load_vector<katch::NodeIterator >(main_name+"_global.fw_gtch_mapper");
    mapper = invert_permutation(mapper);

    vector<pair<unsigned short, unsigned int>> row_dis = _global_cpd->get_row_distribution(row_id);
    vector<unsigned short> gobal_decode_row = _global_cpd->decode_row(row_id);
    for(int i = 0 ; i < 5; i++ ){
        unsigned short symbol = row_dis[i].first;
        vector<unsigned short> ori_decode_row = ori_cpd_vec[0]->decode_row(row_id);
        vector<unsigned >ori_id;
        vector<unsigned >global_id;
        for(int j = 0; j < ori_decode_row.size(); j ++){
            if(ori_decode_row[j] == symbol){
                ori_id.push_back(mapper[j]);
            }
            if(gobal_decode_row[j] == symbol){
                global_id.push_back(mapper[j]);
            }
        }
        f_graph.export_wkt(ori_id,output_dir_name +"/ori_top_"+ to_string(i) +".csv");
        f_graph.export_wkt(global_id,output_dir_name +"/global_top_"+ to_string(i) +".csv");
    }


    vector<unsigned >target_id;
    target_id.push_back(mapper[row_id]);
    f_graph.export_wkt(target_id,output_dir_name +"/target.csv");



}



void check_node_level (const std::string& tpgr_input_file_name, const std::string& btch_input_file_name, const string& output_dir_name){

    string main_name = btch_input_file_name.substr(0,btch_input_file_name.find_last_of("."));
    std::vector<EdgeInfo> edge_list = katch::btch_format::read_edges<EdgeInfo>(btch_input_file_name,period);
    F_Graph f_graph(std::move(edge_list));
    auto ordering = katch::btch_format::read_ordering(btch_input_file_name,period);
    f_graph.test_ordering(ordering);
    std::vector<katch::NodeIterator > mapper = load_vector<katch::NodeIterator >(main_name+"_global.fw_gtch_mapper");
    vector<unsigned> converted_ordering (ordering.size());

    for(int i = 0; i < ordering.size(); i++){
        converted_ordering[mapper[i]]     =        ordering[i];
    }

    std::ofstream myFile(output_dir_name+"/TCH_ordering.csv");
    myFile<<"Ordering\n";
    for(auto l : converted_ordering){
        myFile<< l <<"\n";
    }

//    f_graph.check_graph();
//    std::cout<<"finished"<<std::endl;
//    std::vector<int> node_level = std::vector<int>(f_graph.get_n_nodes(),-1);
//    f_graph.mark_level_down2(node_level, 32213);
//    std::cout<<"finished"<<std::endl;
////
//
//    vector<int> node_level = f_graph.get_node_level();
//
//    std::map<int, unsigned> level_map;
//    int max_level = 0;
//    for(const auto& level : node_level){
//        if(level_map.find(level) != level_map.end()){
//            level_map[level] ++;
//        }else{
//            level_map[level] = 0;
//        }
//        max_level = std::max(max_level,level);
//    }
//    for(int i = 0; i < max_level; i ++){
//        std::cout<<"Level :" << to_string(i) <<" #Nodes: "<< level_map[i] <<std::endl;
//    }

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

    // build original CPD
    construct_min_cpd(tpgr_input_file_name, btch_input_file_name, output_dir_name);

    // build TCPD
    construct_min_tcpd(tpgr_input_file_name, btch_input_file_name, output_dir_name);
    // build TCPD for STCH
    construct_stch_tcpd<&STCH_period>(tpgr_input_file_name, btch_input_file_name, output_dir_name);
    // build TCPD for MTCH
    construct_mtch_tcpd<&one_MTCH_period>(1,1,tpgr_input_file_name, btch_input_file_name, output_dir_name);
    construct_mtch_tcpd<&two_MTCH_period>(1,2,tpgr_input_file_name, btch_input_file_name, output_dir_name);
    construct_mtch_tcpd<&four_MTCH_period>(1,4,tpgr_input_file_name, btch_input_file_name, output_dir_name);

}