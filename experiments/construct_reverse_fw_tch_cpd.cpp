//
// Created by Bojie Shen on 22/10/21.
//

#include <iostream>
#include <stdlib.h>
#include <string>
#include <util/my_timer.h>
#include <unordered_set>
#include <iomanip>

#include "util/misc.h"
#include "io/btch_format.h"
#include "io/edge_info.h"
#include "io/vec_io.h"
#include "datastr/graph/forward_search_graph.h"
#include "cpd_prepro/tch_dijkstra_query.h"
#include "omp.h"
#include "datastr/cpd/rev_cpd.h"
#include <map>
#include "cpd_prepro/reverse_tch_updown_query.h"
static constexpr double period = 864000.0;
static constexpr double h_period_5 = 185000.0;
using S_Graph = katch::StaticForwardSearchGraph;
using F_Graph = katch::ForwardSearchGraph<&period>;
using EdgeInfo = katch::EdgeInfo<katch::ForwardSearchGraph<&period>::TTF>;


void construct_min_reverse_tch_cpd(const std::string& tpgr_input_file_name, const std::string& btch_input_file_name, const string& output_dir_name) {


    string main_name = btch_input_file_name.substr(0,btch_input_file_name.find_last_of("."));

    std::vector<int> target_list =  load_vector<int>(main_name+".target");

    std::vector<EdgeInfo> edge_list = katch::btch_format::read_edges<EdgeInfo>(btch_input_file_name,period);
    F_Graph f_graph(std::move(edge_list));
    KATCH_STATUS("Graph has " << f_graph.get_n_nodes() << " nodes, " << f_graph.get_n_edges() << " edges.\n");
    S_Graph in_graph(f_graph,F_Graph::min_graph);
    in_graph.convert_to_income_graph();


    my_timer timer1 = my_timer();
    timer1.start();


    KATCH_STATUS("Building Min CPD ..." <<"\n");
    std::vector<katch::ReverseTchUpDownQuery*> dijkstra;
    for (int i = 0; i < omp_get_max_threads(); i ++){
        dijkstra.push_back(new katch::ReverseTchUpDownQuery(&in_graph));
    }

    unsigned number_of_nodes = target_list.size();
    auto* cpd = new katch::REV_CPD();
    {
        {
            katch::ReverseTchUpDownQuery dij(&in_graph);
            my_timer t;
            t.start();
            dij.get_first_move(0);
            t.stop();
            double tots = t.elapsed_time_micro()*number_of_nodes / 1000000/ omp_get_max_threads();
            printf("Estimated sequential running time : %fmin\n", tots / 60.0);

        }

        printf("Using %d threads\n", omp_get_max_threads());
        std::vector<katch::REV_CPD>thread_cpd(omp_get_max_threads());

        int progress = 0;

#pragma omp parallel
        {
            const int thread_count = omp_get_num_threads();
            const int thread_id = omp_get_thread_num();
//            const int begin_int  = 0;
            const int node_count = number_of_nodes;

            int node_begin = (node_count*thread_id) / thread_count ;
            int node_end = (node_count*(thread_id+1)) / thread_count ;
            katch::ReverseTchUpDownQuery& thread_dij = *dijkstra[thread_id];

            for(int source_node=node_begin; source_node < node_end; ++source_node){
                const std::vector<unsigned char>& result = thread_dij.get_first_move(target_list[source_node]);
                thread_cpd[thread_id].append_row(result);
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


    std::string fname = main_name+"_min.rev_fw_tch_cpd";
    FILE*f = fopen(fname.c_str(), "wb");
    cpd->save(f);
    fclose(f);
    cout << "done" << endl;
    KATCH_STATUS("Save CPD to " <<fname<<"\n");



    std::ofstream myFile(output_dir_name+"/rev_fw_tch_cpd_preprocessing_min.csv");
    myFile<<"build_time,memory_cost\n";
    myFile<<std::fixed<<setprecision(8)<<timer1.elapsed_time_mins()/target_list.size()*in_graph.get_n_nodes()<<","<<(double)cpd->get_entry_size()/1000000/target_list.size()*in_graph.get_n_nodes()
          <<"\n";
    myFile.close();

    std::cout<<target_list.size()<<std::endl;
    std::cout<< "Finished in: " << timer1.elapsed_time_mins()/target_list.size()*in_graph.get_n_nodes()<<" Mins"<< std::endl;
    std::cout<< "Memory cost: " << (double)cpd->get_entry_size()/1000000/target_list.size()*in_graph.get_n_nodes()<<" MB"<< std::endl;
}





void construct_min_reverse_tch_cpd_eval_time(const std::string& tpgr_input_file_name, const std::string& btch_input_file_name, const string& output_dir_name) {


    string main_name = btch_input_file_name.substr(0,btch_input_file_name.find_last_of("."));


    std::vector<EdgeInfo> edge_list = katch::btch_format::read_edges<EdgeInfo>(btch_input_file_name,period);
    F_Graph f_graph(std::move(edge_list));
    KATCH_STATUS("Graph has " << f_graph.get_n_nodes() << " nodes, " << f_graph.get_n_edges() << " edges.\n");
    S_Graph in_graph(f_graph,F_Graph::min_graph);
    in_graph.convert_to_income_graph();


    my_timer timer1 = my_timer();
    timer1.start();


    KATCH_STATUS("Building Min CPD ..." <<"\n");
    std::vector<katch::ReverseTchUpDownQuery*> dijkstra;
    for (int i = 0; i < omp_get_max_threads(); i ++){
        dijkstra.push_back(new katch::ReverseTchUpDownQuery(&in_graph));
    }

    unsigned number_of_nodes = in_graph.get_n_nodes();
    auto* cpd = new katch::REV_CPD();
    {
        {
            katch::ReverseTchUpDownQuery dij(&in_graph);
            my_timer t;
            t.start();
            dij.get_first_move(0);
            t.stop();
            double tots = t.elapsed_time_micro()*number_of_nodes / 1000000/ omp_get_max_threads();
            printf("Estimated sequential running time : %fmin\n", tots / 60.0);

        }

        printf("Using %d threads\n", omp_get_max_threads());
        std::vector<katch::REV_CPD>thread_cpd(omp_get_max_threads());

        int progress = 0;

#pragma omp parallel
        {
            const int thread_count = omp_get_num_threads();
            const int thread_id = omp_get_thread_num();
//            const int begin_int  = 0;
            const int node_count = number_of_nodes;

            int node_begin = (node_count*thread_id) / thread_count ;
            int node_end = (node_count*(thread_id+1)) / thread_count ;
            katch::ReverseTchUpDownQuery& thread_dij = *dijkstra[thread_id];

            for(int source_node=node_begin; source_node < node_end; ++source_node){
                const std::vector<unsigned char>& result = thread_dij.get_first_move(source_node);
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

    std::ofstream myFile(output_dir_name+"/rev_fw_tch_cpd_time_only.csv");
    myFile<<"build_time,memory_cost\n";
    myFile<<std::fixed<<setprecision(8)<<timer1.elapsed_time_mins()<<","<<0
          <<"\n";
    myFile.close();

    std::cout<< "Finished in: " << timer1.elapsed_time_mins()<<" Mins"<< std::endl;
}





template<const double* h_period>
void construct_reverse_hourly_min_htch_cpd(const std::string& btch_input_file_name, const string& output_dir_name,
                                            const std::vector<int>& target_list , int hour) {

    string main_name = btch_input_file_name.substr(0,btch_input_file_name.find_last_of("."));

    using H_Graph = katch::ForwardSearchGraph<h_period>;
    using HEdgeInfo = katch::EdgeInfo<typename H_Graph::TTF>;
    std::vector<HEdgeInfo> edge_list = katch::btch_format::read_edges<HEdgeInfo>(btch_input_file_name,*h_period);
    H_Graph h_graph(std::move(edge_list));
    KATCH_STATUS("Graph has " << h_graph.get_n_nodes() << " nodes, " << h_graph.get_n_edges() << " edges.\n");
    S_Graph s_graph(h_graph,H_Graph::min_graph);
    s_graph.convert_to_income_graph();


    my_timer timer1 = my_timer();
    timer1.start();


    KATCH_STATUS("Building Min CPD ..." <<"\n");

    std::vector<katch::ReverseTchUpDownQuery*> dijkstra;
    for (int i = 0; i < omp_get_max_threads(); i ++){
        dijkstra.push_back(new katch::ReverseTchUpDownQuery(&s_graph));
    }

    unsigned number_of_nodes = target_list.size();
    auto* cpd = new katch::REV_CPD();
    {
        {
            katch::ReverseTchUpDownQuery dij(&s_graph);
            my_timer t;
            t.start();
            dij.get_first_move(0);
            t.stop();
            double tots = t.elapsed_time_micro()*number_of_nodes/10 /1000000 / omp_get_max_threads();
            printf("Estimated sequential running time : %fmin\n", tots / 60.0);

        }

        printf("Using %d threads\n", omp_get_max_threads());
        std::vector<katch::REV_CPD>thread_cpd(omp_get_max_threads());

        int progress = 0;

#pragma omp parallel
        {
            const int thread_count = omp_get_num_threads();
            const int thread_id = omp_get_thread_num();
//            const int begin_int  = 0;
            const int node_count = number_of_nodes;

            int node_begin = (node_count*thread_id) / thread_count ;
            int node_end = (node_count*(thread_id+1)) / thread_count ;
            katch::ReverseTchUpDownQuery& thread_dij = *dijkstra[thread_id];

            for(int source_node=node_begin; source_node < node_end; ++source_node){
                const std::vector<unsigned char>& result = thread_dij.get_first_move( target_list[source_node]);
                thread_cpd[thread_id].append_row(result);
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


    std::string fname = main_name+"_min.rev_fw_htch_cpd";
    FILE*f = fopen(fname.c_str(), "wb");
    cpd->save(f);
    fclose(f);
    cout << "done" << endl;
    KATCH_STATUS("Save CPD to " <<fname<<"\n");

    std::ofstream myFile(output_dir_name+"/rev_fw_htch_cpd_preprocessing_min_"+to_string(hour)+".csv");
    myFile<<"build_time,memory_cost\n";
    myFile<<std::fixed<<setprecision(8)<<timer1.elapsed_time_mins()/target_list.size()*s_graph.get_n_nodes()<<","<<(double)cpd->get_entry_size()/1000000/target_list.size()*s_graph.get_n_nodes()
          <<"\n";
    myFile.close();

    std::cout<<target_list.size()<<std::endl;
    std::cout<< "Finished in: " << timer1.elapsed_time_mins()/target_list.size()*s_graph.get_n_nodes()<<" Mins"<< std::endl;
    std::cout<< "Memory cost: " << (double)cpd->get_entry_size()/1000000/target_list.size()*s_graph.get_n_nodes()<<" MB"<< std::endl;
    delete cpd;

}


void combine_hourly_reverse_hourly_cpd(const std::string& main_name){
    unsigned  long row_size = 0 ;
    size_t  total_size = 0;
    for(int i = 0; i < 24; i ++){
        string hourly_cpd_name = main_name + "_"+to_string(i) + "_min.rev_fw_htch_cpd";
        FILE*f = fopen(hourly_cpd_name.c_str(), "r");
        if(std::fread(& row_size, sizeof( row_size), 1, f) != 1)
            throw std::runtime_error("std::fread failed");
        size_t s;
        if(std::fread(&s, sizeof(s), 1, f) != 1)
            throw std::runtime_error("std::fread failed");
        total_size += s ;
        fclose(f);
    }

    string output_cpd_name = main_name + "_min.rev_fw_htch_cpd";
    FILE*output_file = fopen(output_cpd_name.c_str(), "wb");
    if(std::fwrite(&row_size, sizeof(row_size), 1, output_file) != 1)
        throw std::runtime_error("std::fwrite failed");
    if(std::fwrite(&total_size, sizeof(total_size), 1, output_file) != 1)
        throw std::runtime_error("std::fwrite failed");


    for(int i = 0; i < 24; i ++){
        std::cout<<"finished writing file: "<<  to_string(i)<< std::endl;
        string hourly_cpd_name = main_name + "_"+to_string(i) + "_min.rev_fw_htch_cpd";

        FILE*f = fopen(hourly_cpd_name.c_str(), "r");
        if(std::fread(& row_size, sizeof( row_size), 1, f) != 1)
            throw std::runtime_error("std::fread failed");
        size_t s;
        if(std::fread(&s, sizeof(s), 1, f) != 1)
            throw std::runtime_error("std::fread failed");

        vector<unsigned char> entry(s);
        if((size_t)std::fread(& entry[0], sizeof(unsigned char), s, f) != s)
            throw std::runtime_error("std::fread failed");

        //write to file
        if(std::fwrite(&entry[0],  sizeof(unsigned char),s,output_file) != s)
            throw std::runtime_error("std::fwrite failed");
        fclose(f);

        std::cout<<"finished writing file: "<<  to_string(i)<< std::endl;
        int status = remove( hourly_cpd_name.c_str() );
        if(status==0)
            cout<<"\nFile Deleted Successfully!";
        else
            cout<<"\nError Occurred!";

    }
    fclose(output_file);

    string map_name = main_name.substr(main_name.find_last_of("/") + 1 );

    // combine preprocessing results.

    double total_build_time = 0;
    double total_memory = 0;
    for(int i = 0; i < 24; i ++){
        string preprocessing_file = "dataset/"+ map_name + "/results/preprocessing/rev_fw_htch_cpd_preprocessing_min_" + to_string(i)+".csv";
        std::ifstream myFile(preprocessing_file);
        std::string line;
        double val;
        int line_id = 0;
        vector<double> value_list;
        while(std::getline(myFile, line))
        {
            if(line_id == 0 ){
                line_id ++;
                continue;
            }
            // Create a stringstream of the current line
            std::stringstream ss(line);

            // Keep track of the current column index
            int colIdx = 0;

            // Extract each integer
            while(ss >> val){

                // Add the current integer to the 'colIdx' column's values vector
                value_list.push_back(val);
                // If the next token is a comma, ignore it and move on
                if(ss.peek() == ',') ss.ignore();

                // Increment the column index
                colIdx++;
            }
            line_id ++;
        }
        total_build_time += value_list[0];
        total_memory  += value_list[1];
        myFile.close();
    }


    std::ofstream myFile("dataset/"+ map_name + "/results/preprocessing/rev_fw_htch_cpd_preprocessing_min.csv");
    myFile<<"build_time,memory_cost\n";
    myFile<<std::fixed<<setprecision(8)<<total_build_time<<","<<total_memory
          <<"\n";
    myFile.close();

}

template<const double* h_period>
void construct_reverse_hourly_htch_cpd(const std::string& tpgr_input_file_name, const std::string& btch_input_file_name, const string& output_dir_name){
    string main_name = btch_input_file_name.substr(0,btch_input_file_name.find_last_of("."));
    std::vector<int> target_list =  load_vector<int>(main_name+".target");
    for(int i = 0; i < 24; i ++){
        // use one mapper only should be fine.
        construct_reverse_hourly_min_htch_cpd<h_period>(main_name+"_"+to_string(i)+".btch", output_dir_name, target_list, i);
    }
    combine_hourly_reverse_hourly_cpd(main_name);
}









template<const double* h_period>
void construct_reverse_hourly_min_htch_cpd_eval_time(const std::string& btch_input_file_name, const string& output_dir_name,int hour) {

    string main_name = btch_input_file_name.substr(0,btch_input_file_name.find_last_of("."));

    using H_Graph = katch::ForwardSearchGraph<h_period>;
    using HEdgeInfo = katch::EdgeInfo<typename H_Graph::TTF>;
    std::vector<HEdgeInfo> edge_list = katch::btch_format::read_edges<HEdgeInfo>(btch_input_file_name,*h_period);
    H_Graph h_graph(std::move(edge_list));
    KATCH_STATUS("Graph has " << h_graph.get_n_nodes() << " nodes, " << h_graph.get_n_edges() << " edges.\n");
    S_Graph s_graph(h_graph,H_Graph::min_graph);
    s_graph.convert_to_income_graph();


    my_timer timer1 = my_timer();
    timer1.start();


    KATCH_STATUS("Building Min CPD ..." <<"\n");

    std::vector<katch::ReverseTchUpDownQuery*> dijkstra;
    for (int i = 0; i < omp_get_max_threads(); i ++){
        dijkstra.push_back(new katch::ReverseTchUpDownQuery(&s_graph));
    }

    unsigned number_of_nodes = s_graph.get_n_nodes();
    auto* cpd = new katch::REV_CPD();
    {
        {
            katch::ReverseTchUpDownQuery dij(&s_graph);
            my_timer t;
            t.start();
            dij.get_first_move(0);
            t.stop();
            double tots = t.elapsed_time_micro()*number_of_nodes/10 /1000000 / omp_get_max_threads();
            printf("Estimated sequential running time : %fmin\n", tots / 60.0);

        }

        printf("Using %d threads\n", omp_get_max_threads());
        std::vector<katch::REV_CPD>thread_cpd(omp_get_max_threads());

        int progress = 0;

#pragma omp parallel
        {
            const int thread_count = omp_get_num_threads();
            const int thread_id = omp_get_thread_num();
//            const int begin_int  = 0;
            const int node_count = number_of_nodes;

            int node_begin = (node_count*thread_id) / thread_count ;
            int node_end = (node_count*(thread_id+1)) / thread_count ;
            katch::ReverseTchUpDownQuery& thread_dij = *dijkstra[thread_id];

            for(int source_node=node_begin; source_node < node_end; ++source_node){
                const std::vector<unsigned char>& result = thread_dij.get_first_move( source_node);
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

    std::ofstream myFile(output_dir_name+"/rev_fw_htch_cpd_time_only_"+to_string(hour)+".csv");
    myFile<<"build_time,memory_cost\n";
    myFile<<std::fixed<<setprecision(8)<<timer1.elapsed_time_mins()<<","<<0
          <<"\n";
    myFile.close();

    std::cout<< "Finished in: " << timer1.elapsed_time_mins()<<" Mins"<< std::endl;
    delete cpd;

}


template<const double* h_period>
void construct_reverse_hourly_htch_cpd_eval_time(const std::string& tpgr_input_file_name, const std::string& btch_input_file_name, const string& output_dir_name){
    string main_name = btch_input_file_name.substr(0,btch_input_file_name.find_last_of("."));
    std::vector<int> target_list =  load_vector<int>(main_name+".target");
    for(int i = 0; i < 24; i ++){
        // use one mapper only should be fine.
        construct_reverse_hourly_min_htch_cpd_eval_time<h_period>(main_name+"_"+to_string(i)+".btch", output_dir_name,i);
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


//    construct_min_reverse_tch_cpd(tpgr_input_file_name, btch_input_file_name, output_dir_name);

//    construct_reverse_hourly_htch_cpd<&h_period_5>(tpgr_input_file_name, btch_input_file_name, output_dir_name);


    construct_min_reverse_tch_cpd_eval_time(tpgr_input_file_name, btch_input_file_name, output_dir_name);
//    construct_reverse_hourly_htch_cpd_eval_time<&h_period_5>(tpgr_input_file_name, btch_input_file_name, output_dir_name);
}