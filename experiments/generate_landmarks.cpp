//
// Created by Bojie Shen on 21/7/21.
//

#include <iostream>
#include <stdlib.h>
#include <string>
#include <io/edge_info.h>
#include <io/tpgr_format.h>
#include "datastr/graph/dynamic_search_graph.h"
#include "cpd_prepro/st_dijkstra_query.h"
#include <iomanip>
#define INF std::numeric_limits<double>::max()
using S_Graph = katch::StDijkstraQuery::Graph;

struct landmark{
    unsigned id;
    std::vector<double> out_distance;
    std::vector<double> in_distance;
};
void generate_landmarks(S_Graph & out_graph,S_Graph & in_graph, int n_landmark, std::string output_file ){
    katch::StDijkstraQuery out_dij(&out_graph);
    katch::StDijkstraQuery in_dij(&in_graph);
    std::vector<landmark> l_list;
    std::cout<<"Building landmark ...."<<std::endl;
    unsigned random =0;
    std::vector<double> d_list =out_dij.get_distance_to_all_nodes(random);
    unsigned max_id =0;
    double  max = 0;
    unsigned cur_id = 0;
    for(double d: d_list){
        if(d > max && d != INF){
            max = d;
            max_id = cur_id;
        }
        cur_id ++;
    }
    // start with the farthest node
    for(int i = 0 ; i < n_landmark; i++){
        std::cout<<"Selected node: " << max_id << std::endl;
        d_list = out_dij.get_distance_to_all_nodes(max_id);
        std::vector<double> in_d_list =in_dij.get_distance_to_all_nodes(max_id);
        landmark l  = landmark{max_id,d_list,in_d_list};
        l_list.push_back(l);

        max_id = 0;
        max = 0;
        for(int j = 0; j < out_graph.get_n_nodes(); j ++){
            // select the node that have maximal sum distance ;
            double cur_distance =0;
            for( const landmark& ll : l_list){
                if(ll.id == j){
                    cur_distance = 0;
                    break;
                }
                if(ll.out_distance[j] != INF){
                    cur_distance = cur_distance + ll.out_distance[j];
                }
            }
            if(cur_distance > max) {
                max= cur_distance;
                max_id = j;
            }
        }
    }
    std::vector<double > merged_l_list = std::vector<double>(n_landmark*out_graph.get_n_nodes()*2);
    for(int i = 0; i< out_graph.get_n_nodes(); i ++){
        for(int j = 0; j < l_list.size();j++){
            merged_l_list[i*n_landmark*2+j*2] =  l_list[j].out_distance[i];
            merged_l_list[i*n_landmark*2+j*2+1] =  l_list[j].in_distance[i];
        }
    }
    save_vector(output_file, merged_l_list );

}




int main(int argc, char **argv) {
    const char* binary_name = argv[0];

    if ( argc != 3 )
    {
        std::cerr
                << std::endl
                << "USAGE: " << binary_name
                << " <.tpgr input file> <n_landmarks>" << std::endl
                << std::endl;

        return EXIT_FAILURE;
    }
    // building landmark

    const std::string tpgr_input_file_name(argv[1]);
    const std::string n_landmarks (argv[2]);

    using Graph = katch::DynamicSearchGraph;
    using EdgeInfo = katch::EdgeInfo<Graph::TTF>;
    std::vector<EdgeInfo> edge_list = katch::tpgr_format::read_edges<EdgeInfo>(tpgr_input_file_name);
    Graph graph(std::move(edge_list));
    std::string main_name = tpgr_input_file_name.substr(0,tpgr_input_file_name.find_last_of("."));
    graph.print_num_of_directed_edge();
    S_Graph out_graph(graph,false);
    S_Graph in_graph(graph,false);
    in_graph.sort_to_incoming_graph();

    generate_landmarks(out_graph, in_graph, stoi(n_landmarks), main_name+".landmark_"+n_landmarks);
    return EXIT_SUCCESS;
}
