//
// Created by Bojie Shen on 26/5/21.
//

#include <iostream>
#include <stdlib.h>
#include <string>
#include <io/tpgr_format_noise.h>
#include <iomanip>

#include "io/edge_info.h"
#include "io/tpgr_format.h"
#include "datastr/graph/dynamic_search_graph.h"
#include "datastr/graph/basic.h"
using Graph = katch::DynamicSearchGraph;
using EdgeInfo = katch::EdgeInfo<Graph::TTF>;
int offset = 5000;

void split_graph_into_n_hours_data(int hour, const Graph& graph, const std::string& output){
    double time_bucket = hour*36000;
    for(int current_hour = 0; current_hour < 24; current_hour ++){
        std::vector<katch::Edge> result_graph (graph.get_n_edges());
        for(katch::NodeIterator it = 0 ; it < graph.get_n_nodes(); it++){
            for(katch::EdgeIterator eit = graph.edges_begin(it); eit != graph.edges_end(it); ++eit ){
                katch::Edge edge ;
                edge.src = it;
                edge.tgt = graph.get_other_node(eit);
                edge.sample_points = graph.get_ttf(eit).get_turning_points({current_hour*36000,current_hour*36000+time_bucket});
                for(auto& p: edge.sample_points){
                    p.x = p.x - current_hour*36000;
                    assert(p.x >= 0 && p.x<= time_bucket );
                }
                result_graph[eit] = edge;
            }
        }

        unsigned n_nodes = 0;
        unsigned n_edges = 0;
        unsigned n_points = 0 ;
        for(const auto&edge : result_graph){
            n_nodes = std::max( n_nodes,edge.tgt);
            n_nodes = std::max( n_nodes,edge.src);
            n_points += edge.sample_points.size();
            n_edges ++;
        }
        n_nodes++;
        std::cout<<n_points<< std::endl;
        std::ofstream myFile(output+ "_" + std::to_string(current_hour)+ ".tpgr");
        myFile<<n_nodes<<" "<<n_edges<<" "<<n_points<<" "<<180000.0 + offset<<"\n";
        for(const auto&edge : result_graph){
            myFile<<edge.src<<" "<<edge.tgt<<" "<<edge.sample_points.size();
            for(const auto& point : edge.sample_points){
                assert(katch::eq((int)point.x,point.x));
                myFile<<std::fixed<<std::setprecision(8)<<" "<< (int)point.x<<" "<< point.y;
            }
            myFile<<"\n";
        }
    }
}


void split_graph_into_hours_data(int hour_durtion, int shifting_hour, const Graph& graph, const std::string& output){
    double time_bucket = shifting_hour*36000 + hour_durtion * 36000;
    for(int current_hour = 0; current_hour < 24; current_hour += hour_durtion){
        std::vector<katch::Edge> result_graph (graph.get_n_edges());
        for(katch::NodeIterator it = 0 ; it < graph.get_n_nodes(); it++){
            for(katch::EdgeIterator eit = graph.edges_begin(it); eit != graph.edges_end(it); ++eit ){
                katch::Edge edge ;
                edge.src = it;
                edge.tgt = graph.get_other_node(eit);
                edge.sample_points = graph.get_ttf(eit).get_turning_points({current_hour*36000,current_hour*36000+time_bucket});
                for(auto& p: edge.sample_points){
                    p.x = p.x - current_hour*36000;
                    assert(p.x >= 0 && p.x<= time_bucket );
                }
                result_graph[eit] = edge;
            }
        }

        unsigned n_nodes = 0;
        unsigned n_edges = 0;
        unsigned n_points = 0 ;
        for(const auto&edge : result_graph){
            n_nodes = std::max( n_nodes,edge.tgt);
            n_nodes = std::max( n_nodes,edge.src);
            n_points += edge.sample_points.size();
            n_edges ++;
        }
        n_nodes++;
        std::cout<<n_points<< std::endl;
        std::ofstream myFile(output+ "_" + std::to_string(hour_durtion) + "_" + std::to_string(shifting_hour) + "_" + std::to_string(current_hour)+ ".tpgr");
        myFile<<n_nodes<<" "<<n_edges<<" "<<n_points<<" "<<time_bucket+offset<<"\n";
        for(const auto&edge : result_graph){
            myFile<<edge.src<<" "<<edge.tgt<<" "<<edge.sample_points.size();
            for(const auto& point : edge.sample_points){
                assert(katch::eq((int)point.x,point.x));
                myFile<<std::fixed<<std::setprecision(8)<<" "<< (int)point.x<<" "<< point.y;
            }
            myFile<<"\n";
        }
    }
}

int main(int argc, char** argv)
{

    const char* binary_name = argv[0];

    if ( argc != 2 )
    {
        std::cerr
                << std::endl
                << "USAGE: " << binary_name
                << "<.tch output file>" << std::endl
                << std::endl;

        return EXIT_FAILURE;
    }

    const std::string tpgr_input_file_name(argv[1]);
    std::vector<EdgeInfo> edge_list = katch::tpgr_format::read_edges<EdgeInfo>(tpgr_input_file_name);

    if ( edge_list.empty() )
    {
        KATCH_ERROR("Empty graph.\n");
        return EXIT_FAILURE;
    }

    Graph graph(std::move(edge_list));
    KATCH_STATUS("Graph has " << graph.get_n_nodes() << " nodes, " << graph.get_n_edges() << " edges.\n");

    // split the graph for building STCH
    split_graph_into_n_hours_data(5, graph, tpgr_input_file_name.substr(0,tpgr_input_file_name.find_last_of(".")));

    // split the graph for building MTCH
    split_graph_into_hours_data(1 , 1, graph, tpgr_input_file_name.substr(0,tpgr_input_file_name.find_last_of(".")));
    split_graph_into_hours_data(1 , 2, graph, tpgr_input_file_name.substr(0,tpgr_input_file_name.find_last_of(".")));
    split_graph_into_hours_data(1 , 4, graph, tpgr_input_file_name.substr(0,tpgr_input_file_name.find_last_of(".")));

    // these are default setting used in the paper
    return EXIT_SUCCESS;
}