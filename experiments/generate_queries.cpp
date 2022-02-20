//
// Created by Bojie Shen on 19/5/21.
//

#include <iostream>
#include <stdlib.h>
#include <string>
#include <io/edge_info.h>
#include <io/tpgr_format.h>
#include "datastr/graph/dynamic_search_graph.h"
#include <iomanip>

int main(int argc, char **argv) {
    const char* binary_name = argv[0];

    if ( argc != 3 )
    {
        std::cerr
                << std::endl
                << "USAGE: " << binary_name
                << " <.tpgr input file> <n_queries>" << std::endl
                << std::endl;

        return EXIT_FAILURE;
    }

    const std::string tpgr_input_file_name(argv[1]);
    const std::string n_queries (argv[2]);

    using Graph = katch::DynamicSearchGraph;
    using EdgeInfo = katch::EdgeInfo<Graph::TTF>;
    std::vector<EdgeInfo> edge_list = katch::tpgr_format::read_edges<EdgeInfo>(tpgr_input_file_name);
    Graph graph(std::move(edge_list));
    std::string main_name = tpgr_input_file_name.substr(0,tpgr_input_file_name.find_last_of("."));
    // need to provide the coordinate of the vertex.
    graph.load_coordinate(main_name + ".coordinate");
    graph.generate_random_queries(stoi(n_queries),main_name);


    return EXIT_SUCCESS;
}
