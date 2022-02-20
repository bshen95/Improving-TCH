//
// Created by Bojie Shen on 11/5/21.
//

#ifndef TIME_DEPENDENT_STATIC_GRAPH_H
#define TIME_DEPENDENT_STATIC_GRAPH_H
#include <algorithm>
#include <iterator>
#include <utility>
#include <vector>
#include <stack>
#include "datastr/base/ttf_wrapper.h"
#include "datastr/base/pwl_ttf.h"
#include "datastr/graph/basic.h"
#include "datastr/graph/dynamic_search_graph.h"
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
namespace geom = boost::geometry;
typedef geom::model::d2::point_xy<double> point_type;
namespace katch
{

    class StaticSearchGraph
    {

    private:

        struct StaticNode
        {
            EdgeIterator _first_edge;
            StaticNode() : _first_edge(INVALID_EDGE_ITERATOR) {}
        };

        static constexpr NodeIterator UNINITIALIZED_NODE_IT = (1 << 28) - 1;

        struct StaticEdge
        {
            NodeIterator _other_node : 28;
            double _free_flow;

            StaticEdge(const StaticEdge&) = delete;
            StaticEdge& operator= (const StaticEdge&) = delete;

            StaticEdge(StaticEdge&& edge)
            {
                _other_node = edge._other_node;
                _free_flow = edge._free_flow;

            }

            StaticEdge& operator= (StaticEdge&& edge)
            {
                if ( this != &edge )
                {
                    _other_node = edge._other_node;
                    _free_flow = edge._free_flow;
                }

                return *this;
            }

            StaticEdge()
                    : _other_node(UNINITIALIZED_NODE_IT),
                      _free_flow (0)

            {}

            ~StaticEdge()
            {
                _free_flow =0;
            }
        };

        std::vector<StaticNode> _nodes;
        std::vector<StaticEdge> _edges;
        std::pair<double,double>_time_period;

    public:
        std::vector< point_type> _coordinate;
        // convert dynamic graph to static graph, use for constructing CPD.
        StaticSearchGraph(DynamicSearchGraph& d_graph, std::pair<double,double> time_period)
                : _nodes(), _edges()
        {

            _time_period = time_period;
            _nodes.resize(d_graph.get_n_nodes()+1);
            _edges.resize(d_graph.get_n_edges());
            NodeIterator node_index = 0;
            EdgeIterator edge_index = d_graph.edges_begin(node_index);
            std::for_each( _nodes.begin(), _nodes.end()-1,
                           [&node_index,&d_graph](StaticNode& n) -> void
                           {
                               n._first_edge = d_graph.edges_begin(node_index);
                               node_index++;
                           } );
            // don't forget to set the lastest node;
            _nodes[node_index]._first_edge = d_graph.edges_end(node_index-1);
            std::for_each( _edges.begin(), _edges.end(),
                           [&edge_index ,&d_graph, &time_period](StaticEdge& e) -> void
                           {
                                e._other_node = d_graph.get_other_node(edge_index );
                               e._free_flow = d_graph.get_free_flow(edge_index,time_period);
                               edge_index++;
                           } );
            // not sure why we have to insert a dummy edge, but for the sake of consistency, just do it.
            _edges.emplace_back();
            assert(get_n_edges() == d_graph.get_n_edges());
            assert(get_n_nodes() == d_graph.get_n_nodes());
            check_free_flow(d_graph);
        }


        StaticSearchGraph(DynamicSearchGraph& d_graph, bool max_graph)
                : _nodes(), _edges()
        {
            _time_period = std::make_pair(0,864000.0);
            _nodes.resize(d_graph.get_n_nodes()+1);
            _edges.resize(d_graph.get_n_edges());
            NodeIterator node_index = 0;
            EdgeIterator edge_index = d_graph.edges_begin(node_index);
            std::for_each( _nodes.begin(),  _nodes.end()-1,
                           [&node_index,&d_graph](StaticNode& n) -> void
                           {
                               n._first_edge = d_graph.edges_begin(node_index);
                               node_index++;
                           } );
            _nodes[node_index]._first_edge = d_graph.edges_end(node_index-1);
            if(max_graph){
                std::for_each( _edges.begin(), _edges.end(),
                               [&edge_index ,&d_graph](StaticEdge& e) -> void
                               {
                                   e._other_node = d_graph.get_other_node(edge_index );
                                   e._free_flow = d_graph.get_ttf(edge_index).get_max();
                                   edge_index++;
                               } );
            }else{
                std::for_each( _edges.begin(), _edges.end(),
                               [&edge_index ,&d_graph](StaticEdge& e) -> void
                               {
                                   e._other_node = d_graph.get_other_node(edge_index );
                                   e._free_flow = d_graph.get_ttf(edge_index).get_min();
                                   edge_index++;
                               } );
            }
            _edges.emplace_back();
            assert(get_n_edges() == d_graph.get_n_edges());
            assert(get_n_nodes() == d_graph.get_n_nodes());
//            if (!max_graph) check_free_flow(d_graph);
        }

        bool check_free_flow(DynamicSearchGraph& d_graph){
            for ( unsigned i  = 0 ; i <  get_n_edges(); i++){
                for(double j  = _time_period.first; j <= _time_period.second; j += 50){
                    if(lt(d_graph.get_ttf((EdgeIterator)i).eval(j), _edges[i]._free_flow)){
                        std::cout<<"error"<<std::endl;
                        return false;
                    }
                }
            }
            std::cout<<"convert to free flow successfully "<<std::endl;
            return true;
        }

        size_t get_n_nodes() const noexcept { return _nodes.size() - 1; }
        size_t get_n_edges() const noexcept { return _edges.size() - 1; }

        NodeIterator nodes_begin() const noexcept { return 0; }
        NodeIterator nodes_end() const noexcept { return _nodes.size() - 1; }

        EdgeIterator edges_begin(const NodeIterator u) const noexcept
        {
            assert( u + 1 < _nodes.size() );
            return _nodes[u]._first_edge;
        }

        EdgeIterator edges_end(const NodeIterator u) const noexcept
        {
            assert( u + 1 < _nodes.size() );
            return _nodes[u+1]._first_edge;
        }

        NodeIterator get_other_node(const EdgeIterator e) const noexcept
        {
            assert( e + 1 < _edges.size() );
            return _edges[e]._other_node;
        }

        double get_free_flow(const EdgeIterator e) const noexcept
        {
            assert( e + 1 < _edges.size() );
            return _edges[e]._free_flow;
        }

        std::vector<double> get_free_flow_vector(){
            std::vector<double> free_flow = std::vector<double>(get_n_edges());
            EdgeIterator edge_index = (EdgeIterator) 0 ;
            std::for_each( _edges.begin(), _edges.end()-1,
                           [&edge_index ,&free_flow](StaticEdge& e) -> void
                           {
                                free_flow[edge_index] = e._free_flow;
                               edge_index++;
                           } );
            assert(free_flow.size() == get_n_edges());
            return free_flow;
        }

        std::vector<NodeIterator> generate_DFS_ordering() {
            std::vector<NodeIterator> DFS_ordering(get_n_nodes(),0);
            // Mark all the vertices as not visited

            std::vector<bool> was_pushed(get_n_nodes(), false);
            std::vector<bool> was_popped(get_n_nodes(), false);
            std::vector<EdgeIterator>next_out = std::vector<EdgeIterator>(get_n_nodes()+1);
            for(int i = 0; i < _nodes.size(); i++){
                next_out[i] = _nodes[i]._first_edge;
            }
            std::stack<NodeIterator>to_visit;
            NodeIterator next_id = 0;
            for(NodeIterator source_node=0; source_node < get_n_nodes(); ++source_node){
                if(!was_pushed[source_node]){

                    to_visit.push(source_node);
                    was_pushed[source_node] = true;

                    while(!to_visit.empty()){
                        NodeIterator x = to_visit.top();
                        to_visit.pop();
                        if(!was_popped[x]){
                            DFS_ordering[x]= next_id++;
                            was_popped[x] = true;
                        }

                        while(next_out[x] != _nodes[x+1]._first_edge){
                            NodeIterator y = _edges[next_out[x]]._other_node;
                            if(was_pushed[y])
                                ++next_out[x];
                            else{
                                was_pushed[y] = true;
                                to_visit.push(x);
                                to_visit.push(y);
                                ++next_out[x];
                                break;
                            }
                        }
                    }
                }
            }
            DFS_ordering = invert_permutation( DFS_ordering);
            return DFS_ordering;
        }


        void sort_to_incoming_graph(){

            std::vector< std::vector<NodeIterator > > tmp_other_node =    std::vector< std::vector<NodeIterator > >(_nodes.size());
            std::vector< std::vector<double > > tmp_free_flow =    std::vector< std::vector<double > >(_nodes.size());
            for(unsigned i = 0; i < get_n_nodes(); i ++){
                for(EdgeIterator arc = edges_begin(i); arc <  edges_end(i); ++arc ){
                    tmp_other_node[_edges[arc]._other_node].push_back(i);
                    tmp_free_flow[_edges[arc]._other_node].push_back(_edges[arc]._free_flow);
                }
            }
            std::vector<StaticNode> tmp_nodes = std::vector<StaticNode>(_nodes.size());
            std::vector<StaticEdge> tmp_edges = std::vector<StaticEdge>(_edges.size());
            unsigned edge_index = 0;
            for(unsigned i = 0; i < get_n_nodes(); i ++){
                tmp_nodes[i]._first_edge=edge_index;
                for( unsigned j = 0 ; j < tmp_other_node[i].size(); j++){
                    tmp_edges[edge_index]._other_node = tmp_other_node[i][j];
                    tmp_edges[edge_index]._free_flow = tmp_free_flow[i][j];
                    edge_index ++;
                }
            }
            tmp_nodes[_nodes.size()-1]._first_edge = edge_index;
            _nodes = tmp_nodes;
            for(unsigned i  =0; i < _edges.size(); i++){
                _edges[i]._free_flow = tmp_edges[i]._free_flow;
                _edges[i]._other_node = tmp_edges[i]._other_node;
            }

        }

        void resort_graph(const std::vector<NodeIterator> & ordering) {
            //resort the graph based on certain ordering
            //This function should be used for building CPD only.

            //inverted mapper, map current v_id to ordering.
            //v_id -> cpd_id
            std::vector<NodeIterator > inv_ordering = invert_permutation(ordering);

            std::vector<StaticNode> tmp_nodes = std::vector<StaticNode>(_nodes.size());
            std::vector<StaticEdge> tmp_edges = std::vector<StaticEdge>(_edges.size());

//            assert(ordering.size() == get_n_nodes());
            EdgeIterator start_index = edges_begin(0);
            for(unsigned i = 0; i < get_n_nodes(); i ++){
                tmp_nodes[i]._first_edge =start_index;
                NodeIterator v_id = ordering[i];
                for(EdgeIterator arc = edges_begin(v_id); arc <  edges_end(v_id); ++arc ){
                    tmp_edges[start_index]._other_node = inv_ordering[_edges[arc]._other_node];
                    tmp_edges[start_index]._free_flow = _edges[arc]._free_flow;
                    start_index++;
                }
            }
            tmp_nodes[_nodes.size()-1]._first_edge = start_index;
            _nodes = tmp_nodes;
            for(unsigned i  =0; i < _edges.size(); i++){
                _edges[i]._free_flow = tmp_edges[i]._free_flow;
                _edges[i]._other_node = tmp_edges[i]._other_node;
            }


        }

        void  load_coordinate(const std::string& coordinate_file){
            std::ifstream inNodeFile(coordinate_file);
            if(!inNodeFile)
            {
                std::cout << "Cannot open coordinate file" << coordinate_file<< std::endl;
                return ;
            }
            std::cout << "Reading coordinate file: " << coordinate_file<< std::endl;

            int node_num;
            double min_lat, max_lat, min_lon, max_lon;
            inNodeFile >> node_num >> min_lat >> max_lat >> min_lon >> max_lon;
            _coordinate.resize(node_num);
            double id;
            double x,y;
            for(int i = 0; i < node_num; i++) {
                inNodeFile >> id >> x >> y;
                _coordinate[i] =  point_type(y,x);
            }
        }
    };



}
#endif //TIME_DEPENDENT_STATIC_GRAPH_H
