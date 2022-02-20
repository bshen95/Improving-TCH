//
// Created by Bojie Shen on 24/5/21.
//

#ifndef KATCH_STATIC_FORWARD_SEARCH_GRAPH_H
#define KATCH_STATIC_FORWARD_SEARCH_GRAPH_H
#include <algorithm>
#include <iterator>
#include <utility>
#include <vector>
#include <stack>
#include "datastr/base/ttf_wrapper.h"
#include "datastr/base/pwl_ttf.h"
#include "datastr/graph/basic.h"
#include "datastr/graph/forward_search_graph.h"
//#include "datastr/graph/hourly_forward_search_graph.h"
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <iostream>
#include <io/vec_io.h>

namespace geom = boost::geometry;
typedef geom::model::d2::point_xy<double> point_type;
namespace katch
{

    class StaticForwardSearchGraph
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
            bool _upward : 1;
            bool _downward : 1;
            double _free_flow;
            unsigned char _first_move;

            StaticEdge(const StaticEdge&) = delete;
            StaticEdge& operator= (const StaticEdge&) = delete;

            StaticEdge(NodeIterator other, bool upward, bool downward, double free_flow, unsigned char _first_move)
            {
                _other_node = other;
                _upward = upward;
                _downward = downward;
                _free_flow = free_flow;
            }

            StaticEdge(StaticEdge&& edge)
            {
                _other_node = edge._other_node;
                _upward = edge._upward;
                _downward = edge._downward;
                _free_flow = edge._free_flow;
                _first_move = edge._first_move;
            }

            StaticEdge& operator= (StaticEdge&& edge)
            {
                if ( this != &edge )
                {
                    _other_node = edge._other_node;
                    _upward = edge._upward;
                    _downward = edge._downward;
                    _free_flow = edge._free_flow;
                    _first_move = edge._first_move;
                }

                return *this;
            }

            StaticEdge()
                    : _other_node(UNINITIALIZED_NODE_IT),
                      _upward(false),
                      _downward(false),
                      _free_flow (0),
                      _first_move(0)

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
        // convert ch search graph to static graph, use for constructing CPD.
        template<const double* period>
        StaticForwardSearchGraph(ForwardSearchGraph<period>& f_graph, std::pair<double,double> time_period)
                : _nodes(), _edges()
        {

            _time_period = time_period;
            _nodes.resize(f_graph.get_n_nodes()+1);
            _edges.resize(f_graph.get_n_edges());
            NodeIterator node_index = 0;
            EdgeIterator edge_index = f_graph.edges_begin(node_index);
            std::for_each( _nodes.begin(), _nodes.end()-1,
                           [&node_index,&f_graph](StaticNode& n) -> void
                           {
                               n._first_edge = f_graph.edges_begin(node_index);
                               node_index++;
                           } );
            // don't forget to set the lastest node;
            _nodes[node_index]._first_edge = f_graph.edges_end(node_index-1);
//            std::vector<double> free_flow_v = f_graph.get_free_flow_vector(time_period);
            std::for_each( _edges.begin(), _edges.end(),
                           [&edge_index, &f_graph, &time_period](StaticEdge& e) -> void
                           {
                               e._other_node = f_graph.get_other_node(edge_index );
                               e._free_flow = f_graph.get_free_flow(edge_index,time_period);
                               e._downward = f_graph.is_directed_downward(edge_index);
                               e._upward = f_graph.is_directed_upward(edge_index);
                               edge_index++;
                           } );
            // not sure why we have to insert a dummy edge, but for the sake of consistency, just do it.
            _edges.emplace_back();
            assert(get_n_edges() == f_graph.get_n_edges());
            assert(get_n_nodes() == f_graph.get_n_nodes());
            check_free_flow(f_graph);
        }

        template<const double* period>
        StaticForwardSearchGraph(ForwardSearchGraph<period>& f_graph, typename ForwardSearchGraph<period>::Type graph_type)
                : _nodes(), _edges()
        {
                    //graph_type 0 = max, 1 = min, 2 = median
            _time_period = std::make_pair(0,*period);
            _nodes.resize(f_graph.get_n_nodes()+1);
            _edges.resize(f_graph.get_n_edges());
            NodeIterator node_index = 0;
            EdgeIterator edge_index = f_graph.edges_begin(node_index);
            std::for_each( _nodes.begin(),  _nodes.end()-1,
                           [&node_index,&f_graph](StaticNode& n) -> void
                           {
                               n._first_edge = f_graph.edges_begin(node_index);
                               node_index++;
                           } );
            _nodes[node_index]._first_edge = f_graph.edges_end(node_index-1);
            {
//                std::vector<double> free_flow_v = f_graph.get_free_flow_vector(max_graph);
                std::for_each(_edges.begin(), _edges.end(),
                              [&edge_index, &f_graph, &graph_type](StaticEdge &e) -> void {
                                  e._other_node = f_graph.get_other_node(edge_index);
                                  e._free_flow = f_graph.get_free_flow(edge_index,graph_type);
                                  e._downward = f_graph.is_directed_downward(edge_index);
                                  e._upward = f_graph.is_directed_upward(edge_index);
                                  edge_index++;
                              });
            }
            _edges.emplace_back();
            assert(get_n_edges() == f_graph.get_n_edges());
            assert(get_n_nodes() == f_graph.get_n_nodes());
//            if (!max_graph) check_free_flow(f_graph);
        }



        template<const double* period>
        bool check_free_flow(ForwardSearchGraph<period>& d_graph){
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

        bool is_directed_upward(const EdgeIterator e) const noexcept
        {
            assert( e < _edges.size() );
            return _edges[e]._upward;
        }

        bool is_directed_downward(const EdgeIterator e) const noexcept
        {
            assert( e + 1 < _edges.size() );
            return _edges[e]._downward;
        }

        unsigned char get_first_move(NodeIterator source_node,  NodeIterator target_node){
            unsigned char fm = 0;
            for ( EdgeIterator e = edges_begin( source_node) ; e != edges_end( source_node) ; ++e )
            {
                if(get_other_node(e) == target_node){
                    return fm;
                }
                fm ++;
            }
            std::cout<<"corresponding edge not found"<<std::endl;
            return 0;
        }
        unsigned char get_first_move(EdgeIterator arc_id){
            bool a  = 0 ;
            return _edges[arc_id]._first_move;
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
                    tmp_edges[start_index]._upward = _edges[arc]._upward;
                    tmp_edges[start_index]._downward = _edges[arc]._downward;
                    start_index++;
                }
            }
            tmp_nodes[_nodes.size()-1]._first_edge = start_index;
            _nodes = tmp_nodes;
            for(unsigned i  =0; i < _edges.size(); i++){
                _edges[i]._free_flow = tmp_edges[i]._free_flow;
                _edges[i]._other_node = tmp_edges[i]._other_node;
                _edges[i]._upward = tmp_edges[i]._upward;
                _edges[i]._downward = tmp_edges[i]._downward;
            }


        }



        void convert_to_income_graph() {
            //convert the outgoing graph to incoming graph.
            struct tmp_arc{
                NodeIterator _other_node : 28;
                bool _upward : 1;
                bool _downward : 1;
                double _free_flow;
                unsigned char _first_move;
            };
            //tmp_graph
            std::vector<std::vector<tmp_arc>> tmp_graph;
            tmp_graph.resize(get_n_nodes());
            for(unsigned i = 0; i < get_n_nodes(); i ++){
                unsigned char first_move =0;
                for(EdgeIterator arc = edges_begin(i); arc <  edges_end(i); ++arc ){
                    unsigned other_node = _edges[arc]._other_node;
                    double  free_flow = _edges[arc]._free_flow;
                    bool upward = !_edges[arc]._upward;
                    bool downward = !_edges[arc]._downward;
                    tmp_arc e = tmp_arc {i,upward,downward,free_flow,first_move};
                    tmp_graph[other_node].push_back(e);
                    first_move++;
                }
            }
            std::vector<StaticNode> tmp_nodes = std::vector<StaticNode>(_nodes.size());
            std::vector<StaticEdge> tmp_edges = std::vector<StaticEdge>(_edges.size());
            EdgeIterator start_index = edges_begin(0);
            for(unsigned i = 0; i < get_n_nodes(); i ++){
                tmp_nodes[i]._first_edge =start_index;
                for(const auto& edge : tmp_graph[i]){
                    tmp_edges[start_index]._other_node = edge._other_node;
                    tmp_edges[start_index]._free_flow = edge._free_flow;
                    tmp_edges[start_index]._upward = edge._upward;
                    tmp_edges[start_index]._downward = edge._downward;
                    tmp_edges[start_index]._first_move = edge._first_move;
                    start_index++;
                }
            }
             tmp_nodes[_nodes.size()-1]._first_edge = start_index;
            _nodes = tmp_nodes;
            for(unsigned i  =0; i < _edges.size(); i++){
                _edges[i]._free_flow = tmp_edges[i]._free_flow;
                _edges[i]._other_node = tmp_edges[i]._other_node;
                _edges[i]._upward = tmp_edges[i]._upward;
                _edges[i]._downward = tmp_edges[i]._downward;
                _edges[i]._first_move = tmp_edges[i]._first_move;
            }

//
//
//            std::vector<StaticNode> tmp_nodes = std::vector<StaticNode>(_nodes.size());
//            std::vector<StaticEdge> tmp_edges = std::vector<StaticEdge>(_edges.size());
//
////            assert(ordering.size() == get_n_nodes());
//            EdgeIterator start_index = edges_begin(0);
//            for(unsigned i = 0; i < get_n_nodes(); i ++){
//                tmp_nodes[i]._first_edge =start_index;
//                NodeIterator v_id = ordering[i];
//                for(EdgeIterator arc = edges_begin(v_id); arc <  edges_end(v_id); ++arc ){
//                    tmp_edges[start_index]._other_node = inv_ordering[_edges[arc]._other_node];
//                    tmp_edges[start_index]._free_flow = _edges[arc]._free_flow;
//                    tmp_edges[start_index]._upward = _edges[arc]._upward;
//                    tmp_edges[start_index]._downward = _edges[arc]._downward;
//                    start_index++;
//                }
//            }
//            tmp_nodes[_nodes.size()-1]._first_edge = start_index;
//            _nodes = tmp_nodes;
//            for(unsigned i  =0; i < _edges.size(); i++){
//                _edges[i]._free_flow = tmp_edges[i]._free_flow;
//                _edges[i]._other_node = tmp_edges[i]._other_node;
//                _edges[i]._upward = tmp_edges[i]._upward;
//                _edges[i]._downward = tmp_edges[i]._downward;
//            }


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

        void convert_graph_to_vector(std::vector<NodeIterator>& first_out,std::vector<NodeIterator>& tail,
                                     std::vector<double>& free_flow, std::vector<bool>& is_downward_arc){
            for(auto node : _nodes){
                first_out.push_back(node._first_edge);
            }
            for(const auto& edge : _edges){
                tail.push_back(edge._other_node);
                free_flow.push_back(edge._free_flow);
                is_downward_arc.push_back(edge._downward);
            }
            first_out.shrink_to_fit();
            tail.shrink_to_fit();
            free_flow.shrink_to_fit();
            is_downward_arc.shrink_to_fit();
        }
    };



}
#endif //KATCH_STATIC_FORWARD_SEARCH_GRAPH_H
