//
// Created by Bojie Shen on 23/5/21.
//

#ifndef KATCH_FORWARD_SEARCH_GRAPH_H
#define KATCH_FORWARD_SEARCH_GRAPH_H
#include <algorithm>
#include <iterator>
#include <utility>
#include <vector>
#include <iostream>
#include <stack>
#include "datastr/base/ttf_wrapper.h"
#include "datastr/base/pwl_ttf.h"
#include "datastr/graph/basic.h"
#include "io/vec_io.h"
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
namespace geom = boost::geometry;
typedef geom::model::d2::point_xy<double> point_type;

namespace katch
{

    template<const double* period>
    class ForwardSearchGraph
    {
    private:

        static constexpr double PERIOD = *period;
        static constexpr int BUCKET_SIZE = 8;

        using TTFImpl = PwlTTF<static_cast<int>(PERIOD), BUCKET_SIZE>;

    public:

        using TTFRef = TTFReference<TTFImpl>;
        using TTF = TTFWrapper<TTFImpl>;

        enum Type
        {
            max_graph,
            min_graph,
            average_graph
        };
    private:

        struct Node
        {
            EdgeIterator _first_edge;
            Node() : _first_edge(INVALID_EDGE_ITERATOR) {}
        };

        static constexpr NodeIterator UNINITIALIZED_NODE_IT = (1 << 28) - 1;

        struct Edge
        {
            NodeIterator _other_node : 28;
            bool _upward : 1;
            bool _downward : 1;
            bool _is_exact : 1;
            bool _is_constant : 1;
            uint32_t _first_middle_node;

            union
            {
                double _constant_value;
                const TTFImpl* _ttf_impl_ptr;
                struct { uint32_t _lower; uint32_t _upper; };
            };

            Edge(const Edge&) = delete;
            Edge& operator= (const Edge&) = delete;

            double get_free_flow(std::pair<double,double> time_period){
                if(_is_constant)
                    return _constant_value;
                else
                    return _ttf_impl_ptr->get_min(time_period);
            }

            Edge(Edge&& edge)
            {
                _other_node = edge._other_node;
                _upward = edge._upward;
                _downward = edge._downward;
                _is_exact = edge._is_exact;
                _is_constant = edge._is_constant;
                _first_middle_node = edge._first_middle_node;

                if ( _is_constant && _is_exact )
                {
                    _constant_value = edge._constant_value;
                }
                else if ( _is_constant && ! _is_exact )
                {
                    _lower = edge._lower;
                    _upper = edge._upper;
                }
                else
                {
                    assert( ! _is_constant );

                    _ttf_impl_ptr = edge._ttf_impl_ptr;
                    edge._ttf_impl_ptr = nullptr;
                }
            }

            Edge& operator= (Edge&& edge)
            {
                if ( this != &edge )
                {
                    if ( ! _is_constant )
                        if ( _ttf_impl_ptr )
                            delete _ttf_impl_ptr;

                    _other_node = edge._other_node;
                    _upward = edge._upward;
                    _downward = edge._downward;
                    _is_exact = edge._is_exact;
                    _is_constant = edge._is_constant;
                    _first_middle_node = edge._first_middle_node;

                    if ( _is_constant && _is_exact )
                    {
                        _constant_value = edge._constant_value;
                    }
                    else if ( _is_constant && ! _is_exact )
                    {
                        _lower = edge._lower;
                        _upper = edge._upper;
                    }
                    else
                    {
                        assert( ! _is_constant );

                        _ttf_impl_ptr = edge._ttf_impl_ptr;
                        edge._ttf_impl_ptr = nullptr;
                    }
                }

                return *this;
            }

            Edge()
                    : _other_node(UNINITIALIZED_NODE_IT),
                      _upward(false),
                      _downward(false),
                      _is_exact(true),
                      _is_constant(false),
                      _first_middle_node(0),
                      _ttf_impl_ptr(nullptr)
            {}

            ~Edge()
            {
                if ( ! _is_constant )
                    if ( _ttf_impl_ptr )
                        delete _ttf_impl_ptr;

                _is_constant = false;
                _ttf_impl_ptr = nullptr;
            }


        };

        std::vector<Node> _nodes;
        std::vector<Node> _downward_nodes;
        std::vector<Node> _upward_nodes;
        std::vector<bool> _is_lowest_level_nodes;
        std::vector<Edge> _edges;
        std::vector<MiddleNodeDescr> _middle_node_data;

        std::vector<unsigned > _ranking;
        double _lower_bound_factor;

        struct tmp_edge{
            NodeIterator  src;
            NodeIterator  tgt;
            bool _upward : 1;
            bool _downward : 1;
            bool _is_exact : 1;
            bool _is_constant : 1;
            std::vector<MiddleNodeDescr> middle_nodes_data;
            union
            {
                double _constant_value;
                const TTFImpl* _ttf_impl_ptr;
                struct { uint32_t _lower; uint32_t _upper; };
            };
        };
    public:
        double min_lat, max_lat, min_lon, max_lon;
        std::vector< point_type> _coordinate;
        template <typename EdgeInfo>
        ForwardSearchGraph(std::vector<EdgeInfo>&& edge_list)
                : _nodes(), _edges(), _middle_node_data(),
                  _lower_bound_factor( 1.0)
        {

            std::for_each( edge_list.begin(), edge_list.end(),
                           [](EdgeInfo& edge_info) -> void
                           {
                                if(edge_info.get_backward()){
                                    //backward arc is incoming arc, so for the forward search flip it.
                                    NodeIterator tmp_src = edge_info.get_source();
                                    edge_info.set_source(edge_info.get_target());
                                    edge_info.set_target(tmp_src);
                                }

                           } );
            size_t max_node_it = 0;
            std::for_each( edge_list.begin(), edge_list.end(),
                           [&max_node_it](const EdgeInfo& edge_info) -> void
                           {
                               max_node_it = std::max(max_node_it, edge_info.get_source());
                               max_node_it = std::max(max_node_it, edge_info.get_target());
                           } );

            _nodes.assign(max_node_it + 2, Node());
            _downward_nodes.assign(max_node_it + 2, Node());
            _upward_nodes.assign(max_node_it + 2, Node());
            std::sort( edge_list.begin(), edge_list.end(),
                       [](const EdgeInfo& lhs, const EdgeInfo& rhs) -> bool
                       {
                           if(lhs.get_source() == rhs.get_source()){
                               if(lhs.get_backward() && rhs.get_forward()) {
                                   return true;
                               }else if (lhs.get_backward() && rhs.get_backward()) {
                                   return lhs.get_target() < rhs.get_target();
                               }else if (lhs.get_forward() && rhs.get_forward()) {
                                   return lhs.get_target() < rhs.get_target();
                               }else{
                                   return false;
                               }
                           }else{
                               return lhs.get_source() < rhs.get_source();
                           }


                       } );
//            unsigned counter = 0;
//            for ( auto it = edge_list.begin() ; it != edge_list.end() ; ++it )
//            {
//                std::cout<<"Source: "<< it->get_source() << "Target "<<it->get_target() << "Upward: "<< it->get_forward()
//                <<std::endl;
//                if(counter > 50) break;
//
//                counter++;
//            }
            assert( _nodes[0]._first_edge == INVALID_EDGE_ITERATOR );
            _nodes[0]._first_edge = 0;
            _downward_nodes[0]._first_edge = 0;

            NodeIterator current_source = 0;
            NodeIterator current_downward_source = 0;
            NodeIterator current_upward_source = 0;
            for ( auto it = edge_list.begin() ; it != edge_list.end() ; ++it )
            {
                assert( it->get_source() + 1 < _nodes.size() );

                while ( current_source <= it->get_source() )
                    _nodes[++current_source]._first_edge = _edges.size();

                while ( current_downward_source <= it->get_source() )
                    _downward_nodes[++current_downward_source]._first_edge = _edges.size();

                while( current_upward_source < it->get_source()){
                    _upward_nodes[current_upward_source]._first_edge = _edges.size();
                    current_upward_source++;
                }
                _edges.push_back(std::move(Edge()));
                _edges.back()._other_node = it->get_target();
                _edges.back()._upward = it->get_forward();
                _edges.back()._downward = it->get_backward();
                _edges.back()._is_exact = true;
                if(it->get_forward()){
                        if( current_upward_source == it->get_source() ){
                            _upward_nodes[current_upward_source]._first_edge = _edges.size()-1;
                            current_upward_source++;
                        }
                }
                _edges.back()._first_middle_node = _middle_node_data.size();

                _middle_node_data.insert(
                        _middle_node_data.end(), it->get_middle_node_data().begin(), it->get_middle_node_data().end());

                if ( it->get_ttf().is_constant() )
                {
                    _edges.back()._is_constant = true;
                    _edges.back()._constant_value = it->get_ttf().get_constant_value();

                }
                else
                {
                    assert( it->get_ttf().get_ttf_impl_ptr() );

                    _edges.back()._is_constant = false;
                    _edges.back()._ttf_impl_ptr = it->get_ttf().get_ttf_impl_ptr();
                    it->get_ttf().release();
                }

                ++_nodes[current_source]._first_edge;
                ++_downward_nodes[current_downward_source]._first_edge;
            }

            while ( current_source <= max_node_it )
                _nodes[++current_source]._first_edge = _edges.size();

            while ( current_downward_source <= max_node_it )
                _downward_nodes[++current_downward_source]._first_edge = _edges.size();

            while ( current_upward_source <= max_node_it )
                _upward_nodes[++current_upward_source]._first_edge = _edges.size();


            assert( _nodes.back()._first_edge == _edges.size() );

            _edges.emplace_back();
            _edges.back()._first_middle_node = _middle_node_data.size();
            edge_list.clear();
            check_fw_bw_iterator();
            check_lowest_level_nodes();

            // set flag for lowest level node;
            _is_lowest_level_nodes = std::vector<bool>(_nodes.size()-1,false);
            unsigned num_vertices = get_n_nodes();
            unsigned num = 0;
            for(unsigned i = 0; i < num_vertices; i ++){
                bool have_downward_arc = false;
                for ( unsigned e = edges_begin(i); e != edges_end(i); ++e ){
                    if(is_directed_downward((EdgeIterator) e)){
                        have_downward_arc = true;
                        break;
                    }
                }
                if(!have_downward_arc){
                    _is_lowest_level_nodes[i] = true;
                }
            }
        }


        bool is_lowest_level_node (NodeIterator i ){
            return _is_lowest_level_nodes[i];
        }
        void check_fw_bw_iterator(){
            unsigned num_vertices = get_n_nodes();
            for(unsigned i = 0; i < num_vertices; i ++){
                unsigned num_edges = edges_end(i)-edges_begin(i);
                unsigned bw_start = downward_begin(i);
                unsigned bw_end = downward_end(i);
                for ( unsigned e = bw_start ; e != bw_end ; ++e ){
                    if(_edges[e]._upward){
                        std::cout<<"bw_iterator_error";
                    }
                    if(!_edges[e]._downward){
                        std::cout<<"bw_iterator_error";
                    }
                }
                unsigned fw_start = upward_begin(i);
                unsigned fw_end = upward_end(i);
                for ( unsigned e = fw_start ; e != fw_end ; ++e ){
                    if(!_edges[e]._upward){
                        std::cout<<"fw_iterator_error";
                    }
                    if(_edges[e]._downward){
                        std::cout<<"fw_iterator_error";
                    }
                }
                if(num_edges != bw_end-bw_start + fw_end-fw_start){
                    std::cout<<"edge missing";
                }
            }


        }
        void check_lowest_level_nodes(){
            unsigned num_vertices = get_n_nodes();
            unsigned num = 0;
            for(unsigned i = 0; i < num_vertices; i ++){
                bool have_downward_arc = false;
                for ( unsigned e = edges_begin(i); e != edges_end(i); ++e ){
                   if(is_directed_downward((EdgeIterator) e)){
                       have_downward_arc = true;
                       break;
                   }
                }
                if(!have_downward_arc){
                    num ++;
                }
            }
            std::cout<<"Total number of lowest level nodes : "<< num <<std::endl;
        }

        std::vector<bool> get_lowest_level_mapper(){
            std::vector<bool> _lowest_nodes (_nodes.size()-1,false);
            unsigned num_vertices = get_n_nodes();
            unsigned num = 0;
            for(unsigned i = 0; i < num_vertices; i ++){
                bool have_downward_arc = false;
                for ( unsigned e = edges_begin(i); e != edges_end(i); ++e ){
                    if(is_directed_downward((EdgeIterator) e)){
                        have_downward_arc = true;
                        break;
                    }
                }
                if(!have_downward_arc){
                    _lowest_nodes[i] = true;
                }
            }
            return _lowest_nodes;

        }

        void set_ranking(std::vector<unsigned> ranking){
            _ranking = ranking;
        }

        const unsigned& get_ranking(NodeIterator i){
           return _ranking[i];
        }

        unsigned get_index_size(){
            unsigned size = 0;
            size += _nodes.size()* sizeof(_nodes[0]._first_edge);
            for(auto & _edge : _edges){
//            NodeIterator _other_node : 28;
//            bool _upward : 1;
//            bool _downward : 1;
//            bool _is_exact : 1;
//            bool _is_constant : 1;
//            uint32_t _first_middle_node;
                // 4 bytes  for NodeIterator; 1 byte for 4 bool; 4 bytes for uni32_t
                size +=  4 +  3 + 4;
                if(_edge._is_constant){
                    // each constant value double is 8 bytes
                    size += 8;
                }else{
                    if(_edge._ttf_impl_ptr != nullptr) {
                        // each point is (double, double), 16 bytes
                        size += _edge._ttf_impl_ptr->get_n_points() * 16;
                        // each bucket is int, 4 bytes
                        size += _edge._ttf_impl_ptr->get_n_bucket() * 4;
                        // bucket shift uint32_t
                        size += 4;
                    }
                };
                if(!_edge._is_exact){
                    std::cout<<"wtf! why!?" <<std::endl;
                }
            }
            //    double _time : 8
            //    NodeIterator _middle_node : 4
            size += _middle_node_data.size()* 12;

            return size;
        }
        size_t get_n_nodes() const noexcept { return _nodes.size() - 1; }
        size_t get_n_edges() const noexcept { return _edges.size() - 1; }

        NodeIterator nodes_begin() const noexcept { return 0; }
        NodeIterator nodes_end() const noexcept { return _nodes.size() - 1; }

        EdgeIterator upward_begin(const NodeIterator u) const noexcept
        {
            assert( u + 1 < _upward_nodes.size() );
            return _upward_nodes[u]._first_edge;
        }

        EdgeIterator upward_end(const NodeIterator u) const noexcept
        {
            assert( u + 1 < _downward_nodes.size() );
            return _downward_nodes[u+1]._first_edge;
        }


        EdgeIterator downward_begin(const NodeIterator u) const noexcept
        {
            assert( u + 1 < _downward_nodes.size() );
            return _downward_nodes[u]._first_edge;
        }

        EdgeIterator downward_end(const NodeIterator u) const noexcept
        {
            assert( u + 1 < _upward_nodes.size() );
            return _upward_nodes[u]._first_edge;
        }

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

        double get_free_flow(const EdgeIterator& e, std::pair<double,double> time_period){
            // period free flow on an edge
            assert( e + 1 < _edges.size() );
            return _edges[e].get_free_flow(time_period);
        }

        double get_free_flow(const EdgeIterator& e, Type graph_type){
            // period free flow on an edge
            assert( e + 1 < _edges.size() );
            if(graph_type == Type::max_graph){
                return get_ttf(e).get_max();
            }else if(graph_type == Type::min_graph){
                return get_ttf(e).get_min();
            }else if (graph_type == Type::average_graph){
                // set to 5 mins;
                return get_ttf(e).get_average(5);
            }
        }

        double get_free_flow(const EdgeIterator& e, bool is_max){
            // period free flow on an edge
            assert( e + 1 < _edges.size() );
            return is_max? get_ttf(e).get_max() : get_ttf(e).get_min();
        }

        std::vector<double> get_free_flow_vector(bool _is_max){
            KATCH_STATUS("Compute the shortcut free flow based on original arc");
            // converted to free flow based on original arc;
            std::vector<double> _edge_free_flow (get_n_edges(),0);
            std::vector<NodeIterator> _edge_middle_nodes (get_n_edges(),INVALID_NODE_ITERATOR);
            std::vector<bool> _is_edge_settled (get_n_edges(), false);
            unsigned unsettled_edge = get_n_edges();
            while(unsettled_edge != 0){
//                KATCH_CONTINUE_STATUS("Edges left:" + std::to_string(unsettled_edge));
                for(NodeIterator u = 0; u < get_n_nodes(); u ++){
                    for(EdgeIterator e = edges_begin(u); e != edges_end(u); e++){
                        if(_is_edge_settled[e]){
                            continue;
                        }
                        NodeIterator v = get_other_node(e);
                        if(get_n_middle_nodes(e) == 0){
                            _is_edge_settled[e] = true;
                            _edge_free_flow[e] = _is_max? get_ttf(e).get_max() : get_ttf(e).get_min();
                            unsettled_edge--;
                        }else{
                            bool are_shortcuts_settled = true;
                            double min_free_flow = std::numeric_limits<double>::max();
                            NodeIterator middle_node = INVALID_NODE_ITERATOR;
                            for(unsigned eit = _edges[e]._first_middle_node; eit != _edges[e+1]._first_middle_node; ++eit ){
                                // shortcut edge
                                NodeIterator  m = _middle_node_data[eit]._middle_node;
                                if(m == INVALID_NODE_ITERATOR){
                                    min_free_flow = _is_max? get_ttf(e).get_max() : get_ttf(e).get_min();
                                    continue;
                                }

                                EdgeIterator e_up = INVALID_EDGE_ITERATOR;
                                EdgeIterator e_down = INVALID_EDGE_ITERATOR;
                                for ( EdgeIterator it = edges_begin(m) ; it != edges_end(m) ; ++it )
                                {
                                    if ( get_other_node(it) == v && is_directed_upward(it) ) e_up = it;
                                    if ( e_up != INVALID_EDGE_ITERATOR ) break;
                                }
                                for ( EdgeIterator it = edges_begin(u) ; it != edges_end(u) ; ++it )
                                {
                                    if ( get_other_node(it) == m && is_directed_downward(it) ) e_down = it;
                                    if ( e_down != INVALID_EDGE_ITERATOR) break;
                                }
                                if (_is_edge_settled[e_down] && _is_edge_settled[e_up]){
                                    middle_node = m;
                                    min_free_flow =  std::min(_edge_free_flow[e_up] + _edge_free_flow[e_down], min_free_flow);
                                }else{
                                    are_shortcuts_settled = false;
                                    break;
                                }
                            }
                            if(are_shortcuts_settled){
                                _is_edge_settled[e] = true;
                                _edge_free_flow[e] = min_free_flow;
                                _edge_middle_nodes[e] = middle_node;
                                unsettled_edge--;
                            }
                        }
                    }
                }
            }
            for( auto b : _is_edge_settled){
                if(!b){
                    std::cout<<"error"<<std::endl;
                }
            }
            for( auto b : _edge_free_flow){
                if( b==0){
                    std::cout<<"error"<<std::endl;
                }
            }
            return _edge_free_flow;
        }

        std::vector<double> get_free_flow_vector(std::pair<double,double> time_period){
            KATCH_STATUS("Compute the shortcut free flow based on original arc");
            // converted to free flow based on original arc;
            std::vector<double> _edge_free_flow (get_n_edges(),0);
            std::vector<bool> _is_edge_settled (get_n_edges(), false);
            unsigned unsettled_edge = get_n_edges();
            while(unsettled_edge != 0){
//                KATCH_CONTINUE_STATUS("Edges left:" + std::to_string(unsettled_edge));
                for(NodeIterator u = 0; u < get_n_nodes(); u ++){
                    for(EdgeIterator e = edges_begin(u); e != edges_end(u); e++){
                        if(_is_edge_settled[e]){
                            continue;
                        }
                        NodeIterator v = get_other_node(e);
                        if(get_n_middle_nodes(e) == 0){
                            _is_edge_settled[e] = true;
                            _edge_free_flow[e] = get_free_flow(e,time_period);
                            unsettled_edge--;
                        }else{
                            bool are_shortcuts_settled = true;
                            double min_free_flow = std::numeric_limits<double>::max();
                            for(unsigned eit = _edges[e]._first_middle_node; eit != _edges[e+1]._first_middle_node; ++eit ){
                                // shortcut edge
                                NodeIterator  m = _middle_node_data[eit]._middle_node;
                                if(m == INVALID_NODE_ITERATOR){
                                    min_free_flow =  get_free_flow(e,time_period);
                                    continue;
                                }

                                EdgeIterator e_up = INVALID_EDGE_ITERATOR;
                                EdgeIterator e_down = INVALID_EDGE_ITERATOR;
                                for ( EdgeIterator it = edges_begin(m) ; it != edges_end(m) ; ++it )
                                {
                                    if ( get_other_node(it) == v && is_directed_upward(it) ) e_up = it;
                                    if ( e_up != INVALID_EDGE_ITERATOR ) break;
                                }
                                for ( EdgeIterator it = edges_begin(u) ; it != edges_end(u) ; ++it )
                                {
                                    if ( get_other_node(it) == m && is_directed_downward(it) ) e_down = it;
                                    if ( e_down != INVALID_EDGE_ITERATOR) break;
                                }
                                if (_is_edge_settled[e_down] && _is_edge_settled[e_up]){
                                    min_free_flow =  std::min(_edge_free_flow[e_up] + _edge_free_flow[e_down],min_free_flow);
                                }else{
                                    are_shortcuts_settled = false;
                                    break;
                                }
                            }
                            if(are_shortcuts_settled){
                                _is_edge_settled[e] = true;
                                _edge_free_flow[e] = min_free_flow;
                                unsettled_edge--;
                            }
                        }
                    }
                }
            }
            for( auto b : _is_edge_settled){
                if(!b){
                    std::cout<<"error"<<std::endl;
                }
            }
            for( auto b : _edge_free_flow){
                if( b==0){
                    std::cout<<"error"<<std::endl;
                }
            }
            return _edge_free_flow;
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

        double get_upper(const EdgeIterator e) const noexcept
        {
            assert( e + 1 < _edges.size() );
            if ( _edges[e]._is_exact )
            {
                if ( _edges[e]._is_constant )
                    return _edges[e]._constant_value;

                assert( _edges[e]._ttf_impl_ptr );
                return _edges[e]._ttf_impl_ptr->get_max();
            }

            if ( _edges[e]._is_constant )
                return double(_edges[e]._upper);

            return _edges[e]._ttf_impl_ptr->get_max();
        }

        double get_lower(const EdgeIterator e) const noexcept
        {
            assert( e + 1 < _edges.size() );
            if ( _edges[e]._is_exact )
            {
                if ( _edges[e]._is_constant )
                    return _edges[e]._constant_value;

                assert( _edges[e]._ttf_impl_ptr );
                return _edges[e]._ttf_impl_ptr->get_min();
            }

            if ( _edges[e]._is_constant )
                return double(_edges[e]._lower);

            return _lower_bound_factor * _edges[e]._ttf_impl_ptr->get_min();
        }

        TTFRef get_ttf(const EdgeIterator e) const noexcept
        {
            assert( e + 1 < _edges.size() );
            assert( _edges[e]._is_exact );

            if ( _edges[e]._is_constant )
                return TTFRef(_edges[e]._constant_value);

            assert( _edges[e]._ttf_impl_ptr );
            return TTFRef(_edges[e]._ttf_impl_ptr);
        }

        TTFRef get_upper_bound(const EdgeIterator e) const noexcept
        {
            assert( e + 1 < _edges.size() );
            assert( ! _edges[e]._is_exact );

            if ( _edges[e]._is_constant )
                return TTFRef(double(_edges[e]._upper));

            assert( _edges[e]._ttf_impl_ptr );
            return TTFRef(_edges[e]._ttf_impl_ptr);
        }

        TTF get_lower_bound(const EdgeIterator e) const noexcept
        {
            assert( e + 1 < _edges.size() );
            assert( ! _edges[e]._is_exact );

            if ( _edges[e]._is_constant )
                return TTF(double(_edges[e]._lower));

            assert( _edges[e]._ttf_impl_ptr );
            return TTF(_edges[e]._ttf_impl_ptr->scale(_lower_bound_factor));
        }

        size_t get_n_middle_nodes(const EdgeIterator e) const noexcept
        {
            assert( e + 1 < _edges.size() );
            return _edges[e+1]._first_middle_node - _edges[e]._first_middle_node;
        }

        MiddleNodeDescr get_middle_node_descr(const EdgeIterator e, const size_t index) const noexcept
        {
            assert( e + 1 < _edges.size() );

            if ( get_n_middle_nodes(e) == 0 )
                return MiddleNodeDescr(0.0, INVALID_NODE_ITERATOR);

            assert( index < get_n_middle_nodes(e) );
            const size_t position = _edges[e]._first_middle_node + index;

            assert( position < _edges[e+1]._first_middle_node );
            return _middle_node_data[position];
        }

        NodeIterator get_middle_node(const EdgeIterator e, const double time) const noexcept
        {
            assert( e + 1 < _edges.size() );

            if ( get_n_middle_nodes(e) == 0 )
                return INVALID_NODE_ITERATOR;

            if ( get_n_middle_nodes(e) == 1 )
                return _middle_node_data[_edges[e]._first_middle_node]._middle_node;

            const double x = util::modulo(time, PERIOD);

            assert( _edges[e]._first_middle_node <= _middle_node_data.size() );
            assert( _edges[e+1]._first_middle_node <= _middle_node_data.size() );
            assert( _edges[e]._first_middle_node <= _edges[e+1]._first_middle_node );

            uint32_t index;
            for ( index = _edges[e]._first_middle_node ; index != _edges[e+1]._first_middle_node ; ++index )
                if ( x < _middle_node_data[index]._time )
                    break;

            if ( index == _edges[e]._first_middle_node )
            {
                assert( _edges[e+1]._first_middle_node > 0 );
                return _middle_node_data[_edges[e+1]._first_middle_node - 1]._middle_node;
            }

            if ( index + 1 == _edges[e+1]._first_middle_node )
            {
                assert( _edges[e+1]._first_middle_node > _edges[e]._first_middle_node );

                if ( le(_middle_node_data[_edges[e+1]._first_middle_node - 1]._time, x) )
                    return _middle_node_data[_edges[e+1]._first_middle_node - 1]._middle_node;

                assert( x < _middle_node_data[_edges[e+1]._first_middle_node - 1]._time );
            }

            assert( index > 0 );
            assert( index - 1 >= _edges[e]._first_middle_node );
            assert( index - 1 < _edges[e+1]._first_middle_node );

            assert( le(_middle_node_data[index - 1]._time, x) );
            assert( KATCH_IMPLIES(index < _edges[e+1]._first_middle_node, le( x, _middle_node_data[index]._time)) );

            return _middle_node_data[index - 1]._middle_node;
        }
        void print_statistic() {
            unsigned n_constant = 0;
            unsigned n_td = 0;
            unsigned n_points = 0;
            for(const auto& edge: _edges){
                if(edge._is_constant){
                    n_constant++;
                }else{
                    if(edge._ttf_impl_ptr != nullptr){
                        n_points += edge._ttf_impl_ptr->get_n_points();
                        n_td ++;
                    }
                }
            }
            std::cout<<"Total number of nodes: "<<this->get_n_nodes() << std::endl;
            std::cout<<"Total number of edges: "<<this->get_n_edges() << std::endl;
            std::cout<<"Total number of constant edges: "<<n_constant <<" Percentage: "<< (double)n_constant/this->get_n_edges() << std::endl;
            std::cout<<"Total number of td edges: "<<n_td <<" Percentage: "<< (double)n_td/this->get_n_edges()<< std::endl;
            std::cout<<"Total number of turning points : "<<n_points<<std::endl;
            std::cout<<"Average turning points per edges: "<< (double)n_points/this->get_n_edges() << std::endl;
            std::cout<<"Average turning points per td edges: "<< (double)n_points/n_td << std::endl;
        }

        void print_edge_traffic_index(const std::string& filename,unsigned id) {
            // show traffic index per 20 mins;
            std::vector<double> traffic(24*12);
            const auto& edge = _edges[(EdgeIterator) id];
            if(!edge._is_constant){
                if(edge._ttf_impl_ptr != nullptr){
                    for(int i = 0; i <24*12; i++){
                        traffic[i] += (double)edge._ttf_impl_ptr->eval(i*600*5);
                    }
                }
            }else{
                for(int i = 0; i <24*12; i++){
                    traffic[i] += (double)edge._constant_value;
                }

            }

            std::ofstream myFile(filename);
            myFile<<"time,traffic_index\n";
            for(int i = 0; i < traffic.size(); i++){
                myFile<<std::fixed<<std::setprecision(8)<<i<<","<<(double)traffic[i]
                      <<"\n";
            }
            myFile.close();
        }

        void print_edge_chaining_traffic_index(const std::string& filename,unsigned id,unsigned id2) {
            // show traffic index per 20 mins;
            std::vector<double> traffic(24*12);
            auto ttf = get_ttf((EdgeIterator)id);
            auto ttf2 = get_ttf((EdgeIterator)id2);
            TTF f_v_new = link(ttf,ttf2);

            for(int i = 0; i <24*12; i++){
                traffic[i] += (double)f_v_new.eval(i*600*5);
            }

            std::ofstream myFile(filename);
            myFile<<"time,traffic_index\n";
            for(int i = 0; i < traffic.size(); i++){
                myFile<<std::fixed<<std::setprecision(8)<<i<<","<<(double)traffic[i]
                      <<"\n";
            }
            myFile.close();
        }

        void check_time_cost(){
            int j = 0;
            for(auto& edge: _edges){
                for(double i = 0 ; i < 864000; i += 12000){
                    if(!edge._is_constant){
                        if(edge._ttf_impl_ptr != nullptr){
                            if(!eq(edge._ttf_impl_ptr->eval(i), edge._ttf_impl_ptr->eval(i+864000))){
                                std::cout<<"Error"<<std::endl;
                            }
                        }
                    }


                }
                j++;
            }

            std::cout<<"Passed"<<std::endl;
        }

        bool is_down_reachable(const NodeIterator s_node, const NodeIterator t_node){

            std::stack<NodeIterator>to_visit;
            std::vector<bool> was_popped(get_n_nodes(), false);

            to_visit.push(s_node);
            while(!to_visit.empty()){
                NodeIterator x = to_visit.top();
                to_visit.pop();
                if(x == t_node){
                    return true;
                }

                if(was_popped[x]){
                    continue;
                }else{
                    was_popped[x] = true;
                }

                for ( EdgeIterator e = edges_begin(x) ; e != edges_end(x) ; ++e ) {
                    if(is_directed_downward(e)){
                        NodeIterator other_node = get_other_node(e);
                        if(!was_popped[other_node]){
                            was_popped[other_node] = true;
                            to_visit.push(other_node);
                        }
                    }
                }

            }
            return  false;
        }
        std::vector<NodeIterator> generate_ch_based_dfs_ordering(std::vector<unsigned> ch_ordering ) {
            // ch - > node_id
            std::vector<unsigned> ch_order_to_node = invert_permutation(ch_ordering);

            // sort next other node;

            std::vector<NodeIterator>next_node= std::vector<NodeIterator>(get_n_edges()+1);
            for(unsigned i = 0; i < next_node.size();i++ ){
                next_node[i] =  _edges[i]._other_node;
            }
            // node with larger ch_ordering first;
            for(unsigned i = 0; i < get_n_nodes();i ++){
                unsigned edge_begin = edges_begin(i);
                unsigned edge_end = edges_end(i);
                std::sort(next_node.begin()+edge_begin, next_node.begin()+ edge_end, [& ch_ordering](const unsigned & id_1, const unsigned & id_2) -> bool
                {
                    return ch_ordering[id_1] > ch_ordering[id_2] ;
                });
                for(unsigned arc = edge_begin; arc  < edge_end-1; arc++){
                    if(ch_ordering[next_node[arc]] < ch_ordering[next_node[arc+1]]){
                        std::cout<<"error"<<std::endl;
                    }
                }
            }




            std::vector<NodeIterator> DFS_ordering(get_n_nodes(),0);
            // Mark all the vertices as not visited
            std::vector<EdgeIterator>next_out = std::vector<EdgeIterator>(get_n_nodes()+1);
            for(int i = 0; i < _nodes.size(); i++){
                next_out[i] = _nodes[i]._first_edge;
            }
            std::vector<bool> was_pushed(get_n_nodes(), false);
            std::vector<bool> was_popped(get_n_nodes(), false);
            std::stack<NodeIterator>to_visit;
            NodeIterator next_id = 0;
            for(int s = get_n_nodes()-1; s >= 0; --s){
                NodeIterator source_node = ch_order_to_node[s];
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
                            assert(next_out[x]+1 < _edges.size());
                            NodeIterator y = next_node[next_out[x]];
                            // downward edge only;
                            if(ch_ordering[y] >= ch_ordering[x]){
                                ++next_out[x];
                                continue;
                            }
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


        std::vector<NodeIterator> generate_upward_ch_based_dfs_ordering(std::vector<unsigned> ch_ordering ) {
            // ch - > node_id
            std::vector<unsigned> ch_order_to_node = invert_permutation(ch_ordering);

            // sort next other node;

            std::vector<NodeIterator>next_node= std::vector<NodeIterator>(get_n_edges()+1);
            for(unsigned i = 0; i < next_node.size();i++ ){
                next_node[i] =  _edges[i]._other_node;
            }
            // node with larger ch_ordering first;
            for(unsigned i = 0; i < get_n_nodes();i ++){
                unsigned edge_begin = edges_begin(i);
                unsigned edge_end = edges_end(i);
                std::sort(next_node.begin()+edge_begin, next_node.begin()+ edge_end, [& ch_ordering](const unsigned & id_1, const unsigned & id_2) -> bool
                {
                    return ch_ordering[id_1] < ch_ordering[id_2] ;
                });
            }




            std::vector<NodeIterator> DFS_ordering(get_n_nodes(),0);
            // Mark all the vertices as not visited
            std::vector<EdgeIterator>next_out = std::vector<EdgeIterator>(get_n_nodes()+1);
            for(int i = 0; i < _nodes.size(); i++){
                next_out[i] = _nodes[i]._first_edge;
            }
            std::vector<bool> was_pushed(get_n_nodes(), false);
            std::vector<bool> was_popped(get_n_nodes(), false);
            std::stack<NodeIterator>to_visit;
            NodeIterator next_id = 0;
            for(int s = 0; s< get_n_nodes(); ++s){
                NodeIterator source_node = ch_order_to_node[s];
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
                            assert(next_out[x]+1 < _edges.size());
                            NodeIterator y = next_node[next_out[x]];
                            // upward edge only;
                            if(ch_ordering[y] < ch_ordering[x]){
                                ++next_out[x];
                                continue;
                            }
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
                            assert(next_out[x]+1 < _edges.size());
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



        void  load_coordinate(const std::string& coordinate_file){
            std::ifstream inNodeFile(coordinate_file);
            if(!inNodeFile)
            {
                std::cout << "Cannot open coordinate file" << coordinate_file<< std::endl;
                return ;
            }
            std::cout << "Reading coordinate file: " << coordinate_file<< std::endl;

            int node_num;
            inNodeFile >> node_num >> min_lat >> max_lat >> min_lon >> max_lon;
            _coordinate.resize(node_num);
            double id;
            double x,y;
            for(int i = 0; i < node_num; i++) {
                inNodeFile >> id >> x >> y;
                _coordinate[i] =  point_type(y,x);
            }
        }

        double get_Euclidean_travel_time(NodeIterator start, NodeIterator target){
            // y is Latitude and x is Longitude;
            //88.5139 is the max speed in the map
            // time unit is 0.1 second, so convert from hours to 0.1 second is 60*60*10;
            double Euclidean_distance_in_km  = util::geo_dist(_coordinate[start].y(), _coordinate[start].x(),
                                                       _coordinate[target].y(), _coordinate[target].x()) /1000;
            return Euclidean_distance_in_km / 88.5139 * 60 * 60 * 10;
//            return 0;
        }

        void check_the_graph(){
            auto num_nodes = get_n_nodes();
            for(unsigned i  =0 ; i < num_nodes ; i++){
                for(EdgeIterator e = edges_begin(i); e != edges_end(i); e++){
                    NodeIterator other_node = get_other_node(e);
                    for(EdgeIterator e1 = edges_begin(i); e1 != edges_end(i); e1++) {
                        if(e != e1){
                            if(other_node == get_other_node(e1)){
                                std::cout<<"redundant edges"<<std::endl;
                            }
                        }
                    }
                }

            }
            std::cout<<"test passed"<<std::endl;


        }

        void export_wkt(const std::vector<NodeIterator>& node_id, const std::string& output_file){
            std::ofstream myFile(output_file);
            myFile<<"wkt\n";
            for(auto& id : node_id){
                myFile<< geom::wkt(_coordinate[id]) << "\n";
            }
            myFile.close();
        }


        // Some useless code, maybe useful later !! through it here !

        void mark_level_down(std::vector<int>& level, const NodeIterator root_node ){

            std::stack<std::pair<NodeIterator,int>>to_visit;
            std::vector<bool> was_popped(get_n_nodes(), false);
            level[root_node] = 0;
            to_visit.push({root_node,0});
            while(!to_visit.empty()){
                auto x = to_visit.top();
                to_visit.pop();


                if(level[x.first] >= x.second){
                    continue;
                }

                if(!was_popped[x.first]){
                    was_popped[x.first] = true;
                    level[x.first] = std::max(level[x.first],x.second);
                }else{
                    if(level[x.first] >= x.second){
                        continue;
                    }else{
                        level[x.first] = x.second;
                    }
                }


                for ( EdgeIterator e = edges_begin(x.first) ; e != edges_end(x.first) ; ++e ) {
                    if(is_directed_downward(e)){
                        NodeIterator other_node = get_other_node(e);
                        if(!was_popped[other_node]){
                            was_popped[other_node] = true;
                            to_visit.push({other_node,level[x.first]+1});
                        }else{
                            if(level[other_node] < level[x.first]+1 ){
                                to_visit.push({other_node,level[x.first]+1});
                            }
                        }
                    }
                }

            }
        }



        void mark_level_down2(std::vector<int>& level, const NodeIterator root_node ){
            std::cout<<"here"<<std::endl;
            std::stack<std::pair<NodeIterator,int>>to_visit;
            std::vector<bool> was_popped(get_n_nodes(), false);
            level[root_node] = 0;
            to_visit.push({root_node,0});
            while(!to_visit.empty()){
                auto x = to_visit.top();
                to_visit.pop();
                level[x.first] = std::max(x.second,level[x.first]);
                for ( EdgeIterator e = edges_begin(x.first) ; e != edges_end(x.first) ; ++e ) {
                    if(is_directed_downward(e)){
                        NodeIterator other_node = get_other_node(e);
                        to_visit.push({other_node,level[x.first]+1});
                    }
                }

            }
        }

        std::vector<int> get_node_level(){
            unsigned num_nodes = get_n_nodes();
            std::vector<int> level = std::vector<int>(num_nodes,-1);
            //mark root nodes first ;
            for(NodeIterator i = 0; i < num_nodes; i++){
                bool is_toppest_node = true;
                for ( EdgeIterator e = edges_begin(i) ; e != edges_end(i) ; ++e ) {
                    if(is_directed_upward(e)){
                        // contains upward nodes !
                        is_toppest_node = false;
                        break;
                    }
                }
                if(is_toppest_node){
                    mark_level_down(level, i );
                }
            }
            // mark rest of nodes ;

            for(NodeIterator i = 0; i < num_nodes; i++){
                if(level[i] == -1){
                    for(NodeIterator j = 0; j < num_nodes; j++) {
                        if(is_down_reachable(j,  i)){
                            std::cout<< "Ordering Error"<<std::endl;
                            std::cout << "downward reachable " << j << " "<< i<<std::endl;
                        }
                    }



                }
                bool is_toppest_node = true;
                for ( EdgeIterator e = edges_begin(i) ; e != edges_end(i) ; ++e ) {
                    if(is_directed_upward(e)){
                        // contains upward nodes !
                        is_toppest_node = false;
                        break;
                    }
                }
                if(is_toppest_node){
                    mark_level_down(level, i );
                }
            }






            for(;;){
                int num_of_unmarked_nodes = 0;
                for(NodeIterator i = 0; i < num_nodes; i++){
                    if(level[i] == -1){
                        num_of_unmarked_nodes++;
                        bool is_toppest_node = true;
                        for ( EdgeIterator e = edges_begin(i) ; e != edges_end(i) ; ++e ) {
                            if(is_directed_upward(e)){
                                // contains upward nodes !
                                is_toppest_node = false;
                                break;
                            }
                        }
                        if(is_toppest_node){
                            mark_level_down(level, i );
                        }
                    }
                }
                std::cout<<num_of_unmarked_nodes<<std::endl;
                if(num_of_unmarked_nodes == 0 ){
                    break;
                }
            }
            return level;
        }

        void test_ordering(const std::vector<unsigned int>& level){

            for(int i = 0; i < get_n_nodes(); i++){
                for ( EdgeIterator e = edges_begin(i) ; e != edges_end(i) ; ++e ) {
                    auto other = get_other_node(e);
                    if(is_directed_downward(e)){
                        // downward other < i
                        if(level[other]>= level[i]){
                            std::cout<<"error"<<std::endl;
                        }
                    }else{
                        // upward  other > i
                        if(level[other]<= level[i]){
                            std::cout<<"error"<<std::endl;
                        }
                    }
                }

            }
            std::cout<<"Passed"<<std::endl;

        }
    };




}


#endif //KATCH_FORWARD_SEARCH_GRAPH_H
