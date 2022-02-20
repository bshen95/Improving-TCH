//
// Created by Bojie Shen on 21/7/21.
//

#ifndef KATCH_FTCH_L_H
#define KATCH_FTCH_L_H
#include <array>
#include <cassert>
#include <cmath>
#include <deque>
#include <limits>
#include <stack>
#include <tuple>
#include <vector>
#include <iostream>
#include <datastr/cpd/cpd.h>
#include <datastr/cpd/bw_cpd.h>

//#include <boost/heap/pairing_heap.hpp>
#include "util/id_queue.h"
#include "datastr/base/double.h"
#include "datastr/base/interval.h"
#include "datastr/graph/basic.h"
#include "datastr/graph/forward_search_graph.h"
#include "context/search_context.h"

#define INF std::numeric_limits<double>::max()

namespace katch
{

    class FTCH_L
    {

    private:
        using SearchNodeId = uint32_t;
        static constexpr double period = 864000.0;
        struct Predecessor{
            EdgeIterator _edge_it;
            SearchNodeId _search_node_id;

            Predecessor()
                    : _edge_it(INVALID_EDGE_ITERATOR), _search_node_id(std::numeric_limits<SearchNodeId>::max()) {}

            Predecessor(const EdgeIterator &e, const SearchNodeId id)
                    : _edge_it(e), _search_node_id(id) {}
        };

        struct SearchNode_tmpl {
            NodeIterator _node_it;
            bool _enqueued;
            double _t_arr;
            bool _is_downward;
            double _heuristic;
            Predecessor _predecessor;

            SearchNode_tmpl(const SearchNode_tmpl &) = delete;

            SearchNode_tmpl &operator=(const SearchNode_tmpl &) = delete;

            SearchNode_tmpl(SearchNode_tmpl &&) = default;

            SearchNode_tmpl &operator=(SearchNode_tmpl &&) = default;

            SearchNode_tmpl()
                    :
                    _node_it(INVALID_NODE_ITERATOR),
                    _enqueued(uint32_t(0)),
                    _t_arr(std::numeric_limits<double>::max()),
                    _is_downward(false),
                    _heuristic(std::numeric_limits<double>::max()),
                    _predecessor() {}
        };

    public:

        using SearchNode = SearchNode_tmpl;

    private:

        ForwardSearchGraph<&period> _graph;

        NodeIterator _start;
        NodeIterator _destination;
        double _t_dep;
        double _t_arr;
        SearchContext<SearchNode_tmpl,SearchNodeId> _context;
        double _upper_bound;
        BW_CPD _bw_cpd;
        vector<unsigned> _cpd_mapper;

        unsigned short _search_id;

        vector<unsigned short> _was_bw_cached;
        vector<bool> _was_bw_reachable;
        vector<unsigned short> _was_popped;
        unsigned _mapped_destination;
        vector<double> _landmarks;
        unsigned _number_of_landmark;

        bool is_down_reachable(const unsigned start){
            if(_was_bw_cached[start] == _search_id){
                return _was_bw_reachable[start];
            }else{
                _was_bw_cached[start] = _search_id;
                _was_bw_reachable[start] = _bw_cpd.get_reachability(start,_mapped_destination);
                number_of_bw_first_move_calls++;
                return _was_bw_reachable[start];
            }

        }

        double get_lower_bound(unsigned start, unsigned destination){
//            return _graph.get_Euclidean_travel_time(start,destination);
//            return 0;
            if(start == destination){
                return 0 ;
            }
            double max = 0 ;
            auto startPtr = _landmarks.begin()+ start*  _number_of_landmark *2;
            auto endPtr =_landmarks.begin()+ destination*  _number_of_landmark *2;
            for(int i = 0 ; i < _number_of_landmark; i++){
                if(*(startPtr) ==INF || *(endPtr) == INF){
                    continue ;
                }
                double landmark_distance = std::max((*(endPtr) -*(startPtr) ) ,(*(startPtr +1) - *(endPtr +1))   );
                if(max < landmark_distance){
                    max = landmark_distance;
                }
                startPtr +=2;
                endPtr+=2;
            }
            return max;
        }


        double search(){
            _search_id++;
            _was_bw_reachable[_destination] = true;
            _was_bw_cached[_destination] = _search_id;

            double h = get_lower_bound(_start,_destination);
            SearchNode& s = _context.insert(_start, _t_dep + h);
            number_of_nodes_generated++;
            s._t_arr = _t_dep;
            s._is_downward = false;
            s._heuristic = h;
            while ( ! _context.pq_empty() )
            {

                SearchNode& u = _context.pq_delete_min();

                number_of_nodes_expanded++;
                const SearchNodeId u_id = _context.get_search_node_id(u);
                const unsigned u_it = u._node_it;
                const double t_arr_u = u._t_arr;
                const bool is_downward = u._is_downward;
                const bool downward_reachable = is_down_reachable(u_it);

//                std::cout<<"Popping nodes: "<<"Node Id: "<<u_it<<" Arrival time: "<<t_arr_u
//                <<" Heuristic: "<<heuristic <<" Downward ? "<<is_downward<< std::endl;

                if(u_it == _destination){
                    break;
                };

                if(downward_reachable){
                    // downward_successors;
                    for ( EdgeIterator e = _graph.downward_begin(u_it) ; e != _graph.downward_end(u_it) ; ++e ) {
                        const NodeIterator v_it = _graph.get_other_node(e);
                        if (!is_down_reachable(v_it)) continue;
                        if ( ! _context.reached(v_it) ) {
                            const double t_arr_v_new = _graph.get_ttf(e).eval(t_arr_u)+ t_arr_u;
                            number_of_TD_calculation ++;
                            h  = get_lower_bound(v_it,_destination);
                            SearchNode& v = _context.insert(v_it, t_arr_v_new + h );
                            v._t_arr = t_arr_v_new;
                            v._predecessor = Predecessor(e, u_id);
                            v._is_downward = true;
                            v._heuristic = h;
                            number_of_nodes_generated++;
//                        std::cout<<"Generating nodes: "<<"Node Id: "<<v_it<<" Arrival time: "<<v._t_arr
//                                 <<" Heuristic: "<<v._heuristic<<" Downward ? "<< v._is_downward<< std::endl;
                        }else{
                            assert( _context.reached(v_it) );
                            SearchNode& v = _context.get_search_node(v_it);
                            if(t_arr_u + _graph.get_ttf(e).get_min() >= v._t_arr ){
                                continue;
                            }
                            const double t_arr_v_new = _graph.get_ttf(e).eval(t_arr_u)+ t_arr_u;
                            number_of_TD_calculation ++;
                            if ( t_arr_v_new >= v._t_arr )
                                continue;
                            v._t_arr = t_arr_v_new;
                            v._predecessor = Predecessor(e, u_id);
                            v._is_downward = true;
                            //CH-CPD is a inconsistent heuristic, because the lower bound of shortcuts is larger than the lower bound of original arcs.
                            if ( _context.pq_contains(v) )
                                _context.pq_decrease(v, t_arr_v_new + v._heuristic);
                            else
                                _context.pq_re_insert(v, t_arr_v_new + v._heuristic);
//                        std::cout<<"Generating nodes: "<<"Node Id: "<<v_it<<" Arrival time: "<<v._t_arr
//                                 <<" Heuristic: "<<v._heuristic<<" Downward ? "<< v._is_downward<< std::endl;
                            number_of_nodes_generated++;
                        }
                    }
                }
                if (!is_downward){
                    //upward successor;
                    for ( EdgeIterator e = _graph.upward_begin(u_it) ; e != _graph.upward_end(u_it) ; ++e ) {
                        const NodeIterator v_it = _graph.get_other_node(e);
                        if ( ! _context.reached(v_it) ) {
                            const double t_arr_v_new = _graph.get_ttf(e).eval(t_arr_u)+ t_arr_u;
                            number_of_TD_calculation ++;
                            h  = get_lower_bound(v_it,_destination);
                            SearchNode& v = _context.insert(v_it, t_arr_v_new + h );
                            v._t_arr = t_arr_v_new;
                            v._predecessor = Predecessor(e, u_id);
                            v._is_downward = false ;
                            v._heuristic = h;
                            number_of_nodes_generated++;
//                        std::cout<<"Generating nodes: "<<"Node Id: "<<v_it<<" Arrival time: "<<v._t_arr
//                                 <<" Heuristic: "<<v._heuristic<<" Downward ? "<< v._is_downward<< std::endl;
                        }else{
                            assert( _context.reached(v_it) );
                            SearchNode& v = _context.get_search_node(v_it);
                            if(t_arr_u + _graph.get_ttf(e).get_min() >= v._t_arr ){
                                continue;
                            }
                            const double t_arr_v_new = _graph.get_ttf(e).eval(t_arr_u)+ t_arr_u;
                            number_of_TD_calculation ++;
                            if ( t_arr_v_new >= v._t_arr )
                                continue;
                            v._t_arr = t_arr_v_new;
                            v._predecessor = Predecessor(e, u_id);
                            v._is_downward = false;
                            if ( _context.pq_contains(v) )
                                _context.pq_decrease(v, t_arr_v_new + v._heuristic);
                            else
                                _context.pq_re_insert(v, t_arr_v_new + v._heuristic);
//                        std::cout<<"Generating nodes: "<<"Node Id: "<<v_it<<" Arrival time: "<<v._t_arr
//                                 <<" Heuristic: "<<v._heuristic<<" Downward ? "<< v._is_downward<< std::endl;
                            number_of_nodes_generated++;
                        }
                    }
                }
            }


            assert( _context.reached(_destination) );
            const SearchNode& d = _context.get_search_node(_destination);
            return d._t_arr;

        }


//        double search(){
//            _search_id++;
//            _was_bw_reachable[_destination] = true;
//            _was_bw_cached[_destination] = _search_id;
//            double h = get_lower_bound(_start,_destination);
//            SearchNode& s = _context.insert(_start, _t_dep + h);
//            number_of_nodes_generated++;
//            s._t_arr = _t_dep;
//            s._is_downward = false;
//            s._heuristic = h;
//
//            while ( ! _context.pq_empty() )
//            {
//
//                SearchNode& u = _context.pq_delete_min();
//
//                number_of_nodes_expanded++;
//                const SearchNodeId u_id = _context.get_search_node_id(u);
//                const NodeIterator u_it = u._node_it;
//                const double t_arr_u = u._t_arr;
//                const bool is_downward = u._is_downward;
//                const double heuristic = u._heuristic;
//                const bool downward_reachable = is_down_reachable(u_it);
////                std::cout<<"Popping nodes: "<<"Node Id: "<<u_it<<" Arrival time: "<<t_arr_u
////                <<" Heuristic: "<<heuristic <<" Downward ? "<<is_downward<< std::endl;
//
//                // always take value here, because emplace back may reallocate the memory;
//                if(u_it == _destination){
//                    break;
//                };
//                for ( EdgeIterator e = _graph.edges_begin(u_it) ; e != _graph.edges_end(u_it) ; ++e )
//                {
//                    // search go down can never go up.
//                    const NodeIterator v_it = _graph.get_other_node(e);
//                    if(_graph.is_directed_upward(e)){
//                        // upward successor;
//                        if(is_downward) continue;
//
//                    }else {
//                        // downward successor
//                        if (!downward_reachable) continue;
//                        if (!is_down_reachable(v_it)) continue;
//                    }
//
//                    if ( ! _context.reached(v_it) ) {
//                        h  = get_lower_bound(v_it,_destination);
//                        const double t_arr_v_new = _graph.get_ttf(e).eval(t_arr_u)+ t_arr_u;
//                        SearchNode& v = _context.insert(v_it, t_arr_v_new + h );
//                        v._t_arr = t_arr_v_new;
//                        v._predecessor = Predecessor(e, u_id);
//                        v._is_downward = _graph.is_directed_downward(e);
//                        v._heuristic = h;
//                        number_of_nodes_generated++;
//
////                        std::cout<<"Generating nodes: "<<"Node Id: "<<v_it<<" Arrival time: "<<v._t_arr
////                                 <<" Heuristic: "<<v._heuristic<<" Downward ? "<< v._is_downward<< std::endl;
//                    }else{
//                        assert( _context.reached(v_it) );
//                        SearchNode& v = _context.get_search_node(v_it);
//                        const double t_arr_v_new = _graph.get_ttf(e).eval(t_arr_u)+ t_arr_u;
//                        if ( t_arr_v_new >= v._t_arr )
//                            continue;
//
//                        v._t_arr = t_arr_v_new;
//                        v._predecessor = Predecessor(e, u_id);
//                        v._is_downward = _graph.is_directed_downward(e);
//                        if ( _context.pq_contains(v) )
//                            _context.pq_decrease(v, t_arr_v_new + v._heuristic);
//                        else
//                            _context.pq_re_insert(v, t_arr_v_new + v._heuristic);
////                        std::cout<<"Generating nodes: "<<"Node Id: "<<v_it<<" Arrival time: "<<v._t_arr
////                                 <<" Heuristic: "<<v._heuristic<<" Downward ? "<< v._is_downward<< std::endl;
//                        number_of_nodes_generated++;
//                    }
//
//                }
//            }
//
//            assert( _context.reached(_destination) );
//            const SearchNode& d = _context.get_search_node(_destination);
////            std::cout<<d._t_arr<<std::endl;
//            return d._t_arr;
//        }

    public:

        unsigned long number_of_nodes_generated;
        unsigned long number_of_nodes_expanded;
        unsigned long number_of_first_move_calls;
        unsigned long number_of_nodes_visited;
        unsigned long number_of_TD_calculation;
        unsigned long number_of_bw_first_move_calls;

        using Graph = ForwardSearchGraph<&period> ;

        FTCH_L(ForwardSearchGraph<&period> && graph)
                : _graph(std::move(graph)),
                  _start(INVALID_NODE_ITERATOR),
                  _destination(INVALID_NODE_ITERATOR),
                  _t_dep(std::numeric_limits<double>::max()),
                  _t_arr(std::numeric_limits<double>::max()),
                  _context(_graph.get_n_nodes()),
                  _upper_bound(std::numeric_limits<double>::max()),
                  _was_popped(_graph.get_n_nodes()),
                  _was_bw_cached(_graph.get_n_nodes()),
                  _was_bw_reachable(_graph.get_n_nodes()),
                  _search_id(0)
        {}

        void load_cpd(const std::string& filename){
            string name = (filename.substr(0,filename.find_last_of(".")));
            _cpd_mapper = load_vector<unsigned>(name.substr(0,name.find_last_of("_"))+".bw_tch_mapper_ch_order");
//            _cpd_free_flow = load_vector<double>(filename.substr(0,filename.find_last_of(".")) + ".fw_tch_free_flow");
            _bw_cpd = BW_CPD();
            FILE*f = fopen(filename.c_str(), "r");
            _bw_cpd.load(f);

        }


        void load_landmark(const std::string& filename,unsigned number_of_landmark){
            _landmarks.clear();
            _landmarks = load_vector<double>(filename);
            _number_of_landmark = number_of_landmark;
            if(_landmarks.size() != _number_of_landmark * _graph.get_n_nodes()*2){
                std::cout<<"Loading landmarks failed" << std::endl;
            }
        }


        double get_index_size(){
            double bw_cpd_size  = (double)_bw_cpd.get_entry().size()*4/1000000;
            std::cout<<"Backward CPD size: "<<bw_cpd_size<<" MB"<<std::endl;
            double total_size =  bw_cpd_size;
            return total_size  ;
        }

        void load_target_row(const double t_dep){
            // only used in RTPD
            return;
        }

        double one_to_one(const NodeIterator start, const NodeIterator destination, const double t_dep)
        {
            number_of_first_move_calls = 0;
            number_of_nodes_generated = 0;
            number_of_nodes_expanded = 0;
            number_of_TD_calculation = 0;
            number_of_bw_first_move_calls = 0;

            _t_arr = std::numeric_limits<double>::max();

            _start = start;
            _destination = destination;
            _t_dep = t_dep;
            _mapped_destination = _cpd_mapper[_destination];
            if ( start == destination )
                return _t_dep;

            _upper_bound = std::numeric_limits<double>::max();
            _context.clear_all();

            _t_arr = search();
            return _t_arr;
        }

        class Path
        {
            friend FTCH_L;

        private:

            std::deque<NodeIterator> _nodes;
            std::deque<EdgeIterator> _edges;
            std::deque<double> _t_arr;

            void append_front(const NodeIterator u, const double t_arr_u, const EdgeIterator e)
            {
                _nodes.push_front(u);
                _edges.push_front(e);
                _t_arr.push_front(t_arr_u);
            }

            void append_back(const EdgeIterator e, const NodeIterator u, const double t_arr_u)
            {
                _nodes.push_back(u);
                _edges.push_back(e);
                _t_arr.push_back(t_arr_u);
            }

            bool dbg_check_path(const ForwardSearchGraph<&period>& graph, const double t_dep)
            {
                double t_current = t_dep;

                if ( _nodes.size() == 0 ) return false;
                if ( _nodes.size() != _t_arr.size() ) return false;
                if ( _nodes.size() != _edges.size() + 1 ) return false;

                if ( fabs(t_current - _t_arr.front()) > 0.00001 ) return false;

                for ( size_t  i = 1 ; i < _nodes.size() ; ++i )
                {
                    const EdgeIterator e = _edges[i-1];
                    t_current += graph.get_ttf(e).eval(t_current);

                    if ( fabs(t_current - _t_arr[i]) > 0.00001 ) return false;
                }

                return true;
            }


            bool dbg_check_path_shortcut_free(const ForwardSearchGraph<&period> & graph, const double t_dep)
            {
                double t_current = t_dep;

                for ( size_t  i = 1 ; i < _nodes.size() ; ++i )
                {
                    const EdgeIterator e = _edges[i-1];
                    if ( graph.get_middle_node(e, t_current) != INVALID_NODE_ITERATOR ) return false;

                    t_current += graph.get_ttf(e).eval(t_current);
                }

                return true;
            }

        public:

            Path(const NodeIterator destination, const double t_arr)
                    : _nodes(1, destination), _edges(), _t_arr(1, t_arr)
            {}

            size_t get_n_edges() const noexcept
            {
                assert( _nodes.size() == _edges.size() + 1 );
                return _edges.size();
            }

            std::tuple<NodeIterator,EdgeIterator,NodeIterator> get_edge(const size_t index) const noexcept
            {
                assert( _nodes.size() == _edges.size() + 1 );
                assert( index < _edges.size() );

                NodeIterator u = _nodes[index];
                EdgeIterator e = _edges[index];
                NodeIterator v = _nodes[index+1];

                return std::make_tuple(u, e, v);
            }

            double get_t_arr(const size_t index) const noexcept
            {
                assert( _nodes.size() == _t_arr.size() );
                assert( index < _t_arr.size() );

                return _t_arr[index];
            }
            void print_nodes() const {

                for(int i =0; i< _nodes.size(); i++){
                    std::cout<<_nodes[i]<<" time: "<<_t_arr[i]<<std::endl;
                }
                double cost = 0;
                for(int i = 0 ; i < _edges.size(); i++){

                    std::cout<<" edge id: "<< _edges[i]<< std::endl;
                }
            }

            void print_nodes(const ForwardSearchGraph<&period> & graph) const {

                for(int i =0; i< _nodes.size(); i++){
                    std::cout<<_nodes[i]<<" time: "<<_t_arr[i]<<std::endl;
                }
                double cost = 0;
                for(int i = 0 ; i < _edges.size(); i++){

                    std::cout<<" edge id: "<< _edges[i]<< "Direction" << graph.is_directed_upward(_edges[i]) << std::endl;
                }
            }
        };

        Path get_up_down_path(NodeIterator node) const noexcept
        {
            if ( ! _context.reached(_start) ) return Path(_start, _t_dep);
            if ( ! _context.reached(node) ) return Path(_start, _t_dep);

            const SearchNode& d = _context.get_search_node(node);
            const SearchNodeId d_id = _context.get_search_node_id(d);

            const SearchNode& s = _context.get_search_node(_start);
            const SearchNodeId s_id = _context.get_search_node_id(s);

            Path result(d._node_it, d._t_arr);

            SearchNodeId u_id = d_id;

            while ( u_id != s_id )
            {
                const SearchNode& u_prev = _context.get_search_node_from_id(u_id);

                const EdgeIterator e = u_prev._predecessor._edge_it;
                assert ( e != INVALID_EDGE_ITERATOR );

                u_id = u_prev._predecessor._search_node_id;
                const SearchNode& u = _context.get_search_node_from_id(u_id);

                assert( eq( _graph.get_ttf(e).eval(u._t_arr) + u._t_arr, u_prev._t_arr ) );
                result.append_front(u._node_it, u._t_arr, e);
            }

            assert( result._nodes.front() == _start );
            assert( result._nodes.back() == _destination );
            assert( result._t_arr.front() == _t_dep );
            assert( fabs(result._t_arr.back() - _t_arr) <= 0.00001 );
            assert( result.dbg_check_path(_graph, _t_dep) );

            return result;
        }

        Path expand_path(const Path& up_down_path) const noexcept
        {
            Path result(_start, _t_dep);
            double t_current = _t_dep;

            for ( size_t i = 0 ; i != up_down_path.get_n_edges() ; ++i )
            {
                assert( result._nodes.back() == up_down_path._nodes[i] );
                assert( fabs(result._t_arr.back() - up_down_path.get_t_arr(i)) <= 0.00001 );

                std::stack< std::tuple<NodeIterator,EdgeIterator,NodeIterator> > S;
                S.push( up_down_path.get_edge(i) );

                while ( ! S.empty() )
                {
                    NodeIterator u = std::get<0>(S.top());
                    EdgeIterator e = std::get<1>(S.top());
                    NodeIterator v = std::get<2>(S.top());
                    S.pop();

                    NodeIterator x = _graph.get_middle_node(e, t_current);

                    if ( x == INVALID_NODE_ITERATOR )
                    {
                        t_current += _graph.get_ttf(e).eval(t_current);
                        result.append_back(e, v, t_current);
                        continue;
                    }

                    EdgeIterator e_up = INVALID_EDGE_ITERATOR;
                    EdgeIterator e_down = INVALID_EDGE_ITERATOR;

                    for ( EdgeIterator it = _graph.edges_begin(x) ; it != _graph.edges_end(x) ; ++it )
                    {
                        if ( _graph.get_other_node(it) == v && _graph.is_directed_upward(it) ) e_up = it;
                        if ( e_up != INVALID_EDGE_ITERATOR ) break;
                    }
                    for ( EdgeIterator it = _graph.edges_begin(u) ; it != _graph.edges_end(u) ; ++it )
                    {
                        if ( _graph.get_other_node(it) == x && _graph.is_directed_downward(it) ) e_down = it;
                        if ( e_down != INVALID_EDGE_ITERATOR) break;
                    }

                    S.push( std::make_tuple(x, e_up, v) );
                    S.push( std::make_tuple(u, e_down, x) );
                }
            }

            assert( result._nodes.front() == _start );
            assert( result._nodes.back() == _destination );
            assert( result._t_arr.front() == _t_dep );
            assert( fabs(result._t_arr.back() - _t_arr) <= 0.00001 );
            assert( result.dbg_check_path(_graph, _t_dep) );
            assert( result.dbg_check_path_shortcut_free(_graph, _t_dep) );

            return result;
        }

        Path get_path(NodeIterator node) const noexcept{
            const auto up_down_path = get_up_down_path(node);
//            up_down_path.print_nodes(_graph);
            return  expand_path( up_down_path)  ;
        }
    };
}
#endif //KATCH_FTCH_L_H
