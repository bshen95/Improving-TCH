//
// Created by Bojie Shen on 3/10/21.
//

#ifndef KATCH_BMTCH_H
#define KATCH_BMTCH_H
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
#include "query/BSTCH.h"
#include "BMTCH_query.h"

#define INF 9999999;

namespace katch
{

    class BMTCH
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


        CPD _cpd;
        vector<unsigned> _cpd_mapper;
        vector<unsigned> _cpd_bw_mapper;



        //times_flag;
        unsigned _search_id;
        vector<unsigned> _was_cached;
        vector<unsigned> _was_bw_cached;

        vector<double> _lower_bound_cache;
        vector<bool> _was_bw_reachable;
        unsigned _mapped_destination;
        unsigned _bw_mapped_destination;
        unsigned _destination_ranking;

        vector<BMTCHQuery_Wrapper*> hourly_query;
        unsigned hourly_query_index;



        double get_upper_bound(const unsigned& start, double departure_time){
            double arrival_time = departure_time;
            unsigned current_node = start;
            for(;;){
                const unsigned char& fm = _cpd.get_first_move(current_node, _mapped_destination);
                number_of_first_move_calls++;
                const unsigned& next_id = _graph.get_other_node( (EdgeIterator)(_graph.edges_begin(current_node) + fm));
                // taking the actual path wont make hourly different. so max should be fine.
                arrival_time = arrival_time + _graph.get_ttf((EdgeIterator)(_graph.edges_begin(current_node) + fm)).get_max();
                number_of_TD_calculation ++;
                if(next_id == _destination){
                    break;
                }else{
                    current_node = next_id;
                }
            }
            return arrival_time -departure_time;
        }

        double get_lower_bound(const unsigned& start){
            // recursion
            if(_was_cached[start] == _search_id){
                return _lower_bound_cache[start];
            }
            const unsigned char& fm = _cpd.get_first_move(start, _mapped_destination);
            number_of_first_move_calls++;
            if(fm == 0XFF){
                _was_cached[start] = _search_id;
                _lower_bound_cache[start] = INF;
                return INF;
            }
            const unsigned& next_id = _graph.get_other_node( (EdgeIterator)(_graph.edges_begin(start) + fm));
            _lower_bound_cache[start] = get_lower_bound(next_id) + _graph.get_ttf((EdgeIterator)(_graph.edges_begin(start) + fm)).get_min();
            _was_cached[start] = _search_id;
            return _lower_bound_cache[start];
        }


        bool is_down_reachable(const unsigned start){
            if(_was_bw_cached[start] == _search_id){
                return _was_bw_reachable[start];
            }else{
                _was_bw_cached[start] = _search_id;
                _was_bw_reachable[start] = _bw_cpd.get_reachability(start, _bw_mapped_destination);
                number_of_bw_first_move_calls ++;
                return _was_bw_reachable[start];
            }

        }

        double search(){
            _search_id++;
            _was_bw_reachable[_destination] = true;
            _was_bw_cached[_destination] = _search_id;

            _lower_bound_cache[_destination] = 0;
            _was_cached[_destination] = _search_id;

            double h = get_lower_bound(_start);
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

                // always take value here, because emplace back may reallocate the memory;
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
                            h  = get_lower_bound(v_it);
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
                            h  = get_lower_bound(v_it);
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
            }


            assert( _context.reached(_destination) );
            const SearchNode& d = _context.get_search_node(_destination);
            return d._t_arr;

        }

    public:

        unsigned long number_of_nodes_generated;
        unsigned long number_of_nodes_expanded;
        unsigned long number_of_nodes_visited;
        unsigned long number_of_TD_calculation;
        unsigned long number_of_first_move_calls;
        unsigned long number_of_bw_first_move_calls;

        unsigned long query_answered_by_1_hour_query;
        unsigned long query_answered_by_2_hour_query;
        unsigned long query_answered_by_4_hour_query;
        unsigned long query_answered_by_tch_query;


        using Graph = ForwardSearchGraph<&period> ;

        BMTCH(ForwardSearchGraph<&period> && graph)
                : _graph(std::move(graph)),
                  _start(INVALID_NODE_ITERATOR),
                  _destination(INVALID_NODE_ITERATOR),
                  _t_dep(std::numeric_limits<double>::max()),
                  _t_arr(std::numeric_limits<double>::max()),
                  _context(_graph.get_n_nodes()),
                  _upper_bound(std::numeric_limits<double>::max()),
                  _was_cached(_graph.get_n_nodes()),
                  _lower_bound_cache(_graph.get_n_nodes()),
                  _was_bw_cached(_graph.get_n_nodes()),
                  _was_bw_reachable(_graph.get_n_nodes()),
                  _search_id(0)
        {
            query_answered_by_1_hour_query = 0;
            query_answered_by_2_hour_query = 0;
            query_answered_by_4_hour_query = 0;
            query_answered_by_tch_query = 0;

        }

        void set_hourly_queries(vector<BMTCHQuery_Wrapper*> h_queries){
            hourly_query  = h_queries;
        }

        void load_cpd(const std::string& filename){
            std::cout<<filename << std::endl;
            _cpd = CPD();
            FILE*f = fopen(filename.c_str(), "r");
            _cpd.load(f);
            string name = (filename.substr(0,filename.find_last_of(".")));
            _cpd_mapper = load_vector<unsigned>(name.substr(0,name.find_last_of("_"))+".fw_tch_mapper");
            _cpd_bw_mapper = load_vector<unsigned>(name.substr(0,name.find_last_of("_"))+".bw_tch_mapper_ch_order");
            _bw_cpd = BW_CPD();
            std::cout<<name+".bw_tch_cpd_ch_order" << std::endl;
            FILE*f1 = fopen((name+".bw_tch_cpd_ch_order").c_str(), "r");
            _bw_cpd.load(f1);
        }



        double get_index_size(){
            double fw_cpd_size  = (double)_cpd.get_entry().size()*4/1000000;
            std::cout<<"Forward CPD size: "<<fw_cpd_size<<" MB"<<std::endl;

            double bw_cpd_size  = (double)_bw_cpd.get_entry().size()*4/1000000;
            std::cout<<"Backward CPD size: "<<bw_cpd_size<<" MB"<<std::endl;

            double total_size =  fw_cpd_size + bw_cpd_size;
            return total_size  ;
        }

        bool fit_into_hourly_query(const double t_dep, const double upper_bound){
            for(unsigned i = 0; i < hourly_query.size(); i ++){
                if(hourly_query[i]->can_ts_tch_query_answer_the_query(t_dep,upper_bound)){
                    hourly_query_index = i;
                    return true;
                }
            }
            return false;
        }
        void load_target_row(const double t_dep){
            // only used in RTPD
            return;
        }
        double one_to_one(const NodeIterator start, const NodeIterator destination, const double t_dep)
        {
            number_of_first_move_calls = 0;
            number_of_bw_first_move_calls = 0;
            number_of_nodes_generated = 0;
            number_of_nodes_expanded = 0;
            number_of_TD_calculation = 0;
            hourly_query_index = 9999;
            _start = start;
            _destination = destination;
            _t_dep = t_dep;
            _mapped_destination = _cpd_mapper[_destination];

            if ( start == destination )
                return _t_dep;

            double upper_bound  = get_upper_bound(_start,t_dep);

            if(fit_into_hourly_query(t_dep,upper_bound)){
                if(hourly_query_index == 0){
                    query_answered_by_1_hour_query ++;
                }else if (hourly_query_index == 1 ){
                    query_answered_by_2_hour_query ++;
                }else{
                    query_answered_by_4_hour_query ++;
                }
                _t_arr = hourly_query[hourly_query_index]->one_to_one(start, destination,t_dep);
                number_of_first_move_calls += hourly_query[hourly_query_index]->number_of_first_move_calls;
                number_of_bw_first_move_calls = hourly_query[hourly_query_index]->number_of_bw_first_move_calls;
                number_of_nodes_generated = hourly_query[hourly_query_index]-> number_of_nodes_generated;
                number_of_nodes_expanded = hourly_query[hourly_query_index]-> number_of_nodes_expanded;
                number_of_TD_calculation = hourly_query[hourly_query_index]-> number_of_TD_calculation;
                return _t_arr;
            }else{
                query_answered_by_tch_query++;
                _t_arr = std::numeric_limits<double>::max();
                _bw_mapped_destination = _cpd_bw_mapper[_destination];
                _context.clear_all();
                _t_arr = search();
                return _t_arr;
            }
        }

        void check_up_down_path(NodeIterator node)
        {
            vector<unsigned> path;
            vector<unsigned> pre_path;
            if ( ! _context.reached(_start) ) {
                path.push_back(_start);
            }else if ( ! _context.reached(node) ) {
                path.push_back(_start);
            }else{
                const SearchNode& d = _context.get_search_node(node);
                const SearchNodeId d_id = _context.get_search_node_id(d);

                const SearchNode& s = _context.get_search_node(_start);
                const SearchNodeId s_id = _context.get_search_node_id(s);

                path.push_back(d._node_it);

                SearchNodeId u_id = d_id;

                while ( u_id != s_id )
                {
                    const SearchNode& u_prev = _context.get_search_node_from_id(u_id);

                    const EdgeIterator e = u_prev._predecessor._edge_it;
                    assert ( e != INVALID_EDGE_ITERATOR );

                    u_id = u_prev._predecessor._search_node_id;
                    const SearchNode& u = _context.get_search_node_from_id(u_id);

                    assert( eq( _graph.get_ttf(e).eval(u._t_arr) + u._t_arr, u_prev._t_arr ) );
                    path.push_back(u._node_it);
                    pre_path.push_back(_graph.is_directed_downward(e));
                }
            }
            get_path(node);
            vector<unsigned> ranking;
            for(auto i : path) {
                ranking.push_back(_graph.get_ranking(i) );
            }
            for(auto i : path){
                if(_graph.get_ranking(i) > _graph.get_ranking(_destination) ){
                    unsigned aaa = _graph.get_ranking(i);
                    if(_bw_cpd.get_reachability(i,  _bw_mapped_destination) ==0XFF){
                        std::cout<<"well that is not true"<<std::endl;
                    }
                }
            }
        }

        class Path
        {
            friend BMTCH;

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

            Path(BMTCH_Path ts_path ){
                _nodes = ts_path._nodes;
                _edges = ts_path._edges;
                _t_arr = ts_path._t_arr;
            }
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

//            void check_up_down_path (ForwardSearchGraph<&period> graph, NodeIterator destination){
//                for(int i =0; i< _nodes.size(); i++){
//                   if(graph.get_ranking(_nodes[i]) > graph.get_ranking(destination)){
//                       if(_bw_cpd.get_first_move(_cpd_bw_mapper[_nodes[i]], _cpd_bw_mapper[destination]) ==0XFF){
//                           std::cout<<"well that is not true"<<std::endl;
//                       }
//                   }
//                }
//            }
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
            if(hourly_query_index != 9999){
                return Path(hourly_query[hourly_query_index]->get_path(node));
            }else{
                const auto up_down_path = get_up_down_path(node);
//            up_down_path.print_nodes(_graph);
                return  expand_path( up_down_path)  ;
            }
        }
    };
}
#endif //KATCH_BMTCH_H
