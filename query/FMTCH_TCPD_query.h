//
// Created by Bojie Shen on 19/10/21.
//

#ifndef KATCH_FMTCH_TCPD_QUERY_H
#define KATCH_FMTCH_TCPD_QUERY_H
#include <array>
#include <cassert>
#include <cmath>
#include <deque>
#include <limits>
#include <stack>
#include <tuple>
#include <vector>
#include <iostream>
#include <unordered_map>

#include <datastr/cpd/cpd.h>
#include <datastr/cpd/bw_cpd.h>
#include "util/id_queue.h"
#include "datastr/base/double.h"
#include "datastr/base/interval.h"
#include "datastr/graph/basic.h"
#include "datastr/graph/forward_search_graph.h"
#include "context/search_context.h"

#define INF 9999999;

namespace katch
{
    class FMTCH_TCPD_Path
    {

    public:

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

        FMTCH_TCPD_Path(const NodeIterator destination, const double t_arr)
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
                std::cout<<_nodes[i]<<"time:"<<_t_arr[i]<<std::endl;
            }
        }
    };



    class FMTCH_TCPDQuery_Wrapper {
    public:
        // pure virtual function providing interface framework.
//        virtual int getArea() = 0;
        unsigned long number_of_nodes_generated;
        unsigned long number_of_nodes_expanded;
        unsigned long number_of_first_move_calls;
        unsigned long number_of_TD_calculation;
        unsigned long number_of_bw_first_move_calls;

        virtual bool can_ts_tch_query_answer_the_query(const double& departure_time, const double& upper_duration) = 0;

        virtual double one_to_one(const NodeIterator& start, const NodeIterator& destination, const double& t_dep) = 0;

        virtual FMTCH_TCPD_Path get_path(const NodeIterator& node) const noexcept = 0;

    };



    template<const double* period>
    class FMTCH_TCPDQuery : public FMTCH_TCPDQuery_Wrapper
    {

    private:
        using SearchNodeId = uint32_t;
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

        ForwardSearchGraph<period>* _graph;
        std::vector<ForwardSearchGraph<period>> _graph_vector;
        NodeIterator _start;
        NodeIterator _destination;
        double _t_dep;
        double _t_arr;
        SearchContext<SearchNode_tmpl,SearchNodeId> _context;
        double _upper_bound;
        BW_CPD _bw_cpd;
        CPD _cpd;
        vector<unsigned> _cpd_mapper;
        vector<unsigned> _bw_cpd_mapper;


        std::unordered_map<int,unsigned> _cpd_start_index_mapper;
        unsigned _cpd_start_index;


        vector<double> _lower_bound_cache;
        vector<unsigned short> _was_cached;
        unsigned short _search_id;

        vector<unsigned short> _was_bw_cached;
        vector<bool> _was_bw_reachable;
        vector<unsigned short> _was_popped;
        unsigned _mapped_destination;
        unsigned _bw_mapped_destination;

//        double get_lower_bound(const unsigned start){
//            // recursion
//            unsigned mapped_source = _cpd_mapper[ start];
//            if(mapped_source == _mapped_destination){
//                return 0;
//            }
//            if(_was_cached[start] == _search_id){
//                return _lower_bound_cache[start];
//            }
//            unsigned char fm =  _cpd.get_first_move(mapped_source +_cpd_start_index, _mapped_destination);
//            number_of_first_move_calls++;
//            if(fm == 0XFF){
//                _was_cached[start] = _search_id;
//                _lower_bound_cache[start] = INF;
//                return INF;
//            }
//            unsigned next_id = _graph->get_other_node( (EdgeIterator)(_graph->edges_begin(start) + fm));
//            double h = get_lower_bound(next_id);
//            _lower_bound_cache[start] = h + _graph->get_ttf((EdgeIterator)(_graph->edges_begin(start) + fm)).get_min();
//            _was_cached[start] = _search_id;
//            return _lower_bound_cache[start];
//        }

//        bool is_down_reachable(const unsigned start){
//            if(_was_bw_cached[start] == _search_id){
//                return _was_bw_reachable[start];
//            }else{
//                _was_bw_cached[start] = _search_id;
//                unsigned char fm = _bw_cpd.get_first_move(_cpd_mapper[ start]+_cpd_start_index, _mapped_destination);
//                number_of_first_move_calls ++;
//                if(fm ==0XFF){
//                    _was_bw_reachable[start] = false;
//                    return false;
//                }
//                _was_bw_reachable[start] = true;
//                return true;
//            }
//
//        }


        double get_lower_bound(const unsigned& start){
            // recursion
            if(_was_cached[start] == _search_id){
                return _lower_bound_cache[start];
            }
            const unsigned char& fm = _cpd.get_first_move(start+_cpd_start_index, _mapped_destination);
            number_of_first_move_calls++;
            if(fm == 0XFF){
                _was_cached[start] = _search_id;
                _lower_bound_cache[start] = INF;
                return INF;
            }
            const unsigned& next_id = _graph->get_other_node( (EdgeIterator)(_graph->edges_begin(start) + fm));
            _lower_bound_cache[start] = get_lower_bound(next_id) + _graph->get_ttf((EdgeIterator)(_graph->edges_begin(start) + fm)).get_min();
            _was_cached[start] = _search_id;
            return _lower_bound_cache[start];
        }



        bool is_down_reachable(const unsigned start){
//            if(_graph.get_ranking(start) < _destination_ranking){
//                return false;
//            }
            if(_was_bw_cached[start] == _search_id){
                return _was_bw_reachable[start];
            }else{
                _was_bw_cached[start] = _search_id;
                _was_bw_reachable[start] = _bw_cpd.get_reachability(start+_cpd_start_index, _bw_mapped_destination);
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
                    for ( EdgeIterator e = _graph->downward_begin(u_it) ; e != _graph->downward_end(u_it) ; ++e ) {
                        const NodeIterator v_it = _graph->get_other_node(e);
                        if (!is_down_reachable(v_it)) continue;
                        if ( ! _context.reached(v_it) ) {
                            const double t_arr_v_new = _graph->get_ttf(e).eval(t_arr_u)+ t_arr_u;
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
                            if(t_arr_u + _graph->get_ttf(e).get_min() >= v._t_arr ){
                                continue;
                            }
                            const double t_arr_v_new = _graph->get_ttf(e).eval(t_arr_u)+ t_arr_u;
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
                    for ( EdgeIterator e = _graph->upward_begin(u_it) ; e != _graph->upward_end(u_it) ; ++e ) {
                        const NodeIterator v_it = _graph->get_other_node(e);
                        if ( ! _context.reached(v_it) ) {
                            const double t_arr_v_new = _graph->get_ttf(e).eval(t_arr_u)+ t_arr_u;
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
                            if(t_arr_u + _graph->get_ttf(e).get_min() >= v._t_arr ){
                                continue;
                            }
                            const double t_arr_v_new = _graph->get_ttf(e).eval(t_arr_u)+ t_arr_u;
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
        int hour_duration;
        using Graph = ForwardSearchGraph<period> ;

        FMTCH_TCPDQuery(std::vector<ForwardSearchGraph<period>>&& graph)
                : _graph_vector(std::move(graph)),
                  _graph(nullptr),
                  _start(INVALID_NODE_ITERATOR),
                  _destination(INVALID_NODE_ITERATOR),
                  _t_dep(std::numeric_limits<double>::max()),
                  _t_arr(std::numeric_limits<double>::max()),
                  _context(_graph_vector[0].get_n_nodes()),
                  _upper_bound(std::numeric_limits<double>::max()),
                  _was_cached(_graph_vector[0].get_n_nodes()),
                  _lower_bound_cache(_graph_vector[0].get_n_nodes()),
                  _was_popped(_graph_vector[0].get_n_nodes()),
                  _was_bw_cached(_graph_vector[0].get_n_nodes()),
                  _was_bw_reachable(_graph_vector[0].get_n_nodes()),
                  _search_id(0)
        {
            hour_duration = 24/ _graph_vector.size();

        }


        void load_cpd(const std::string& main_name, int shifting_hour){
            _cpd = CPD();
            _bw_cpd = BW_CPD();
            for(int hour = 0; hour < 24 ; hour += hour_duration){
                string cpd_name = main_name+ "_" + std::to_string(hour_duration) + "_" + std::to_string(shifting_hour) + "_" + std::to_string(hour)+ "_min.fw_ts_tch_cpd";
                string bw_cpd_name = main_name + "_"+ std::to_string(hour_duration) + "_" + std::to_string(shifting_hour) + "_" + std::to_string(hour)+ "_min.bw_ts_tch_cpd";

                CPD cur_fw_cpd = CPD();
                FILE*f = fopen( cpd_name.c_str(), "r");
                cur_fw_cpd.load(f);

                BW_CPD cur_bw_cpd = BW_CPD();
                FILE*bw_f = fopen( bw_cpd_name.c_str(), "r");
                cur_bw_cpd.load(bw_f);

                string bw_cpd_mapper_name = main_name + "_" +
                        std::to_string(hour_duration) + "_" + std::to_string(shifting_hour) + "_" + std::to_string(hour)+"_min.bw_ts_tch_mapper";
                vector<unsigned> bw_mapper = load_vector<unsigned>(bw_cpd_mapper_name );

                assert(cur_fw_cpd.node_count() ==_graph_vector[hour].get_n_nodes());
                assert(cur_bw_cpd.node_count() == _graph_vector[hour].get_n_nodes());

                _cpd.append_rows(cur_fw_cpd);
                _bw_cpd.append_rows(cur_bw_cpd);
                _bw_cpd_mapper.reserve( _bw_cpd_mapper.size() +  bw_mapper.size());
                _bw_cpd_mapper.insert(  _bw_cpd_mapper.end(), bw_mapper.begin(), bw_mapper.end() );
            }
            string mapper_name = main_name + "_" + std::to_string(hour_duration) + "_" + std::to_string(shifting_hour) + ".fw_ts_tch_mapper";
            _cpd_mapper = load_vector<unsigned>(mapper_name );
            for(int i = 0; i < 24 ; i ++){
                _cpd_start_index_mapper.insert({i,_graph_vector[i].get_n_nodes()*i});
            }
        }


        bool can_ts_tch_query_answer_the_query(const double& departure_time, const double& upper_duration) override{
            int expect_starting_hour  = (int) (departure_time/36000.0 / hour_duration ) * hour_duration ;
            double arrival_time  = departure_time + upper_duration;
            return expect_starting_hour * 36000.0 <= departure_time &&
                   arrival_time < expect_starting_hour * 36000.0 + *period - 5000;
        }

        double get_index_size(){
            double graph_size = 0;
            for(ForwardSearchGraph<period>& g: _graph_vector){
                graph_size += (double)g.get_index_size()/1000000;
            }
            std::cout<<"HTCH graph size: "<<graph_size<<" MB"<<std::endl;
            double fw_cpd_size  = (double)_cpd.get_entry().size()*4/1000000;
            std::cout<<"Forward CPD size: "<<fw_cpd_size<<" MB"<<std::endl;

            double bw_cpd_size  = (double)_bw_cpd.get_entry().size()*4/1000000;
            std::cout<<"Backward CPD size: "<<bw_cpd_size<<" MB"<<std::endl;

            double total_size = graph_size + fw_cpd_size + bw_cpd_size;
            return total_size  ;
        }



        double one_to_one(const NodeIterator& start, const NodeIterator& destination, const double& t_dep) override
        {
            number_of_first_move_calls = 0;
            number_of_nodes_generated = 0;
            number_of_nodes_expanded = 0;
            number_of_TD_calculation = 0;
            number_of_bw_first_move_calls = 0;

            if ( start == destination )
                return _t_dep;
            _t_arr = std::numeric_limits<double>::max();

            _start = start;
            _destination = destination;
            _t_dep = katch::util::modulo(t_dep,36000.0)/hour_duration;
            double _t_base = t_dep - _t_dep;

            _graph = &_graph_vector[katch::util::modulo((int)(t_dep/36000.0),24)/hour_duration];
            int cpd_index = util::modulo(t_dep,864000.0)/ 36000;
            _cpd_start_index = _cpd_start_index_mapper[cpd_index];

            _mapped_destination = _cpd_mapper[_destination];
            _bw_mapped_destination = _bw_cpd_mapper[_cpd_start_index + _destination];

            _upper_bound = std::numeric_limits<double>::max();
            _context.clear_all();

            _t_arr = search();
            return _t_arr+_t_base;
        }

        FMTCH_TCPD_Path get_up_down_path(NodeIterator node) const noexcept
        {
            if ( ! _context.reached(_start) ) return FMTCH_TCPD_Path(_start, _t_dep);
            if ( ! _context.reached(node) ) return FMTCH_TCPD_Path(_start, _t_dep);

            const SearchNode& d = _context.get_search_node(node);
            const SearchNodeId d_id = _context.get_search_node_id(d);

            const SearchNode& s = _context.get_search_node(_start);
            const SearchNodeId s_id = _context.get_search_node_id(s);

            FMTCH_TCPD_Path result(d._node_it, d._t_arr);

            SearchNodeId u_id = d_id;

            while ( u_id != s_id )
            {
                const SearchNode& u_prev = _context.get_search_node_from_id(u_id);

                const EdgeIterator e = u_prev._predecessor._edge_it;
                assert ( e != INVALID_EDGE_ITERATOR );

                u_id = u_prev._predecessor._search_node_id;
                const SearchNode& u = _context.get_search_node_from_id(u_id);

                assert( eq( _graph->get_ttf(e).eval(u._t_arr) + u._t_arr, u_prev._t_arr ) );
                result.append_front(u._node_it, u._t_arr, e);
            }

            assert( result._nodes.front() == _start );
            assert( result._nodes.back() == _destination );
            assert( result._t_arr.front() == _t_dep );
            assert( fabs(result._t_arr.back() - _t_arr) <= 0.00001 );
//            assert( result.dbg_check_path(*_graph, _t_dep) );

            return result;
        }

        FMTCH_TCPD_Path expand_path(const FMTCH_TCPD_Path& up_down_path) const noexcept
        {
            FMTCH_TCPD_Path result(_start, _t_dep);
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

                    NodeIterator x = _graph->get_middle_node(e, t_current);

                    if ( x == INVALID_NODE_ITERATOR )
                    {
                        t_current += _graph->get_ttf(e).eval(t_current);
                        result.append_back(e, v, t_current);
                        continue;
                    }

                    EdgeIterator e_up = INVALID_EDGE_ITERATOR;
                    EdgeIterator e_down = INVALID_EDGE_ITERATOR;

                    for ( EdgeIterator it = _graph->edges_begin(x) ; it != _graph->edges_end(x) ; ++it )
                    {
                        if ( _graph->get_other_node(it) == v && _graph->is_directed_upward(it) ) e_up = it;
                        if ( e_up != INVALID_EDGE_ITERATOR ) break;
                    }
                    for ( EdgeIterator it = _graph->edges_begin(u) ; it != _graph->edges_end(u) ; ++it )
                    {
                        if ( _graph->get_other_node(it) == x && _graph->is_directed_downward(it) ) e_down = it;
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
//            assert( result.dbg_check_path(*_graph, _t_dep) );
//            assert( result.dbg_check_path_shortcut_free(*_graph, _t_dep) );

            return result;
        }


        FMTCH_TCPD_Path get_path(const NodeIterator& node) const noexcept override{
            auto up_down_path = get_up_down_path(node);
            return  expand_path( up_down_path) ;
        }
    };
}
#endif //KATCH_FMTCH_TCPD_QUERY_H
