//
// Created by Bojie Shen on 11/5/21.
//

#ifndef KATCH_BTCH_H
#define KATCH_BTCH_H
#include <array>
#include <cassert>
#include <cmath>
#include <deque>
#include <limits>
#include <stack>
#include <tuple>
#include <vector>

//#include <boost/heap/pairing_heap.hpp>
#include "util/id_queue.h"
#include "datastr/base/double.h"
#include "datastr/base/interval.h"
#include "datastr/graph/basic.h"
#include "datastr/graph/search_graph.h"
#include "context/bi_search_context.h"
namespace katch
{

    class BTCH
    {

    private:
        using SearchNodeId = uint32_t;
        static constexpr double period = 864000.0;
        struct Predecessor {
            EdgeIterator _edge_it;
            SearchNodeId _search_node_id;

            Predecessor()
                    : _edge_it(INVALID_EDGE_ITERATOR), _search_node_id(std::numeric_limits<SearchNodeId>::max()) {}

            Predecessor(const EdgeIterator &e, const SearchNodeId id)
                    : _edge_it(e), _search_node_id(id) {}
        };

        struct SearchNode_tmpl {
            NodeIterator _node_it;
            std::bitset<2> _enqueued;
            double _t_arr;
            double _t_stall;
            Interval _interval;
            double _upper_stall;
            Predecessor _fw_predecessor;
            std::vector <Predecessor> _bw_predecessors;

            SearchNode_tmpl(const SearchNode_tmpl &) = delete;

            SearchNode_tmpl &operator=(const SearchNode_tmpl &) = delete;

            SearchNode_tmpl(SearchNode_tmpl &&) = default;

            SearchNode_tmpl &operator=(SearchNode_tmpl &&) = default;

            SearchNode_tmpl()
                    :
                    _node_it(INVALID_NODE_ITERATOR),
                    _enqueued(uint32_t(0)),
                    _t_arr(std::numeric_limits<double>::max()),
                    _t_stall(std::numeric_limits<double>::max()),
                    _interval(),
                    _upper_stall(std::numeric_limits<double>::max()),
                    _fw_predecessor(),
                    _bw_predecessors() {}
        };
        static constexpr size_t FORWARD = BiSearchContext<SearchNode_tmpl,SearchNodeId>::FORWARD;
        static constexpr size_t BACKWARD = BiSearchContext<SearchNode_tmpl,SearchNodeId>::BACKWARD;

    public:

        using SearchNode = SearchNode_tmpl;

    private:

        static size_t flip(const size_t direction)
        {
            if ( direction == FORWARD )
                return BACKWARD;

            return FORWARD;
        }

        double downward_search()
        {
            if ( _upper_bound == std::numeric_limits<double>::max() )
                return std::numeric_limits<double>::max();

            _context.clear_pq(FORWARD);

            for ( const auto& search_node_id : _candidates )
            {
                SearchNode& u = _context.get_search_node_from_id(search_node_id);

                if ( u._t_arr + u._interval.get_lower() <= _upper_bound + 0.01 )
                    if ( ! _context.pq_contains(u, FORWARD) ){
                        _context.pq_re_insert(u, u._t_arr, FORWARD);
                        number_of_nodes_generated++;
                    }

            }


            while ( ! _context.pq_empty(FORWARD) )
            {
                const SearchNode& u = _context.pq_delete_min(FORWARD);
                number_of_nodes_expanded++;
                if ( u._node_it == _destination ) break;

                for ( auto& pred : u._bw_predecessors )
                {
                    const EdgeIterator e = pred._edge_it;
                    SearchNode& v = _context.get_search_node_from_id(pred._search_node_id);
                    assert( _graph.is_directed_downward(e) );
                    assert( _graph.get_other_node(e) == u._node_it );

                    const double t_arr_u = u._t_arr;

                    if ( v._t_arr == std::numeric_limits<double>::max() )
                    {
                        const double t_arr_v_new = _graph.get_ttf(e).eval(t_arr_u) + t_arr_u;
                        number_of_TD_calculation ++;
                        _context.pq_re_insert(v, t_arr_v_new, FORWARD);
                        number_of_nodes_generated++;
                        v._t_arr = t_arr_v_new;
                        v._fw_predecessor = Predecessor(e, _context.get_search_node_id(u));
                    }
                    else
                    {
                        if(t_arr_u + _graph.get_ttf(e).get_min() >= v._t_arr ){
                            continue;
                        }
                        const double t_arr_v_new = _graph.get_ttf(e).eval(t_arr_u)+ t_arr_u;
                        number_of_TD_calculation ++;

                        if ( t_arr_v_new >= v._t_arr ) continue;

                        v._t_arr = t_arr_v_new;
                        v._fw_predecessor = Predecessor(e, _context.get_search_node_id(u));

                        if ( _context.pq_contains(v, FORWARD) )
                            _context.pq_decrease(v, t_arr_v_new, FORWARD);
                        else
                            _context.pq_re_insert(v, t_arr_v_new, FORWARD);
                        number_of_nodes_generated++;
                    }
                }
            }

            assert( _context.reached(_destination) );
            const SearchNode& d = _context.get_search_node(_destination);

            return d._t_arr;
        }

        void bidirectional_search()
        {
            assert( _start != _destination );

            size_t direction = BACKWARD;
            std::array<bool,2> finished = { false, false };

            SearchNode& s = _context.insert(_start, _t_dep, FORWARD);
            number_of_nodes_generated++;
            s._t_arr = _t_dep;

            SearchNode& d = _context.insert(_destination, 0.0, BACKWARD);
            number_of_nodes_generated++;
            d._interval = Interval(0.0, 0.0);

            while ( ! _context.pq_empty(FORWARD) || ! _context.pq_empty(BACKWARD) )
            {
                if ( finished[FORWARD] && finished[BACKWARD] ) break;
                if ( ! _context.pq_empty(flip(direction)) ) direction = flip(direction);
                if ( finished[FORWARD] ) direction = BACKWARD;
                if ( finished[BACKWARD] ) direction = FORWARD;

                if ( _context.pq_empty(direction) )
                {
                    finished[direction] = true;
                    continue;
                }

                if ( direction == FORWARD )
                {
                    SearchNode& u = _context.pq_delete_min(FORWARD);
                    number_of_nodes_expanded++;
                    if ( _upper_bound != std::numeric_limits<double>::max() && u._t_arr > _upper_bound + 0.001 )
                    {
                        finished[FORWARD] = true;
                        continue;
                    }

                    bool stalled = false;
                    for ( EdgeIterator e = _graph.edges_begin(u._node_it) ; e != _graph.edges_end(u._node_it) ; ++e )
                    {
                        //stall on demand pruning here ;
                        if ( ! _graph.is_directed_downward(e) ) continue;

                        const NodeIterator w_it = _graph.get_other_node(e);
                        if ( ! _context.reached(w_it) ) continue;

                        SearchNode& w = _context.get_search_node(w_it);
                        if ( w._t_arr == std::numeric_limits<double>::max() ) continue;

                        const double t_stall_w = std::min(w._t_arr, w._t_stall);
                        if(u._t_stall > (_graph.get_ttf(e).get_min() + t_stall_w) ){
                            u._t_stall = std::min(u._t_stall, _graph.get_ttf(e).eval(t_stall_w) + t_stall_w);
                            number_of_TD_calculation ++;
                        }
                        if ( u._t_stall + 0.001 < u._t_arr )
                        {
                            stalled = true;
                            break;
                        }
                    }

                    if ( stalled ) continue;

                    if ( u._interval != Interval::INFTY )
                    {
                        if (
                                _upper_bound == std::numeric_limits<double>::max() ||
                                u._t_arr + u._interval.get_lower() <= _upper_bound + 0.001
                                )
                            _candidates.push_back(_context.get_search_node_id(u));
                        // update upper bound
                        _upper_bound = std::min( _upper_bound, u._t_arr + u._interval.get_upper() );
                    }

                    const SearchNodeId u_id = _context.get_search_node_id(u);
                    const NodeIterator u_it = u._node_it;
                    const double t_arr_u = u._t_arr;

                    for ( EdgeIterator e = _graph.edges_begin(u_it) ; e != _graph.edges_end(u_it) ; ++e )
                    {
                        if ( ! _graph.is_directed_upward(e) ) continue;

                        const NodeIterator v_it = _graph.get_other_node(e);
                        if ( ! _context.reached(v_it) )
                        {
                            const double t_arr_v_new = _graph.get_ttf(e).eval(t_arr_u)+ t_arr_u;
                            number_of_TD_calculation ++;
                            SearchNode& v = _context.insert(v_it, t_arr_v_new, FORWARD);
                            v._t_arr = t_arr_v_new;
                            v._fw_predecessor = Predecessor(e, u_id);
                            number_of_nodes_generated++;
                        }
                        else
                        {
                            assert( _context.reached(v_it) );
                            SearchNode& v = _context.get_search_node(v_it);

                            if ( v._t_arr != std::numeric_limits<double>::max() ){
                                if(t_arr_u + _graph.get_ttf(e).get_min() >= v._t_arr ){
                                    continue;
                                }
                            }
                            const double t_arr_v_new = _graph.get_ttf(e).eval(t_arr_u)+ t_arr_u;
                            number_of_TD_calculation ++;

                            if ( v._t_arr != std::numeric_limits<double>::max() )
                                if ( t_arr_v_new >= v._t_arr )
                                    continue;
                            v._t_arr = t_arr_v_new;
                            v._fw_predecessor = Predecessor(e, u_id);
                            // search still follows FIFO, but just they maintain search node in one place for both side
                            // therefore, some nodes are reached but not in the queue
                            if ( _context.pq_contains(v, FORWARD) )
                                _context.pq_decrease(v, t_arr_v_new, FORWARD);
                            else
                                _context.pq_re_insert(v, t_arr_v_new, FORWARD);
                            number_of_nodes_generated++;
                        }
                    }
                }
                else
                {
                    SearchNode& u = _context.pq_delete_min(BACKWARD);
                    number_of_nodes_expanded++;
                    if ( _upper_bound != std::numeric_limits<double>::max() && u._interval.get_lower() >  _upper_bound - _t_dep + 0.001 )
                    {
                        finished[BACKWARD] = true;
                        continue;
                    }

                    bool stalled = false;
                    for ( EdgeIterator e = _graph.edges_begin(u._node_it) ; e != _graph.edges_end(u._node_it) ; ++e )
                    {
                        if ( ! _graph.is_directed_upward(e) ) continue;

                        const NodeIterator w_it = _graph.get_other_node(e);
                        if ( ! _context.reached(w_it) ) continue;

                        SearchNode& w = _context.get_search_node(w_it);
                        if ( w._interval == Interval::INFTY ) continue;

                        u._upper_stall = std::min( u._upper_stall, w._interval.get_upper() + _graph.get_ttf(e).get_max() );
                        u._upper_stall = std::min( u._upper_stall, w._upper_stall + _graph.get_ttf(e).get_max() );

                        if ( u._upper_stall + 0.001 < u._interval.get_lower() )
                        {
                            stalled = true;
                            break;
                        }
                    }

                    if ( stalled ) continue;

                    if ( u._t_arr != std::numeric_limits<double>::max() )
                    {
                        if (
                                _upper_bound == std::numeric_limits<double>::max() ||
                                u._t_arr + u._interval.get_lower() <= _upper_bound + 0.001
                                )
                            _candidates.push_back(_context.get_search_node_id(u));

                        _upper_bound = std::min( _upper_bound, u._t_arr + u._interval.get_upper() );
                    }

                    const SearchNodeId u_id = _context.get_search_node_id(u);
                    const NodeIterator u_it = u._node_it;
                    const Interval interval_u = u._interval;

                    for ( EdgeIterator e = _graph.edges_begin(u_it) ; e != _graph.edges_end(u_it) ; ++e )
                    {
                        if ( ! _graph.is_directed_downward(e) ) continue;

                        const NodeIterator v_it = _graph.get_other_node(e);
                        const double min_e = _graph.get_ttf(e).get_min();
                        const double max_e = _graph.get_ttf(e).get_max();

                        const Interval interval_v_new(interval_u.get_lower() + min_e, interval_u.get_upper() + max_e);

                        if ( ! _context.reached(v_it) )
                        {
                            SearchNode& v = _context.insert(v_it, interval_v_new.get_lower(), BACKWARD);
                            v._interval = interval_v_new;
                            v._bw_predecessors.emplace_back(e, u_id);
                            number_of_nodes_generated++;
                        }
                        else
                        {
                            assert( _context.reached(v_it) );
                            SearchNode& v = _context.get_search_node(v_it);

                            if ( v._interval != Interval::INFTY )
                            {
                                if ( interval_v_new.get_lower() > v._interval.get_upper() + 0.001 )
                                    continue;

                                if ( interval_v_new.get_upper() + 0.001 < v._interval.get_lower() )
                                    v._bw_predecessors.clear();
                            }

                            v._bw_predecessors.emplace_back(e, u_id);

                            if ( v._interval != Interval::INFTY )
                                if ( interval_v_new.get_lower() >= v._interval.get_lower() )
                                    if( interval_v_new.get_upper() >= v._interval.get_upper() )
                                        continue;

                            v._interval = merge(v._interval, interval_v_new);

                            if ( _context.pq_contains(v, BACKWARD) )
                                _context.pq_decrease(v, v._interval.get_lower(), BACKWARD);
                            else
                                _context.pq_re_insert(v, v._interval.get_lower(), BACKWARD);
                            number_of_nodes_generated++;
                        }
                    }
                }
            }
        }

        SearchGraph<&period> _graph;

        NodeIterator _start;
        NodeIterator _destination;
        double _t_dep;
        double _t_arr;
        BiSearchContext<SearchNode_tmpl,SearchNodeId> _context;
        std::vector<SearchNodeId> _candidates;
        double _upper_bound;

    public:

        unsigned long number_of_nodes_generated;
        unsigned long number_of_nodes_expanded;
        unsigned long number_of_first_move_calls;
        unsigned long number_of_TD_calculation;
        unsigned long number_of_bw_first_move_calls;

        using Graph = SearchGraph<&period>;

        BTCH(SearchGraph<&period> && graph)
                : _graph(std::move(graph)),
                  _start(INVALID_NODE_ITERATOR),
                  _destination(INVALID_NODE_ITERATOR),
                  _t_dep(std::numeric_limits<double>::max()),
                  _t_arr(std::numeric_limits<double>::max()),
                  _context(_graph.get_n_nodes()),
                  _upper_bound(std::numeric_limits<double>::max())
        {}

        double get_index_size(){
            double graph_size = _graph.get_index_size()/1000000;
            std::cout<<"TCH graph size: "<<graph_size<<" MB"<<std::endl;
            return graph_size;
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

            if ( start == destination )
                return _t_dep;

            _upper_bound = std::numeric_limits<double>::max();
            _candidates.clear();
            _context.clear_all();

            bidirectional_search();
            _t_arr = downward_search();

            return _t_arr;
        }

        class Path
        {
            friend BTCH;

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

            bool dbg_check_path(const SearchGraph<&period> & graph, const double t_dep)
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


            bool dbg_check_path_shortcut_free(const SearchGraph<&period> & graph, const double t_dep)
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
                    std::cout<<_nodes[i]<<"time:"<<_t_arr[i]<<std::endl;
                }
            }
            void print_nodes(const SearchGraph<&period> & graph) const {

                for(int i =0; i< _nodes.size(); i++){
                    std::cout<<_nodes[i]<<" time: "<<_t_arr[i]<<std::endl;
                }
                double cost = 0;
                for(int i = 0 ; i < _edges.size(); i++){

                    std::cout<<" edge id: "<< _edges[i]<< "Direction" << graph.is_directed_upward(_edges[i])<<graph.is_directed_downward(_edges[i]) << std::endl;
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

                const EdgeIterator e = u_prev._fw_predecessor._edge_it;
                assert ( e != INVALID_EDGE_ITERATOR );

                u_id = u_prev._fw_predecessor._search_node_id;
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
                        if ( _graph.get_other_node(it) == u && _graph.is_directed_downward(it) ) e_down = it;
                        if ( _graph.get_other_node(it) == v && _graph.is_directed_upward(it) ) e_up = it;

                        if ( e_down != INVALID_EDGE_ITERATOR && e_up != INVALID_EDGE_ITERATOR ) break;
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
            const katch::BTCH::Path up_down_path = get_up_down_path(node);
//            up_down_path.print_nodes(_graph);
            return  expand_path( up_down_path) ;
//            up_down_path.print_nodes(_graph);
//            return get_up_down_path(node);
        }
    };
}
#endif //KATCH_BTCH_H
