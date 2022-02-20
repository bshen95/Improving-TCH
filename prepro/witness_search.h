/*
 * katch/datastr/base/witness_search.h
 *
 *
 * This file is part of
 *
 * KaTCH -- Karlsruhe Time-Dependent Contraction Hierarchies
 *
 * Copyright (C) 2015
 *
 * Institut fuer Theroretische Informatik,
 * Karlsruher Institut fuer Technology (KIT),
 * 76131 Karlsruhe, Germany
 *
 * Author: Gernot Veit Batz
 *
 *
 *
 * KaTCH is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Affero General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaTCH is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with Contraction Hierarchies; see the file COPYING;
 * if not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef KATCH_WITNESS_SEARCH_H_
#define KATCH_WITNESS_SEARCH_H_

#include <algorithm>
#include <cstdint>
#include <iterator>
#include <limits>
#include <unordered_map>
#include <vector>

//#include <boost/heap/pairing_heap.hpp>

#include "datastr/graph/dynamic_graph.h"
#include "datastr/graph/basic.h"
#include "datastr/graph/order_edge_data.h"
#include "prepro/edge_description.h"
#include "util/id_queue.h"

namespace katch
{
    template<const double* period>
    class WitnessSearch
    {

    public:

        static constexpr unsigned int UNDECIDED = 0;
        static constexpr unsigned int NECESSARY = 1;
        static constexpr unsigned int NOT_NECESSARY = 2;

        using ShortcutStatus = unsigned int;

    private:

        using Graph = DynamicGraph<OrderEdgeData<period>>;
        using EdgeData = typename Graph::EdgeData;

        using SearchNodeId = uint32_t;
        using TTF = typename OrderEdgeData<period>::TTF;
        using UndercutDescriptor = typename OrderEdgeData<period>::TTF::UndercutDescriptor;

        struct Predecessor
        {
            BwEdgeIterator _bw_edge_it;
            SearchNodeId _search_node_id;

            Predecessor()
                    : _bw_edge_it(INVALID_BW_EDGE_ITERATOR),
                      _search_node_id(std::numeric_limits<SearchNodeId>::max())
            {}

            Predecessor(const BwEdgeIterator& e, const SearchNodeId& id)
                    : _bw_edge_it(e), _search_node_id(id)
            {}
        };

        struct SearchNode_tmpl
        {
            NodeIterator _node_it;
            bool _enqueued : 8;
            uint8_t _n_hops_sample;
            uint8_t _n_hops_interval;
            uint8_t _n_hops_profile;
            double _t_arr;
            Interval _interval;
            TTF _ttf;
            std::pair<Predecessor,Predecessor> _predecessors;

            SearchNode_tmpl(const SearchNode_tmpl&) = delete;
            SearchNode_tmpl& operator=(const SearchNode_tmpl&) = delete;

            SearchNode_tmpl(SearchNode_tmpl&&) = default;
            SearchNode_tmpl& operator=(SearchNode_tmpl&&) = default;

            SearchNode_tmpl()
                    :
                    _node_it(INVALID_NODE_ITERATOR), _enqueued(false),
                    _n_hops_sample(0), _n_hops_interval(0), _n_hops_profile(0),
                    _t_arr(std::numeric_limits<double>::max()), _interval(), _ttf(), _predecessors()
            {}
        };

        class SearchContext
        {

        private:

            using HashTable = std::unordered_map<NodeIterator, SearchNodeId>;

            MinIDQueue _heap;
            HashTable _hash_table;
            std::vector<SearchNode_tmpl> _search_nodes;

        public:

            using SearchNode = SearchNode_tmpl;

            SearchContext(const size_t n_nodes)
                    : _heap(n_nodes), _hash_table(), _search_nodes()
            { }

            bool pq_empty() const
            {
                return _heap.empty();
            }

            bool pq_contains(const NodeIterator& node_iterator) const
            {
                assert( _hash_table.count(node_iterator) < 2 );

                const auto it = _hash_table.find(node_iterator);
                if ( it == _hash_table.end() ) return false;

                return _search_nodes[it->second]._enqueued;
            }

            bool pq_contains(const SearchNode& search_node) const
            {
                return search_node._enqueued;
            }

            bool reached(const NodeIterator& node_iterator) const
            {
                return _hash_table.find(node_iterator) != _hash_table.end();
            }

            const SearchNode& get_search_node(const NodeIterator& node_iterator) const
            {
                assert( _hash_table.find(node_iterator) != _hash_table.end() );
                return _search_nodes[_hash_table.find(node_iterator)->second];
            }

            SearchNode& get_search_node(const NodeIterator& node_iterator)
            {
                assert( _hash_table.find(node_iterator) != _hash_table.end() );
                return _search_nodes[_hash_table.find(node_iterator)->second];
            }

            SearchNode& get_search_node_from_id(const SearchNodeId& search_node_id)
            {
                assert( search_node_id < _search_nodes.size() );
                return _search_nodes[search_node_id];
            }

            SearchNodeId get_search_node_id(const SearchNode& search_node) const
            {
                assert( reached(search_node._node_it) );

                const SearchNodeId search_node_id = std::distance(&(*(_search_nodes.begin())), &search_node);
                return search_node_id;
            }

            const SearchNode& get_min() const
            {
                assert( ! pq_empty() );
                return _search_nodes[_heap.peek().id];
            }

            const double get_min_priority() const
            {
                assert( ! pq_empty() );
                return _heap.peek().key;
            }

            SearchNode& insert(const NodeIterator& node_iterator, const double& priority)
            {
                assert( ! reached(node_iterator) );
                assert( ! pq_contains(node_iterator) );

                const SearchNodeId search_node_id(_search_nodes.size());

                _search_nodes.emplace_back();
                _heap.push({search_node_id,priority});

                _hash_table[node_iterator] = search_node_id;
                _search_nodes.back()._node_it = node_iterator;
                _search_nodes.back()._enqueued = true;

                return _search_nodes.back();
            }

            void pq_re_insert(SearchNode& search_node, const double& priority)
            {
                assert( reached(search_node._node_it) );
                assert( ! pq_contains(search_node._node_it) );

                const SearchNodeId search_node_id = get_search_node_id(search_node);
                _heap.push({search_node_id,priority});
                search_node._enqueued = true;

            }

            void pq_decrease(SearchNode& search_node, const double& priority)
            {
                assert( reached(search_node._node_it) );
                assert( pq_contains(search_node._node_it) );

                const SearchNodeId search_node_id = get_search_node_id(search_node);
                _heap.decrease_key({search_node_id,priority});
            }

            SearchNode& pq_delete_min()
            {
                assert( ! pq_empty() );

                SearchNode& search_node = _search_nodes[_heap.peek().id];
                assert( _hash_table.find(search_node._node_it) != _hash_table.end() );
                _heap.pop();

                search_node._enqueued = false;
                return search_node;
            }

            void clear_pq()
            {
                for(auto& search_node:_search_nodes )
                    search_node._enqueued =false;
                _heap.clear();
            }

            void clear_all()
            {
                _heap.clear();
                _hash_table.clear();
                _search_nodes.clear();
            }
        };

        using SearchNode = typename SearchContext::SearchNode;

        const Graph* _graph;
        SearchContext _context;

        void backward_interval_search(const NodeIterator destination_it, const NodeIterator start_it, const size_t& hop_limit = 16)
        {
            SearchNode& d = _context.insert(destination_it, 0.0);
            d._n_hops_interval = 0;
            d._interval = Interval(0.0, 0.0);

            while ( ! _context.pq_empty() )
            {
                if ( _context.reached(start_it) )
                {
                    const SearchNode& s = _context.get_search_node(start_it);

                    if ( s._interval != Interval::INFTY )
                        if ( ge( _context.get_min_priority(), s._interval.get_upper() ) )
                            return;
                }

                SearchNode& u = _context.pq_delete_min();

                if ( u._n_hops_interval >= hop_limit ) continue;

                const NodeIterator u_it = u._node_it;
                const double u_lower = u._interval.get_lower();
                const double u_upper = u._interval.get_upper();
                const uint32_t u_hops = u._n_hops_interval;
                const SearchNodeId u_id = _context.get_search_node_id(u);

                for ( auto e = _graph->in_edges_begin(u_it) ; e != _graph->in_edges_end(u_it) ; ++e )
                {
                    const NodeIterator v_it = _graph->get_source(e);
                    const EdgeData& edge_vu = _graph->get_edge_data(e);

                    const double min_e = edge_vu.get_ttf().get_min();
                    const double max_e = edge_vu.get_ttf().get_max();

                    const Interval interval_v_new(u_lower + min_e, u_upper + max_e);

                    if ( ! _context.reached(v_it) )
                    {
                        SearchNode& v = _context.insert(v_it, interval_v_new.get_lower());

                        v._n_hops_interval = u_hops + 1;
                        v._interval = interval_v_new;
                        v._predecessors.first = Predecessor(e, u_id);
                        v._predecessors.second = Predecessor(e, u_id);
                    }
                    else
                    {
                        assert( _context.reached(v_it) );

                        SearchNode& v = _context.get_search_node(v_it);
                        assert( v._interval != Interval::INFTY );

                        if (
                                ge( interval_v_new.get_lower(), v._interval.get_lower() ) &&
                                ge( interval_v_new.get_upper(), v._interval.get_upper() )
                                ) continue;

                        if ( lt( interval_v_new.get_upper(), v._interval.get_upper() ) )
                            v._predecessors.first = Predecessor(e, u_id);

                        if ( lt( interval_v_new.get_lower(), v._interval.get_lower() ) )
                            v._predecessors.second = Predecessor(e, u_id);

                        v._n_hops_interval = std::max(uint8_t(u_hops + 1), v._n_hops_interval);

                        const Interval interval_merged = merge(v._interval, interval_v_new);
                        v._interval = interval_merged;

                        if ( ! _context.pq_contains(v) )
                            _context.pq_re_insert(v, interval_merged.get_lower());
                        else
                            _context.pq_decrease(v, interval_merged.get_lower());
                    }
                }
            }

            return;
        }

        void sample_search
                (
                        const NodeIterator& start_it,
                        const NodeIterator& destination_it,
                        const double& t_start,
                        const size_t& hop_limit = 16
                )
        {
            assert( _context.reached(start_it) );

            SearchNode& s = _context.get_search_node(start_it);
            s._n_hops_sample = 0;
            s._t_arr = t_start;

            _context.pq_re_insert(s, t_start + s._interval.get_lower());

            while ( ! _context.pq_empty() )
            {
                assert( _context.reached(destination_it) );
                const SearchNode& d = _context.get_search_node(destination_it);

                if ( d._t_arr != std::numeric_limits<double>::max() )
                    if ( ge(_context.get_min_priority(), d._t_arr) )
                        return;

                SearchNode& u = _context.pq_delete_min();

                if ( u._n_hops_sample >= hop_limit ) continue;

                if ( u._predecessors.first._bw_edge_it == u._predecessors.second._bw_edge_it )
                    u._predecessors.second._bw_edge_it = INVALID_BW_EDGE_ITERATOR;

                for ( auto& pred : { u._predecessors.first, u._predecessors.second } )
                {
                    const BwEdgeIterator e = pred._bw_edge_it;
                    if ( e == INVALID_BW_EDGE_ITERATOR ) continue;

                    SearchNode& v = _context.get_search_node_from_id(pred._search_node_id);
                    assert ( _context.reached(v._node_it) );

                    const EdgeData& edge_uv = _graph->get_edge_data(e);
                    const double t_arr_u = u._t_arr;
                    const double t_arr_v_new = edge_uv.get_ttf().eval(t_arr_u) + t_arr_u;


                    if ( v._t_arr == std::numeric_limits<double>::max() )
                    {
                        _context.pq_re_insert(v, t_arr_v_new + v._interval.get_lower());

                        v._n_hops_sample = u._n_hops_sample + 1;
                        v._t_arr = t_arr_v_new;
                    }
                    else
                    {
                        if ( ge(t_arr_v_new, v._t_arr) ) continue;

                        v._n_hops_sample = u._n_hops_sample + 1;
                        v._t_arr = t_arr_v_new;

                        if ( _context.pq_contains(v) )
                            _context.pq_decrease(v, t_arr_v_new + v._interval.get_lower());
                        else
                            _context.pq_re_insert(v, t_arr_v_new + v._interval.get_lower());
                    }
                }
            }

            return;
        }

        void profile_search(const NodeIterator& start_it, const NodeIterator& destination_it, const size_t& hop_limit = 16)
        {
            SearchNode& s = _context.get_search_node(start_it);
            s._ttf = TTF(0.0);
            s._n_hops_profile = 0;

            _context.pq_re_insert(s, s._interval.get_lower());

            const SearchNode& d = _context.get_search_node(destination_it);

            while ( ! _context.pq_empty() )
            {
                if ( d._ttf != TTF::INFTY )
                    if ( gt( _context.get_min_priority(), d._ttf.get_max() ) )
                        return;

                SearchNode& u = _context.pq_delete_min();

                if ( u._n_hops_profile >= hop_limit ) continue;

                if ( u._predecessors.first._bw_edge_it == u._predecessors.second._bw_edge_it )
                    u._predecessors.second._bw_edge_it = INVALID_BW_EDGE_ITERATOR;

                for ( auto pred : { u._predecessors.first, u._predecessors.second } )
                {
                    BwEdgeIterator e = pred._bw_edge_it;
                    if ( e == INVALID_BW_EDGE_ITERATOR ) continue;

                    SearchNode& v = _context.get_search_node_from_id(pred._search_node_id);
                    assert( v._interval != Interval::INFTY );

                    const EdgeData& edge_uv = _graph->get_edge_data(e);

                    if ( v._ttf != TTF::INFTY )
                    {
                        if ( ge(u._ttf.get_min() + edge_uv.get_ttf().get_min() , v._ttf.get_max()) )
                            continue;

                        if (
                                gt(
                                        u._ttf.get_min() + edge_uv.get_ttf().get_min() + v._interval.get_lower(),
                                        v._ttf.get_max() + v._interval.get_upper()
                                )
                                )
                            continue;
                    }

                    TTF f_v_new = link(edge_uv.get_ttf(), u._ttf);

                    if ( v._ttf == TTF::INFTY )
                    {
                        _context.pq_re_insert(v, f_v_new.get_min() + v._interval.get_lower());
                        v._ttf = std::move(f_v_new);
                        v._n_hops_profile = u._n_hops_profile + 1;
                    }
                    else
                    {
                        if ( ge(f_v_new.get_min(), v._ttf.get_max()) ) continue;

                        UndercutDescriptor descr;
                        merge(f_v_new.add(0.00001), v._ttf, descr);
                        if ( ! descr._f_undercuts_strictly ) continue;

                        TTF f_v_merge = merge(f_v_new, v._ttf);
                        double priority_v = f_v_merge.get_min();

                        v._n_hops_profile = std::max(v._n_hops_profile, uint8_t(u._n_hops_profile + 1));
                        v._ttf = std::move(f_v_merge);

                        if ( ! _context.pq_contains(v) )
                            _context.pq_re_insert(v, priority_v + v._interval.get_lower());
                        else
                            _context.pq_decrease(v, priority_v + v._interval.get_lower());
                    }
                }
            }

            return;
        }

    public:

        WitnessSearch(const Graph* graph)
                : _graph(graph), _context(_graph->get_n_nodes())
        {}

        ShortcutStatus run
                (
                        const NodeIterator& u_it,
                        const EdgeData& edge_ux,
                        const NodeIterator& x_it,
                        const EdgeData& edge_xv,
                        const NodeIterator& v_it,
                        EdgeDescription<period>& insert_edge
                )
        {
            _context.clear_all();

            //
            // backward interval search
            //
            backward_interval_search(v_it, u_it);
            if ( ! _context.reached(u_it) ) return NECESSARY;

            const SearchNode& u = _context.get_search_node(u_it);
            const Interval& interval_uv = u._interval;

            if ( interval_uv == Interval::INFTY ) return NECESSARY;

            if ( interval_uv.get_upper() + 0.01 < edge_ux.get_ttf().get_min() + edge_xv.get_ttf().get_min() )
                return NOT_NECESSARY;

            insert_edge._ttf = link(edge_xv.get_ttf(), edge_ux.get_ttf());
            insert_edge._complexity = insert_edge._ttf.size();

            const TTF& f_uxv = insert_edge._ttf;
            if ( interval_uv.get_upper() + 0.01 < f_uxv.get_min() )
                return NOT_NECESSARY;

            _context.clear_pq();

            //
            // forward single departure search in thinned predecessor graph of backward interval search
            //
            double t_start = TTF::PERIOD() / 2.0;
            sample_search(u_it, v_it, t_start);

            if ( ! _context.reached(v_it) ) return NECESSARY;

            double t_arr_ux = edge_ux.get_ttf().eval(t_start) + t_start;
            double t_arr_uxv = edge_xv.get_ttf().eval(t_arr_ux) + t_arr_ux;

            const SearchNode& v = _context.get_search_node(v_it);

            if ( ge( v._t_arr, t_arr_uxv) )
                return NECESSARY;

            _context.clear_pq();

            //
            // forward profile search in thinned predecessor graph of backward interval search
            //
            profile_search(u_it, v_it);

            const TTF& f_uv = v._ttf;

            if ( f_uv == TTF::INFTY ) return NECESSARY;
            if ( f_uxv.get_min() > f_uv.get_max() + 0.01 )
                return NOT_NECESSARY;

            UndercutDescriptor descr;
            merge(f_uv.add(0.01), f_uxv, descr);

            if ( descr._g_undercuts_strictly )
                return NECESSARY;

            return NOT_NECESSARY;
        }
    };













}

#endif /* KATCH_WITNESS_SEARCH_H_ */
