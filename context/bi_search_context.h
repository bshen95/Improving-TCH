//
// Created by Bojie Shen on 22/5/21.
//

#ifndef KATCH_BI_SEARCH_CONTEXT_H
#define KATCH_BI_SEARCH_CONTEXT_H
#include <array>
#include <cassert>
#include <cmath>
#include <limits>
#include <vector>
#include <bitset>

#include "util/id_queue.h"
#include "datastr/base/interval.h"


namespace katch {
    template<typename SearchNode_tmpl,typename SearchNodeId>
    class BiSearchContext {
    public:
        static constexpr size_t FORWARD = 0;
        static constexpr size_t BACKWARD = 1;

    private:

        static constexpr SearchNodeId INVALID_SEARCH_NODE_ID() { return std::numeric_limits<uint32_t>::max(); };
        std::array<MinIDQueue, 2> _heap;
        //_table is similar to the times stamps
        std::vector <SearchNodeId> _table;
        std::vector <SearchNode_tmpl> _search_nodes;

    public:

        using SearchNode = SearchNode_tmpl;

        BiSearchContext(const size_t n_nodes)
                :
                _heap(),
                _table(n_nodes, INVALID_SEARCH_NODE_ID()),
                _search_nodes() {
            _heap[0] = MinIDQueue(n_nodes);
            _heap[1] = MinIDQueue(n_nodes);
        }

        bool pq_empty(const size_t direction) const {
            assert(direction < _heap.size());
            return _heap[direction].empty();
        }

        bool pq_contains(const NodeIterator &node_iterator, const size_t direction) const {
            SearchNodeId search_node_id = _table[node_iterator];
            if (search_node_id == INVALID_SEARCH_NODE_ID()) return false;

            return _search_nodes[search_node_id]._enqueued[direction];
        }

        bool pq_contains(const SearchNode &search_node, const size_t direction) const {
            return search_node._enqueued[direction];
        }

        bool reached(const NodeIterator &node_iterator) const {
            return _table[node_iterator] != INVALID_SEARCH_NODE_ID();
        }

        const SearchNode &get_search_node(const NodeIterator &node_iterator) const {
            assert(_table[node_iterator] != INVALID_SEARCH_NODE_ID());
            return _search_nodes[_table[node_iterator]];
        }

        SearchNode &get_search_node(const NodeIterator &node_iterator) {
            assert(_table[node_iterator] != INVALID_SEARCH_NODE_ID());
            return _search_nodes[_table[node_iterator]];
        }

        SearchNode &get_search_node_from_id(const SearchNodeId search_node_id) {
            assert(search_node_id < _search_nodes.size());
            return _search_nodes[search_node_id];
        }

        const SearchNode &get_search_node_from_id(const SearchNodeId search_node_id) const {
            assert(search_node_id < _search_nodes.size());
            return _search_nodes[search_node_id];
        }

        SearchNodeId get_search_node_id(const SearchNode &search_node) const {
            assert(reached(search_node._node_it));

            const SearchNodeId search_node_id = std::distance(&(*(_search_nodes.begin())), &search_node);
            return search_node_id;
        }

        const SearchNode &get_min(const size_t direction) const {
            assert(!pq_empty(direction));
            return _search_nodes[_heap[direction].peek().id];
        }

        const double get_min_priority(const size_t direction) const {
            assert(!pq_empty(direction));
            return _heap[direction].peek().key;
        }

        SearchNode &insert(const NodeIterator &node_iterator, const double &priority, const size_t direction) {
            assert(!reached(node_iterator));
            assert(!pq_contains(node_iterator, direction));

            const SearchNodeId search_node_id(_search_nodes.size());

            _search_nodes.emplace_back();
            _heap[direction].push({search_node_id, priority});
            _table[node_iterator] = search_node_id;
            _search_nodes.back()._node_it = node_iterator;
            _search_nodes.back()._enqueued[direction] = true;

            return _search_nodes.back();
        }

        void pq_re_insert(SearchNode &search_node, const double &priority, const size_t direction) {
            assert(reached(search_node._node_it));
            assert(!pq_contains(search_node, direction));

            const SearchNodeId search_node_id = get_search_node_id(search_node);
            _heap[direction].push({search_node_id, priority});
            search_node._enqueued[direction] = true;
        }

        void pq_decrease(SearchNode &search_node, const double &priority, const size_t direction) {
            assert(reached(search_node._node_it));
            assert(pq_contains(search_node._node_it, direction));

            const SearchNodeId search_node_id = get_search_node_id(search_node);

//                assert( (*heap_handle)._priority >= priority );
            _heap[direction].decrease_key({search_node_id, priority});
        }

        SearchNode &pq_delete_min(const size_t direction) {
            assert(!pq_empty(direction));

            SearchNode &search_node = _search_nodes[_heap[direction].peek().id];
            assert(_table[search_node._node_it] != INVALID_SEARCH_NODE_ID());
            _heap[direction].pop();

            search_node._enqueued[direction] = false;
            return search_node;
        }

        void clear_pq(const size_t direction) {
//                for ( const auto& heap_item : _heap[direction] )
//                    _search_nodes[heap_item._search_node_id]._enqueued[direction] = false;
            for (auto &search_node:_search_nodes)
                search_node._enqueued[direction] = false;
            _heap[direction].clear();
        }

        void clear_all() {
            _heap[FORWARD].clear();
            _heap[BACKWARD].clear();

            for (const auto &search_node : _search_nodes)
                _table[search_node._node_it] = INVALID_SEARCH_NODE_ID();

            _search_nodes.clear();
        }
    };
}
#endif //KATCH_BI_SEARCH_CONTEXT_H
