//
// Created by Bojie Shen on 22/5/21.
//

#ifndef KATCH_SEARCH_CONTEXT_H
#define KATCH_SEARCH_CONTEXT_H


#include <array>
#include <cassert>
#include <cmath>
#include <limits>
#include <vector>
#include <bitset>

#include "util/id_queue.h"


namespace katch {
    template<typename SearchNode_tmpl, typename SearchNodeId>
    class SearchContext {

    public:
        static constexpr SearchNodeId INVALID_SEARCH_NODE_ID() { return std::numeric_limits<SearchNodeId>::max(); };

    private:
        MinIDQueue _heap;
        //_table is similar to the times stamps
        std::vector<SearchNodeId> _table;
        std::vector<SearchNode_tmpl> _search_nodes;
    public:
        using SearchNode = SearchNode_tmpl;

        SearchContext(const size_t n_nodes)
                : _heap(n_nodes),
                  _table(n_nodes, INVALID_SEARCH_NODE_ID()),
                  _search_nodes() {}

        bool pq_empty() const {
            return _heap.empty();
        }

        bool pq_contains(const NodeIterator &node_iterator) const {
            SearchNodeId search_node_id = _table[node_iterator];
            if (search_node_id == INVALID_SEARCH_NODE_ID()) return false;

            return _search_nodes[search_node_id]._enqueued;
        }

        bool pq_contains(const SearchNode &search_node) const {
            return search_node._enqueued;
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

        const SearchNode &get_min() const {
            assert(!pq_empty());
            return _search_nodes[_heap.peek().id];
        }

        const double &get_min_priority() const {
            assert(!pq_empty());
            return _heap.peek().key;
        }

        SearchNode &insert(const NodeIterator &node_iterator, const double &priority) {
            assert(!reached(node_iterator));
            assert(!pq_contains(node_iterator));

            const SearchNodeId search_node_id(_search_nodes.size());

            _search_nodes.emplace_back();
            _heap.push({search_node_id, priority});
            _table[node_iterator] = search_node_id;
            _search_nodes.back()._node_it = node_iterator;
            _search_nodes.back()._enqueued = true;

            return _search_nodes.back();
        }

        void pq_re_insert(SearchNode &search_node, const double &priority) {
            assert(reached(search_node._node_it));
            assert(!pq_contains(search_node));

            const SearchNodeId search_node_id = get_search_node_id(search_node);
            _heap.push({search_node_id, priority});
            search_node._enqueued = true;
        }

        void pq_decrease(SearchNode &search_node, const double &priority) {
            assert(reached(search_node._node_it));
            assert(pq_contains(search_node._node_it));

            const SearchNodeId search_node_id = get_search_node_id(search_node);

            _heap.decrease_key({search_node_id, priority});
        }

        void pq_increase(SearchNode &search_node, const double &priority) {
            assert(reached(search_node._node_it));
            assert(pq_contains(search_node._node_it));

            const SearchNodeId search_node_id = get_search_node_id(search_node);

            _heap.increase_key({search_node_id, priority});
        }

        SearchNode &pq_delete_min() {
            assert(!pq_empty());
            SearchNode &search_node = _search_nodes[_heap.peek().id];
            assert(_table[search_node._node_it] != INVALID_SEARCH_NODE_ID());
            _heap.pop();

            search_node._enqueued = false;
            return search_node;
        }

        size_t pq_size() {
            return _heap.size();
        }

        void clear_pq() {
            _heap.clear();
        }

        void clear_all() {
            _heap.clear();
//
            for (const auto &search_node : _search_nodes)
                _table[search_node._node_it] = INVALID_SEARCH_NODE_ID();
            _search_nodes.clear();
//            _search_nodes = std::vector<SearchNode_tmpl>() ;
        }
    };
}
#endif //KATCH_SEARCH_CONTEXT_H
