//
// Created by Bojie Shen on 24/10/21.
//

#ifndef KATCH_TCH_UPWARD_DIJKSTRA_QUERY_H
#define KATCH_TCH_UPWARD_DIJKSTRA_QUERY_H
#include <array>
#include <cassert>
#include <cmath>
#include <deque>
#include <limits>
#include <stack>
#include <tuple>
#include <utility>
#include <vector>


#include "datastr/base/double.h"
#include "datastr/base/interval.h"
#include "datastr/graph/basic.h"
#include "datastr/graph/static_forward_search_graph.h"
#include "util/id_queue.h"
#include <fstream>

namespace katch
{

    class TchUpwardDijkstraQuery {

    private:
        StaticForwardSearchGraph *_graph;

        MinIDQueue _queue;
        unsigned _search_id;
        std::vector<unsigned > _was_pushed;
        std::vector<double> _tentative_arr_t;
        std::vector<bool> _reachability;
    public:
        using Graph = StaticForwardSearchGraph;

        TchUpwardDijkstraQuery(StaticForwardSearchGraph *graph )
                : _graph(graph),
                  _queue(_graph->get_n_nodes()),
                  _was_pushed(_graph->get_n_nodes()),
                  _tentative_arr_t(_graph->get_n_nodes()),
                  _reachability(_graph->get_n_nodes())
        {
            _search_id = 0;
        }

        const std::vector<bool> &get_reachability(const NodeIterator source) {
            _search_id++;
            _queue.clear();
            fill( _reachability.begin(),  _reachability.end(), false);


            _tentative_arr_t[source] = 0;
            _reachability[source] = true;
            _was_pushed[source] = _search_id;

            for ( EdgeIterator e = _graph->edges_begin(source) ; e != _graph->edges_end(source) ; ++e )
            {
                if(_graph->is_directed_downward(e)){
                    continue;
                }
                double next_cost = _graph->get_free_flow(e);
                NodeIterator next = _graph->get_other_node(e);
                _tentative_arr_t[next] = next_cost;
                _reachability[next] = true;
                _was_pushed[next] = _search_id;
                _queue.push({next, next_cost});
            }

            while (!_queue.empty()) {
                auto x = _queue.pop();
                for ( EdgeIterator e = _graph->edges_begin(x.id) ; e != _graph->edges_end(x.id) ; ++e )
                {
                    if( _graph->is_directed_downward(e)){
                        // prune down successor;
                        continue;
                    }
                    double next_cost = x.key + _graph->get_free_flow(e);
                    NodeIterator next = _graph->get_other_node(e);
                    if (_was_pushed[next] == _search_id) {
                        // reached
                        if (next_cost <=_tentative_arr_t[next]) {
                            _tentative_arr_t[next] = next_cost;
                            _queue.decrease_key({next, next_cost});
                            _reachability[next] =  _reachability[x.id];
                        }
                    } else {
                        _queue.push({next, next_cost});
                        _tentative_arr_t[next] = next_cost;
                        _was_pushed[next] = _search_id;
                        _reachability[next] = _reachability[x.id];
                    }
                }
            }
            return  _reachability;
        }

        unsigned get_number_of_reachable_nodes(){
            unsigned reachable =0 ;
            for(int i = 0; i < _was_pushed.size(); i++){
                if(_was_pushed[i] == _search_id){
                    reachable ++;
                }
            }
            return reachable;
        }
    };
}
#endif //KATCH_TCH_UPWARD_DIJKSTRA_QUERY_H
