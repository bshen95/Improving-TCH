//
// Created by Bojie Shen on 12/10/21.
//

#ifndef KATCH_TCH_UPDOWN_QUERY_H
#define KATCH_TCH_UPDOWN_QUERY_H
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

    class TchUpDownQuery {
    private:
        StaticForwardSearchGraph *_graph;

        MinIDQueue _queue;
        unsigned _search_id;
        unsigned number_of_nodes;
        std::vector<unsigned> _was_pushed;

        std::vector<double> _tentative_arr_t;
        std::vector<unsigned short> _first_move;
        std::vector<unsigned short> _true_first_move;
        std::vector<bool> _was_popped;

    public:
        using Graph = StaticForwardSearchGraph;

        TchUpDownQuery(StaticForwardSearchGraph *graph )
        : _graph(graph),
        _queue(_graph->get_n_nodes()*2),
        _was_pushed(_graph->get_n_nodes()*2),
        _was_popped(_graph->get_n_nodes()*2),
        _tentative_arr_t(_graph->get_n_nodes()*2),
        _first_move(_graph->get_n_nodes()*2),
        _true_first_move(_graph->get_n_nodes())
        {
            number_of_nodes  = _graph->get_n_nodes();
            _search_id = 0;
        }

        const std::vector<unsigned short> &get_first_move(const NodeIterator source) {
            // up successor < number_of_nodes; down successory > number_of_nodes;
            _search_id++;

            _queue.clear();
            fill(_first_move.begin(), _first_move.end(), 0xFF);
            _tentative_arr_t[source] = 0;
            _tentative_arr_t[source + number_of_nodes] = 0;
            _first_move[source] = 0XFE;
            _first_move[source + number_of_nodes] = 0XFE;
            _was_pushed[source] = _search_id;
            _was_pushed[source + number_of_nodes] = _search_id;


            // add first move
            unsigned short fm = 0;

            for ( EdgeIterator e = _graph->edges_begin(source) ; e != _graph->edges_end(source) ; ++e )
            {
                double next_cost = _graph->get_free_flow(e);
                NodeIterator next = _graph->get_other_node(e);
                if( _graph->is_directed_downward(e) ){
                    next = next + number_of_nodes;
                }
                _tentative_arr_t[next] = next_cost;
                _first_move[next] = fm;
                _was_pushed[next] = _search_id;
                _queue.push({next, next_cost});
                fm++;
            }
            if(fm >256){
                std::cout<<"CPD error: not enough first move symbols " <<std::endl;
            }


            while (!_queue.empty()) {
                auto x = _queue.pop();
//                if(_was_popped[x.id] == _search_id){
//                    std::cout<<"Node should only be popped once"<<std::endl;
//                }
//                _was_popped[x.id] = _search_id;
                bool is_downward_search = (x.id >= number_of_nodes);
                unsigned popped_node_id = is_downward_search ?  x.id -number_of_nodes : x.id;
                for ( EdgeIterator e = _graph->edges_begin(popped_node_id) ; e != _graph->edges_end(popped_node_id) ; ++e )
                {

                    if(is_downward_search && _graph->is_directed_upward(e)){
                        // prune down up successor;
                        continue;
                    }

                    double next_cost = x.key + _graph->get_free_flow(e);
                    NodeIterator next = _graph->get_other_node(e);
                    if( _graph->is_directed_downward(e) ){
                        next = next + number_of_nodes;
                    }


                    if (_was_pushed[next] == _search_id) {
                        if (next_cost < _tentative_arr_t[next]) {
                            _tentative_arr_t[next] = next_cost;
                            _queue.decrease_key({next, next_cost});
                            _first_move[next] = _first_move[x.id];
                        }
                    }else {
                        _queue.push({next, next_cost});
                        _tentative_arr_t[next] = next_cost;
                        _was_pushed[next] = _search_id;
                        _first_move[next] = _first_move[x.id];
                    }
                }
            }

            for(unsigned i = 0; i < number_of_nodes; i ++){
                if( _was_pushed[i] != _search_id && _was_pushed[i + number_of_nodes] != _search_id){
                    //non_reachable
                    _true_first_move[i] = 0xFF;
                }else{
                    double best_up = _was_pushed[i] == _search_id ? _tentative_arr_t[i] : std::numeric_limits<double>::max();
                    double best_down = _was_pushed[i + number_of_nodes] == _search_id ? _tentative_arr_t[i + number_of_nodes] : std::numeric_limits<double>::max();
                    _true_first_move[i] = best_up <= best_down ? _first_move[i] :  _first_move[i + number_of_nodes];
                }
            }
            _true_first_move[source] = 0XFE;
            return _true_first_move;
        }
    };
}
#endif //KATCH_TCH_UPDOWN_QUERY_H
