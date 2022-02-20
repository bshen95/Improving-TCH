//
// Created by Bojie Shen on 22/10/21.
//

#ifndef KATCH_REVERSE_TCH_UPDOWN_QUERY_H
#define KATCH_REVERSE_TCH_UPDOWN_QUERY_H
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

    class ReverseTchUpDownQuery {
    private:
        StaticForwardSearchGraph *_in_graph;

        StaticForwardSearchGraph *_out_graph;

        MinIDQueue _queue;
        unsigned _search_id;
        unsigned number_of_nodes;
        std::vector<unsigned> _was_pushed;

        std::vector<double> _tentative_arr_t;
        std::vector<unsigned> _predecessor_arc;
        std::vector<unsigned> _true_predecessor_arc;
        std::vector<unsigned char> _first_move;
        std::vector<bool> _was_popped;

    public:
        using Graph = StaticForwardSearchGraph;
        
        ReverseTchUpDownQuery(StaticForwardSearchGraph *in_graph )
                : _in_graph(in_graph),
                  _queue(_in_graph->get_n_nodes()*2),
                  _was_pushed(_in_graph->get_n_nodes()*2),
                  _was_popped(_in_graph->get_n_nodes()*2),
                  _tentative_arr_t(_in_graph->get_n_nodes()*2),
                  _predecessor_arc(_in_graph->get_n_nodes()*2),
                  _true_predecessor_arc(_in_graph->get_n_nodes()),
                  _first_move(_in_graph->get_n_nodes())
        {
            number_of_nodes  = _in_graph->get_n_nodes();
            _search_id = 0;
        }


        const std::vector<double> get_distance(const NodeIterator source) {
            get_first_move(source);
            vector<double> distance = vector<double>(_first_move.size());
            for(unsigned i = 0; i < number_of_nodes; i ++){
                if(i  == source){
                    continue;
                }
                if( _was_pushed[i] != _search_id && _was_pushed[i + number_of_nodes] != _search_id){
                    //non_reachable
                    _true_predecessor_arc[i] = 999999999;
                }else{
                    double best_up = _was_pushed[i] == _search_id ? _tentative_arr_t[i] : std::numeric_limits<double>::max();
                    double best_down = _was_pushed[i + number_of_nodes] == _search_id ? _tentative_arr_t[i + number_of_nodes] : std::numeric_limits<double>::max();
                    _true_predecessor_arc[i] = best_up <= best_down ? _predecessor_arc[i] :  _predecessor_arc[i + number_of_nodes];
                }
            }


        }


        const std::vector<unsigned char> &get_first_move(const NodeIterator source) {
            // up successor < number_of_nodes; down successory > number_of_nodes;
            _search_id++;

            _queue.clear();
            fill(_predecessor_arc.begin(), _predecessor_arc.end(), 999999999);
            _tentative_arr_t[source] = 0;
            _tentative_arr_t[source + number_of_nodes] = 0;
            _predecessor_arc[source] = 999999999;
            _predecessor_arc[source + number_of_nodes] = 999999999;
            _was_pushed[source] = _search_id;
            _was_pushed[source + number_of_nodes] = _search_id;



            for ( EdgeIterator e = _in_graph->edges_begin(source) ; e != _in_graph->edges_end(source) ; ++e )
            {
                double next_cost = _in_graph->get_free_flow(e);
                NodeIterator next = _in_graph->get_other_node(e);
                if( _in_graph->is_directed_downward(e) ){
                    next = next + number_of_nodes;
                }
                _tentative_arr_t[next] = next_cost;
                _predecessor_arc[next] = e;
                _was_pushed[next] = _search_id;
                _queue.push({next, next_cost});
            }


            while (!_queue.empty()) {
                auto x = _queue.pop();
//                if(_was_popped[x.id] == _search_id){
//                    std::cout<<"Node should only be popped once"<<std::endl;
//                }
//                _was_popped[x.id] = _search_id;
                bool is_downward_search = (x.id >= number_of_nodes);
                unsigned popped_node_id = is_downward_search ?  x.id -number_of_nodes : x.id;
                for ( EdgeIterator e = _in_graph->edges_begin(popped_node_id) ; e != _in_graph->edges_end(popped_node_id) ; ++e )
                {

                    if(is_downward_search && _in_graph->is_directed_upward(e)){
                        // prune down up successor;
                        continue;
                    }

                    double next_cost = x.key + _in_graph->get_free_flow(e);
                    NodeIterator next = _in_graph->get_other_node(e);
                    if( _in_graph->is_directed_downward(e) ){
                        next = next + number_of_nodes;
                    }


                    if (_was_pushed[next] == _search_id) {
                        if (next_cost < _tentative_arr_t[next]) {
                            _tentative_arr_t[next] = next_cost;
                            _queue.decrease_key({next, next_cost});
                            _predecessor_arc[next] =e;
                        }
                    }else {
                        _queue.push({next, next_cost});
                        _tentative_arr_t[next] = next_cost;
                        _was_pushed[next] = _search_id;
                        _predecessor_arc[next] = e;
                    }
                }
            }

            for(unsigned i = 0; i < number_of_nodes; i ++){
                if(i  == source){
                    continue;
                }
                if( _was_pushed[i] != _search_id && _was_pushed[i + number_of_nodes] != _search_id){
                    //non_reachable
                    _true_predecessor_arc[i] = 999999999;
                }else{
                    double best_up = _was_pushed[i] == _search_id ? _tentative_arr_t[i] : std::numeric_limits<double>::max();
                    double best_down = _was_pushed[i + number_of_nodes] == _search_id ? _tentative_arr_t[i + number_of_nodes] : std::numeric_limits<double>::max();
                    _true_predecessor_arc[i] = best_up <= best_down ? _predecessor_arc[i] :  _predecessor_arc[i + number_of_nodes];
                }
            }
            _true_predecessor_arc[source] = 0XFE;

            for(unsigned i = 0 ; i < number_of_nodes; i ++){
                if(i  == source){
                    _first_move[i] = 0xFE;
                }else if(_true_predecessor_arc[i] ==999999999){
                    _first_move[i] = 0xFF;
                }else{
                    // found first move;
                    _first_move[i] =  _in_graph->get_first_move((EdgeIterator)_true_predecessor_arc[i]);
                }
            }
            return _first_move;
        }
    };
}
#endif //KATCH_REVERSE_TCH_UPDOWN_QUERY_H
