//
// Created by Bojie Shen on 22/5/21.
//

#ifndef KATCH_TCH_DOWNWARD_DIJKSTRA_QUERY_H
#define KATCH_TCH_DOWNWARD_DIJKSTRA_QUERY_H
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

    class TchDownardDijkstraQuery {

    private:
        StaticForwardSearchGraph *_graph;

        MinIDQueue _queue;
        unsigned _search_id;
        std::vector<unsigned > _was_pushed;
        std::vector<double> _tentative_arr_t;
        std::vector<unsigned short> _first_move;
        const std::vector<unsigned>*  _cpd_2_CH_ranking_mapper;
        std::vector<bool> _reachability;
    public:
        using Graph = StaticForwardSearchGraph;

        TchDownardDijkstraQuery(StaticForwardSearchGraph *graph )
                : _graph(graph),
                  _queue(_graph->get_n_nodes()),
                  _was_pushed(_graph->get_n_nodes()),
                  _tentative_arr_t(_graph->get_n_nodes()),
                  _first_move(_graph->get_n_nodes()),
                  _reachability(_graph->get_n_nodes())
        {
            _search_id = 0;
        }

        TchDownardDijkstraQuery(StaticForwardSearchGraph *graph,  std::vector<unsigned>* cpd_2_CH_ranking_mapper )
                : _graph(graph),
                  _cpd_2_CH_ranking_mapper(cpd_2_CH_ranking_mapper),
                  _queue(_graph->get_n_nodes()),
                  _was_pushed(_graph->get_n_nodes()),
                  _tentative_arr_t(_graph->get_n_nodes()),
                  _first_move(_graph->get_n_nodes()),
                  _reachability(_graph->get_n_nodes())
        {
            _search_id = 0;
        }


        double get_max_distance(const NodeIterator source){
            fill(_tentative_arr_t.begin(), _tentative_arr_t.end(), 0);
            get_first_move(source);
            double max = 0;
            for(auto t : _tentative_arr_t){
                if(max < t){
                    max = t;
                }
            }
            return max;
        }
        double get_distance(const NodeIterator source,const NodeIterator target){
            get_first_move(source);
            return _tentative_arr_t[target];
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
                if(_graph->is_directed_upward(e)){
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
                    if( _graph->is_directed_upward(e)){
                        // prune up successor;
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


        const std::vector<unsigned short>& get_reachability_with_wild_card(const NodeIterator source) {
            _search_id++;
            _queue.clear();
            fill( _first_move.begin(),  _first_move.end(), 0);


            _tentative_arr_t[source] = 0;
            _first_move[source] = 2;
            _was_pushed[source] = _search_id;

            for ( EdgeIterator e = _graph->edges_begin(source) ; e != _graph->edges_end(source) ; ++e )
            {
                if(_graph->is_directed_upward(e)){
                    continue;
                }
                double next_cost = _graph->get_free_flow(e);
                NodeIterator next = _graph->get_other_node(e);
                _tentative_arr_t[next] = next_cost;
                _first_move[next] = 1;
                _was_pushed[next] = _search_id;
                _queue.push({next, next_cost});
            }

            while (!_queue.empty()) {
                auto x = _queue.pop();
                for ( EdgeIterator e = _graph->edges_begin(x.id) ; e != _graph->edges_end(x.id) ; ++e )
                {
                    if( _graph->is_directed_upward(e)){
                        // prune up successor;
                        continue;
                    }
                    double next_cost = x.key + _graph->get_free_flow(e);
                    NodeIterator next = _graph->get_other_node(e);
                    if (_was_pushed[next] == _search_id) {
                        // reached
                        if (next_cost <=_tentative_arr_t[next]) {
                            _tentative_arr_t[next] = next_cost;
                            _queue.decrease_key({next, next_cost});
                            _first_move[next] =  _first_move[x.id];
                        }
                    } else {
                        _queue.push({next, next_cost});
                        _tentative_arr_t[next] = next_cost;
                        _was_pushed[next] = _search_id;
                        _first_move[next] = _first_move[x.id];
                    }
                }
            }
            unsigned source_ranking = (*_cpd_2_CH_ranking_mapper)[source];

            for ( unsigned i  = 0 ; i < _cpd_2_CH_ranking_mapper->size(); i ++){
                if((*_cpd_2_CH_ranking_mapper)[i] >= source_ranking){
                    _first_move[i] = 2;
                }

            }
            return  _first_move;
        }


        const std::vector<unsigned short> &get_first_move(const NodeIterator source) {
            _search_id++;

            _queue.clear();
            fill(_first_move.begin(), _first_move.end(), 0xFF);


            _tentative_arr_t[source] = 0;
            _first_move[source] = 0XFE;
            _was_pushed[source] = _search_id;


            // add first move
            unsigned short fm = 0;

            for ( EdgeIterator e = _graph->edges_begin(source) ; e != _graph->edges_end(source) ; ++e )
            {
                if(_graph->is_directed_upward(e)){
                    continue;
                }
                double next_cost = _graph->get_free_flow(e);
                NodeIterator next = _graph->get_other_node(e);
                _tentative_arr_t[next] = next_cost;
                _first_move[next] = fm;
                _was_pushed[next] = _search_id;
                _queue.push({next, next_cost});
            }
            if(fm >256){
                std::cout<<"CPD error: not enough first move symbols " <<std::endl;
            }

            while (!_queue.empty()) {
                auto x = _queue.pop();
                for ( EdgeIterator e = _graph->edges_begin(x.id) ; e != _graph->edges_end(x.id) ; ++e )
                {
                    if( _graph->is_directed_upward(e)){
                        // prune up successor;
                        continue;
                    }
                    double next_cost = x.key + _graph->get_free_flow(e);
                    NodeIterator next = _graph->get_other_node(e);
                    if (_was_pushed[next] == _search_id) {
                        // reached
                        if (next_cost <=_tentative_arr_t[next]) {
                            _tentative_arr_t[next] = next_cost;
                            _queue.decrease_key({next, next_cost});
                            _first_move[next] = _first_move[x.id];
                        }
                    } else {
                        _queue.push({next, next_cost});
                        _tentative_arr_t[next] = next_cost;
                        _was_pushed[next] = _search_id;
                        _first_move[next] = _first_move[x.id];
                    }
                }
            }
            return _first_move;
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
#endif //KATCH_TCH_DOWNWARD_DIJKSTRA_QUERY_H
