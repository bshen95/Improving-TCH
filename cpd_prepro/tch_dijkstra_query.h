//
// Created by Bojie Shen on 24/5/21.
//

#ifndef KATCH_TCH_DIJKSTRA_QUERY_H
#define KATCH_TCH_DIJKSTRA_QUERY_H
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

    class TchDijkstraQuery {

    private:
        StaticForwardSearchGraph *_graph;

        MinIDQueue _queue;
        unsigned _search_id;
        std::vector<unsigned> _was_pushed;
        std::vector<unsigned> _predecessor;
        std::vector<unsigned> _path_length;
        std::vector<double> _tentative_arr_t;
        std::vector<unsigned short> _first_move;
        std::vector<bool> _is_downward_search;
        std::vector<bool> _is_upward_search;
        const std::vector<unsigned>*  _cpd_2_CH_ranking_mapper;
    public:
        using Graph = StaticForwardSearchGraph;

        TchDijkstraQuery(StaticForwardSearchGraph *graph )
                : _graph(graph),
                  _queue(_graph->get_n_nodes()),
                  _was_pushed(_graph->get_n_nodes()),
                  _tentative_arr_t(_graph->get_n_nodes()),
                  _first_move(_graph->get_n_nodes()),
                  _is_downward_search(_graph->get_n_nodes()),
                  _is_upward_search(_graph->get_n_nodes()),
                  _predecessor(_graph->get_n_nodes()),
                  _path_length(_graph->get_n_nodes())
        {
            _search_id = 0;
        }



        TchDijkstraQuery(StaticForwardSearchGraph *graph , std::vector<unsigned>* cpd_2_CH_ranking_mapper)
                : _graph(graph),
                  _queue(_graph->get_n_nodes()),
                  _was_pushed(_graph->get_n_nodes()),
                  _tentative_arr_t(_graph->get_n_nodes()),
                  _first_move(_graph->get_n_nodes()),
                  _is_downward_search(_graph->get_n_nodes()),
                  _cpd_2_CH_ranking_mapper(cpd_2_CH_ranking_mapper)
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
            unsigned current  = target;
            while(_predecessor[current] !=source){
                std::cout<<current<< " Dist: "<<_tentative_arr_t[current]<< std::endl;
                current = _predecessor[current];
            }

            std::cout<<_tentative_arr_t[1503] <<" "<<_predecessor[1503]<<std::endl;
            std::cout<<_tentative_arr_t[23591] <<" "<<_predecessor[23591]<<std::endl;
            return _tentative_arr_t[target];
        }
//
//        const std::vector<unsigned short> &get_first_move(const NodeIterator source) {
//            _search_id++;
//
//            _queue.clear();
//            fill(_first_move.begin(), _first_move.end(), 0xFF);
//
//
//            _tentative_arr_t[source] = 0;
//            _first_move[source] = 0XFE;
//            _was_pushed[source] = _search_id;
//
//
//            // add first move
//            unsigned short fm = 0;
//
//            for ( EdgeIterator e = _graph->edges_begin(source) ; e != _graph->edges_end(source) ; ++e )
//            {
//                double next_cost = _graph->get_free_flow(e);
//                NodeIterator next = _graph->get_other_node(e);
//                _tentative_arr_t[next] = next_cost;
//                _first_move[next] = fm;
//                _was_pushed[next] = _search_id;
//                _is_downward_search[next] = _graph->is_directed_downward(e);
//                _queue.push({next, next_cost});
//                _predecessor[next] = source;
//                fm++;
//            }
//            if(fm >256){
//                std::cout<<"CPD error: not enough first move symbols " <<std::endl;
//            }
//
//            while (!_queue.empty()) {
//                auto x = _queue.pop();
////                std::cout<<"Pop nodes "<<x.id << "cost"<< x.key<<std::endl;
//                for ( EdgeIterator e = _graph->edges_begin(x.id) ; e != _graph->edges_end(x.id) ; ++e )
//                {
//                    if(_is_downward_search[x.id] && _graph->is_directed_upward(e)){
//                        // prune down up successor;
//                        continue;
//                    }
//
//                    double next_cost = x.key + _graph->get_free_flow(e);
//                    NodeIterator next = _graph->get_other_node(e);
//                    if (_was_pushed[next] == _search_id) {
//                        // reached
//                        if (next_cost < _tentative_arr_t[next]) {
//                            _tentative_arr_t[next] = next_cost;
//                            _queue.decrease_key({next, next_cost});
//                            _first_move[next] = _first_move[x.id];
//                            _is_downward_search[next] = _graph->is_directed_downward(e);
//                            _predecessor[next] = x.id;
////                            std::cout<<"Generating node "<<next << "cost"<< next_cost<<std::endl;
//                        }else if( next_cost == _tentative_arr_t[next]){
//                            if(_graph->is_directed_upward(e)){
//                                _first_move[next] = _first_move[x.id];
//                            }
//                            _is_downward_search[next] = _is_downward_search[next] && _graph->is_directed_downward(e);
//                        }
//                    } else {
//                        _queue.push({next, next_cost});
//                        _tentative_arr_t[next] = next_cost;
//                        _was_pushed[next] = _search_id;
//                        _first_move[next] = _first_move[x.id];
//                        _is_downward_search[next] = _graph->is_directed_downward(e);
//                        _predecessor[next] = x.id;
////                        std::cout<<"Generating node "<<next << "cost"<< next_cost<<std::endl;
//                    }
//
//                }
//
//            }
//
//            return _first_move;
//        }
        const std::vector<unsigned short> &get_first_move(const NodeIterator source) {
            _search_id++;

            _queue.clear();
            fill(_first_move.begin(), _first_move.end(), 0xFF);


            _tentative_arr_t[source] = 0;
            _first_move[source] = 0XFE;
            _was_pushed[source] = _search_id;
            _path_length[source] = 0;

            // add first move
            unsigned short fm = 0;

            for ( EdgeIterator e = _graph->edges_begin(source) ; e != _graph->edges_end(source) ; ++e )
            {
                double next_cost = _graph->get_free_flow(e);
                NodeIterator next = _graph->get_other_node(e);
                _tentative_arr_t[next] = next_cost;
                _first_move[next] = fm;
                _was_pushed[next] = _search_id;
                _path_length[next] = 1;
                _queue.push({next, next_cost});
                fm++;
            }
            if(fm >256){
                std::cout<<"CPD error: not enough first move symbols " <<std::endl;
            }

            while (!_queue.empty()) {
                auto x = _queue.pop();
                for ( EdgeIterator e = _graph->edges_begin(x.id) ; e != _graph->edges_end(x.id) ; ++e )
                {
                    double next_cost = x.key + _graph->get_free_flow(e);
                    NodeIterator next = _graph->get_other_node(e);
                    if (_was_pushed[next] == _search_id) {
                        // reached
                        if (next_cost < _tentative_arr_t[next]) {
                            _tentative_arr_t[next] = next_cost;
                            _queue.decrease_key({next, next_cost});
                            _first_move[next] = _first_move[x.id];
                            _path_length[next] = _path_length[x.id] +1;
                        }else if( eq(next_cost,_tentative_arr_t[next])){
                            if(_path_length[next] > _path_length[x.id] +1){
                                _path_length[next] = _path_length[x.id] +1;
                                _first_move[next] = _first_move[x.id];
                            }
                        }
                    } else {
                        _queue.push({next, next_cost});
                        _tentative_arr_t[next] = next_cost;
                        _was_pushed[next] = _search_id;
                        _first_move[next] = _first_move[x.id];
                        _path_length[next] = _path_length[x.id] +1;
                    }
                }
            }
            return _first_move;
        }



        const std::vector<unsigned short> &get_first_move_down_only(const NodeIterator source) {
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
                double next_cost = _graph->get_free_flow(e);
                NodeIterator next = _graph->get_other_node(e);
                _tentative_arr_t[next] = next_cost;
                _first_move[next] = fm;
                _was_pushed[next] = _search_id;
                _is_downward_search[next] = _graph->is_directed_downward(e);
                _queue.push({next, next_cost});
                fm++;
            }
            if(fm >256){
                std::cout<<"CPD error: not enough first move symbols " <<std::endl;
            }

            while (!_queue.empty()) {
                auto x = _queue.pop();
                for ( EdgeIterator e = _graph->edges_begin(x.id) ; e != _graph->edges_end(x.id) ; ++e )
                {
                    if(_is_downward_search[x.id] && _graph->is_directed_upward(e)){
                        // prune down up successor;
                        continue;
                    }
                    double next_cost = x.key + _graph->get_free_flow(e);
                    NodeIterator next = _graph->get_other_node(e);
                    if (_was_pushed[next] == _search_id) {
                        // reached
                        if (next_cost < _tentative_arr_t[next]) {
                            _tentative_arr_t[next] = next_cost;
                            _queue.decrease_key({next, next_cost});
                            _first_move[next] = _first_move[x.id];
                            _is_downward_search[next] = _graph->is_directed_downward(e);
                        }else if( next_cost == _tentative_arr_t[next]){
                            if(_graph->is_directed_upward(e)){
                                _first_move[next] = _first_move[x.id];
                            }
                            _is_downward_search[next] = _is_downward_search[next] && _graph->is_directed_downward(e);
                        }
                    } else {
                        _queue.push({next, next_cost});
                        _tentative_arr_t[next] = next_cost;
                        _was_pushed[next] = _search_id;
                        _first_move[next] = _first_move[x.id];
                        _is_downward_search[next] = _graph->is_directed_downward(e);
                    }

                }

            }
            unsigned source_ranking = (*_cpd_2_CH_ranking_mapper)[source];

            for ( unsigned i  = 0 ; i < _cpd_2_CH_ranking_mapper->size(); i ++){
                if((*_cpd_2_CH_ranking_mapper)[i] >= source_ranking){
                    _first_move[i] = 0XFE;
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
#endif //KATCH_TCH_DIJKSTRA_QUERY_H
