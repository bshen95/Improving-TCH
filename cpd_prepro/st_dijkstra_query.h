//
// Created by Bojie Shen on 12/5/21.
//

#ifndef KATCH_ST_DIJKSTRA_QUERY_H
#define KATCH_ST_DIJKSTRA_QUERY_H

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
#include "datastr/graph/static_search_graph.h"
#include "util/id_queue.h"
#include <boost/geometry.hpp>
namespace geom = boost::geometry;
namespace katch
{

    class StDijkstraQuery {

    private:
        StaticSearchGraph *_graph;

        MinIDQueue _queue;
        unsigned short _search_id;
        std::vector<unsigned short> _was_pushed;
        std::vector<double> _tentative_arr_t;
        std::vector<unsigned short> _first_move;

    public:
        using Graph = StaticSearchGraph;

        StDijkstraQuery(StaticSearchGraph *graph )
                : _graph(graph),
                  _queue(_graph->get_n_nodes()),
                  _was_pushed(_graph->get_n_nodes()),
                  _tentative_arr_t(_graph->get_n_nodes()),
                  _first_move(_graph->get_n_nodes())
                  {
            _search_id = 0;
        }
        unsigned get_number_of_non_reachable(const NodeIterator source){
            get_first_move(source);
            unsigned int num =0;
            std::ofstream myFile("dataset/NY_wkt.csv");
            myFile<<"road_id;wkt\n";
            int index = 0 ;
            for(int i = 0; i < _first_move.size() ; i ++ ){
                if(_first_move[i] != 0XFF){
                    myFile<<std::fixed<<std::setprecision(8)<<index<<";"<<boost::geometry::wkt(_graph->_coordinate[i]) <<"\n";
                    index++;
                }else{
                    num ++;
                }
            }
            myFile.close();
            return num ;
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


        std::vector<double> get_distance_to_all_nodes(const NodeIterator source){
            _search_id++;
            _queue.clear();
            fill(_tentative_arr_t.begin(), _tentative_arr_t.end(), std::numeric_limits<double>::max());
            _tentative_arr_t[source] = 0;
            _was_pushed[source] = _search_id;
            _queue.push({source, 0});

            while (!_queue.empty()) {
                auto x = _queue.pop();
                for ( EdgeIterator e = _graph->edges_begin(x.id) ; e != _graph->edges_end(x.id) ; ++e )
                {
                    double next_cost = x.key + _graph->get_free_flow(e);
                    NodeIterator next = _graph->get_other_node(e);
                    if (_was_pushed[next] == _search_id) {
                        // reached
                        if (lt(next_cost, _tentative_arr_t[next])) {
                            _tentative_arr_t[next] = next_cost;
                            _queue.decrease_key({next, next_cost});
                        }
                    } else {
                        _queue.push({next, next_cost});
                        _tentative_arr_t[next] = next_cost;
                        _was_pushed[next] = _search_id;
                    }

                }

            }
            return _tentative_arr_t;
        }


        const std::vector<unsigned short> &get_first_move(const NodeIterator source) {
//            fill(_tentative_arr_t.begin(), _tentative_arr_t.end(), std::numeric_limits<double>::max());
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
                _queue.push({next, next_cost});
                fm++;
            }


            while (!_queue.empty()) {
                auto x = _queue.pop();
                for ( EdgeIterator e = _graph->edges_begin(x.id) ; e != _graph->edges_end(x.id) ; ++e )
                {
                    double next_cost = x.key + _graph->get_free_flow(e);
                    NodeIterator next = _graph->get_other_node(e);
                    if (_was_pushed[next] == _search_id) {
                        // reached
                        if (lt(next_cost, _tentative_arr_t[next])) {
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
    };
}

#endif //KATCH_ST_DIJKSTRA_QUERY_H
